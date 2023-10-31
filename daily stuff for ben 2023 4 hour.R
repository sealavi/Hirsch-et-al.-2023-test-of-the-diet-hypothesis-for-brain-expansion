library(lubridate)

RV1=read.csv("revisits_12min_threshold with vedba2023_4_hour.csv")
RV1=RV1[,-1]

studaydays2=unique(RV1$day)
IDs=unique(RV1$ID)
RV2=read.csv("FFT Daily Path Length 2023.csv")
RV2=RV2[,-1]




UDs=list.files("C:/Users/salavi/Documents/Shapefiles",pattern=".shp", full.names=TRUE)

crowns= readOGR("C:/Users/salavi/Downloads/documents-export-2022-03-23/BCI_Dipteryx_Master_final.shp")
crowns = spTransform(crowns, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
centroids=geosphere::centroid(crowns)
centroids=data.frame(centroids)
centroids$pkuid=crowns$pkuid

library(rgdal)
library(raster)

dip_in_range=c()
for(i in 1:46){
  id=IDs[i]
  if(length(which(as.character(RV1$ID)==id))==0){next}
  if(length(which(grepl(id, UDs, fixed = TRUE)==TRUE))==0){next}
  UD=readOGR(UDs[[which(grepl(id, UDs, fixed = TRUE)==TRUE)]])
  UD=spTransform(UD, CRS('+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs'))
  if(unique(RV1$study_period[which(as.character(RV1$ID)==id)])=="First"){
    crowns2=crowns[which(crowns$in15==1),]
  }else{
    if(unique(RV1$study_period[which(as.character(RV1$ID)==id)])=="Second"){
      crowns2=crowns[which(crowns$in17==1),]
    }
  }
  
  track_trees=crop(crowns2,UD)
  dip_crowns=length(track_trees)
  x=data.frame(cbind(id,dip_crowns))
  dip_in_range[i]=list(x)
}
dip_in_range2=data.table::rbindlist(dip_in_range)
dip_in_range2=data.frame(dip_in_range2)
colnames(dip_in_range2)=c("ID","dip_in_range")

cleaned_data=read.csv("FFT_Cleaned_Final.csv")
cleaned_data$timestamp<-as.POSIXct(as.character(cleaned_data$timestamp), format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="UTC")

cleaned_data2 <- SpatialPointsDataFrame(coords = cleaned_data[,c(3,4)], data = cleaned_data,
                                        proj4string=CRS("+proj=longlat +datum=WGS84"))
cleaned_data2 <- spTransform(cleaned_data2, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
cleaned_data2=data.frame(cleaned_data2)
cleaned_data2$Z=cleaned_data2$location.long.2+1i*cleaned_data2$location.lat.2
cleaned_data2$loopID=paste(cleaned_data2$individual.local.identifier,cleaned_data2$day)
clean_split=split(cleaned_data2,as.factor(cleaned_data2$loopID))

days_2015=seq(from=as.POSIXct("2015-12-11 10:00:00", format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="UTC"),
              to = as.POSIXct("2016-04-19 10:00:00", format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="UTC"),
              by = "24 hours")
days_20152=as.numeric(as.factor(as.character(days_2015)))

nights_2015=seq(from=as.POSIXct("2015-12-11 23:00:00", format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="UTC"),
                to = as.POSIXct("2016-04-19 23:00:00", format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="UTC"),
                by = "24 hours")
nights_20152=as.numeric(as.factor(as.character(nights_2015)))

days_2017=seq(from=as.POSIXct("2017-12-01 10:00:00", format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="UTC"),
              to = as.POSIXct("2018-06-15 10:00:00", format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="UTC"),
              by = "24 hours")
days_20172=as.numeric(as.factor(as.character(days_2017)))

nights_2017=seq(from=as.POSIXct("2017-12-01 23:00:00", format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="UTC"),
                to = as.POSIXct("2018-06-15 23:00:00", format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="UTC"),
                by = "24 hours")
nights_20172=as.numeric(as.factor(as.character(nights_2017)))

Gridded = split(cleaned_data2,as.factor(cleaned_data2$individual.local.identifier))

for(i in 1:length(Gridded)){
  Gridded[[i]]$day=NA
  if(unique(Gridded[[i]]$individual.taxon.canonical.name)=="Potos flavus"){
    if(min(unique(lubridate::year(Gridded[[i]]$timestamp)),na.rm=TRUE)<2017){
      Days=nights_2015
      Days2=nights_20152
    }else{
      Days=nights_2017
      Days2=nights_20172
    }
  }else{
    if(min(unique(lubridate::year(Gridded[[i]]$timestamp)),na.rm=TRUE)<2017){
      Days=days_2015
      Days2=days_20152
    }else{
      Days=days_2017
      Days2=days_20172
    }
  }
  for(j in 2:length(Days)){
    Gridded[[i]]$day[which(Gridded[[i]]$timestamp>=Days[j-1] & Gridded[[i]]$timestamp<Days[j])]=Days2[j-1]
  }
  
}



cleaned_data2=data.table::rbindlist(Gridded)
cleaned_data2=cleaned_data2[-which(cleaned_data2$individual.local.identifier=="Atlas 4673"),]
cleaned_data2=cleaned_data2[-which(cleaned_data2$individual.local.identifier=="Merk 4665"),]
cleaned_data2=cleaned_data2[-which(cleaned_data2$individual.local.identifier=="Serge 4670"),]
cleaned_data2=cleaned_data2[-which(cleaned_data2$individual.local.identifier=="Sahti 4693"),]
cleaned_data2=cleaned_data2[-which(cleaned_data2$individual.local.identifier=="Chloe 4052"),]


ranges=read.csv("ctmm_ranges (1).csv")
ranges$low[which(ranges$units=="area (square kilometers)")]=ranges$low[which(ranges$units=="area (square kilometers)")]*100
ranges$est[which(ranges$units=="area (square kilometers)")]=ranges$est[which(ranges$units=="area (square kilometers)")]*100
ranges$high[which(ranges$units=="area (square kilometers)")]=ranges$high[which(ranges$units=="area (square kilometers)")]*100

colnames(ranges)=c("X","ID","Home_range_low","Home_range_estimate","Home_range_high","units")
RV2$Home_range_low=NA
RV2$Home_range_estimate=NA
RV2$Home_range_high=NA


for(i in 1:nrow(ranges)){
  if(length(which(RV2$ID==ranges$ID[i]))==0){
    next
  }
  RV2$Home_range_low[which(RV2$ID==ranges$ID[i])]=ranges$Home_range_low[i]
  RV2$Home_range_estimate[which(RV2$ID==ranges$ID[i])]=ranges$Home_range_estimate[i]
  RV2$Home_range_high[which(RV2$ID==ranges$ID[i])]=ranges$Home_range_high[i]
  
}

RV1$exit_time<-as.POSIXct(as.character(RV1$exit_time), format="%Y-%m-%d %H:%M:%S",origin="01-01-1900", tz="America/Panama")
RV1$entry_time<-as.POSIXct(as.character(RV1$entry_time), format="%Y-%m-%d %H:%M:%S",origin="01-01-1900", tz="America/Panama")

hr=read.csv("FFT metadata.csv")


RV1$Weight=NA
RV1$Sex=NA
for(i in 1:nrow(hr)){
  RV1$Weight[which(RV1$ID==hr$ID[i])]=hr$weight..kg.[i]
  RV1$Sex[which(RV1$ID==hr$ID[i])]=hr$sex[i]
}

id_table=unique(RV2[,c(1,3)])
id_table$Day=as.integer(id_table$Study_day)

firstday2015=lubridate::date(min(RV1$entry_time[which(RV1$year==2015)],na.rm=TRUE))-lubridate::days(5)
lastday2015=lubridate::date(max(RV1$entry_time[which(RV1$year==2016)],na.rm=TRUE))

firstday2017=lubridate::date(min(RV1$entry_time[which(RV1$year==2017& RV1$ID!="Merk 4665")],na.rm=TRUE))-lubridate::days(4)
lastday2017=lubridate::date(max(RV1$entry_time[which(RV1$year==2018& RV1$ID!="Merk 4665")],na.rm=TRUE))

weeks2015=seq(from=firstday2015, to=lastday2015,by="1 day")
weeks2017=seq(from=firstday2017, to=lastday2017,by="1 day")

weeks2015_2=sort(rep(seq(1:20),7))[1:length(weeks2015)]
weeks2017_2=sort(rep(seq(1:28),7))[1:length(weeks2017)]

weeks=cbind(c(as.character(weeks2015),as.character(weeks2017)),c(weeks2015_2,weeks2017_2))
weeks = data.frame(weeks)
weeks$X1=lubridate::as_datetime(weeks$X1)

RV2$Number_Dip_visits=NA
RV2$Number_Dip_visits_above_12=NA
RV2$active_minutes=NA
RV2$Weight=NA
RV2$Sex=NA
RV2$Year=NA
RV2$study_period=NA
RV2$week=NA
RV2$Week_year=NA

for(i in 1:nrow(id_table)){
  temp=RV1[which((as.character(RV1$ID)==id_table$ID[i])&(RV1$day==id_table$Day[i])),]
  if(nrow(temp)==0){next}
  visits=nrow(temp)
  visits_above_12=nrow(temp[which(temp$visit_duration>=720),])
  if(unique(temp$Species)=="Potos flavus"){
    temp$vedba_above_35[which(is.nan(temp$vedba_above_35)==TRUE)]=0.895 
  } else{
    if(unique(temp$Species)=="Nasua narica"){
      temp$vedba_above_35[which(is.nan(temp$vedba_above_35)==TRUE)]=0.846  
    } else{
      if(unique(temp$Species)=="Cebus capucinus"){
        temp$vedba_above_35[which(is.nan(temp$vedba_above_35)==TRUE)]=0.834   
      } else{
        if(unique(temp$Species)=="Ateles geoffroyi"){
          temp$vedba_above_35[which(is.nan(temp$vedba_above_35)==TRUE)]=0.588   
        }
      }
    }
  }
  active_minutes=(sum(temp$visit_duration*temp$vedba_above_35))/60
  Weight=unique(temp$Weight)
  Sex=unique(temp$Sex)
  Year=unique(temp$year)
  week=unique(temp$week)
  Study_day=unique(temp$day)
  study_period = unique(temp$study_period)
  Week_year=unique(paste(temp$week,temp$year, sep = "_"))
  
  RV2$Number_Dip_visits[which((RV2$ID==id_table$ID[i])&(RV2$Study_day==id_table$Study_day[i]))]=visits
  RV2$Number_Dip_visits_above_12[which((RV2$ID==id_table$ID[i])&(RV2$Study_day==id_table$Study_day[i]))]=visits_above_12
  RV2$active_minutes[which((RV2$ID==id_table$ID[i])&(RV2$Study_day==id_table$Study_day[i]))]=active_minutes
  
  RV2$Weight[which((RV2$ID==id_table$ID[i])&(RV2$Study_day==id_table$Study_day[i]))]=Weight
  RV2$Sex[which((RV2$ID==id_table$ID[i])&(RV2$Study_day==id_table$Study_day[i]))]=Sex
  RV2$Year[which((RV2$ID==id_table$ID[i])&(RV2$Study_day==id_table$Study_day[i]))]=Year
  RV2$week[which((RV2$ID==id_table$ID[i])&(RV2$Study_day==id_table$Study_day[i]))]=week
  RV2$Study_day[which((RV2$ID==id_table$ID[i])&(RV2$Study_day==id_table$Study_day[i]))]=Study_day
  RV2$Week_year[which((RV2$ID==id_table$ID[i])&(RV2$Study_day==id_table$Study_day[i]))]=Week_year
  RV2$study_period[which((RV2$ID==id_table$ID[i])&(RV2$Study_day==id_table$Study_day[i]))]=study_period
}



for(i in 1:nrow(hr)){
  RV2$Weight[which(RV2$ID==hr$ID[i] & is.na(RV2$Weight))]=hr$weight..kg.[i]
  RV2$Sex[which(RV2$ID==hr$ID[i] & is.na(RV2$Sex))]=hr$sex[i]
  RV2$Year[which(RV2$ID==hr$ID[i] & is.na(RV2$Year))]=hr$Year[i]
}

RV2=RV2[
  with(RV2, order(ID, Study_day)),
]

RV2$date=lubridate::as_datetime(RV2$date)

for(i in 1:nrow(weeks)){
  RV2$week[which(RV2$date==weeks[i,1])]=weeks[i,2]
  
}
RV2$Week_year=paste(RV2$week,RV2$Year, sep = "_")

RV2=RV2[complete.cases(RV2$Interpolated_DJL),]

RV2$study_period[which(RV2$Year<2017)]="First"
RV2$study_period[which(RV2$Year>2016)]="Second"

RV1=RV1[
  with(RV1, order(ID, day)),
]



write.csv(RV1,file = "revisits_12min_threshold with vedba2023_4_hour.csv")

library(dplyr)
RV1b=RV1 %>%
  group_by(ID, day) %>% 
  summarise(across(matches(c("visit_duration")), list(sum = ~sum(.x, na.rm = TRUE))))
RV1b=RV1b[,-4]


RV2_sum <- merge(RV2, RV1b, by.x = c("ID","Study_day"), by.y = c("ID","day"), all.x = TRUE)

RV2_sum$dip_in_range=NA
for(i in 1:nrow(dip_in_range2)){
  RV2_sum$dip_in_range[which(as.character(RV2_sum$ID)==dip_in_range2$ID[i])]=dip_in_range2$dip_in_range[i]
}

write.csv(RV2_sum,file = "Daily_stats_for_Ben_2023_4_hour_updated.csv")

