library(sp)
library(rgdal)
library(lubridate)
library(ctmm)
library(sp)

load("~/onetofour.RData")
load("~/fivetoeight.RData")
load("~/ninetotwelve.RData")
load("~/thirteentosexteen.RData")
load("~/seventeentoeighteen.RData")
load("~/nonteentotwentytwo.RData")
load("~/twentythree.RData")
load("~/twentyfour.RData")
load("~/twentyfive.RData")
load("~/twentysix.RData")
load("~/twentyseventothirty.RData")
load("~/thirtyonetothirtyfour.RData")
load("~/thirtyninetofourtytwo.RData")
load("~/fourty.RData")
load("~/fourtyone.RData")
load("~/fourtytwo.RData")
load("~/fourtythreetofourtyfour.RData")

thirtyninetofourtytwo=thirtyninetofourtytwo[-which(thirtyninetofourtytwo$ID == "Tony Stark 4659"),]

compileddata = list(onetofour,fivetoeight,ninetotwelve,
                    thirteentosexteen,seventeentoeighteen,
                    nonteentotwentytwo,twentythree, twentyfour, twentyfive, 
                    twentysix,twentyseventothirty,
                    thirtyonetothirtyfour,thirtyninetofourtytwo,fourty, fourtyone,
                    fourtytwo, fourtythreetofourtyfour
)

gc()

rm(list = c("onetofour","fivetoeight","ninetotwelve",
            "thirteentosexteen","seventeentoeighteen",
            "nonteentotwentytwo","twentythree", "twentyfour", "twentyfive", 
            "twentysix","twentyseventothirty",
            "thirtyonetothirtyfour","thirtyninetofourtytwo","fourty", "fourtyone",
            "fourtytwo", "fourtythreetofourtyfour"
))

gc()
compileddata = do.call(rbind, compileddata)

compileddata$loopID = paste(compileddata$ID, compileddata$Study_day, sep = "_")

compileddata$loopID = as.factor(compileddata$loopID)

compileddata$timestamp = lubridate::with_tz(compileddata$timestamp, tzone = "America/Panama")
gc()

datasplit=split(compileddata,compileddata$loopID)
for(i in 1:length(datasplit)){
  mintime = datasplit[[i]]$timestamp[1]+lubridate::hours(4)
  datasplit[[i]]=datasplit[[i]][which(datasplit[[i]]$timestamp<=mintime),]
}
compileddata=data.table::rbindlist(datasplit)
rm(datasplit)
gc()

compileddata <- SpatialPointsDataFrame(coords = compileddata[,c(4,5)], data = compileddata,
                                proj4string=CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
#data2 <- spTransform(data2, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))



crowns= readOGR("C:/Users/salavi/Downloads/documents-export-2022-03-23/BCI_Dipteryx_Master_final.shp")
crowns = spTransform(crowns, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
centroids=geosphere::centroid(crowns)
centroids=data.frame(centroids)
centroids$pkuid=crowns$pkuid
dat=over(compileddata,crowns)
gc()
dat=data.frame(dat)
gc()
dat2=cbind(data.frame(compileddata),dat)
rm(compileddata)
rm(dat)

gc()
dat2$visit=NA

#datasplit=split(dat2,dat2$loopID)
gc()

datasplit2=split(dat2,as.factor(dat2$ID))
rm(dat2)
gc()

count=1
visitdata=c()
for(i in 1: length(datasplit2)){
  print(i)
  #print(gc())
  temp=data.frame(datasplit2[[i]])
  if(max(unique(lubridate::year(temp$timestamp)), na.rm=TRUE)<2017){
    temp$visit[which(is.na(temp$PatchID) | temp$in15==0)]="Out"
    temp$visit[which(!is.na(temp$PatchID) & temp$in15==1)]="In"
  } else{
    temp$visit[which(is.na(temp$PatchID) | temp$in17==0)]="Out"
    temp$visit[which(!is.na(temp$PatchID) & temp$in17==1)]="In"
  }

  consec=rle(temp$visit)
  consec=data.frame(cbind(consec[[1]],consec[[2]]))
  consec$id=unique(temp$ID)
  consec$species=unique(temp$Species)
  consec2=consec[which(consec[,2]=="Out"),]
  consec2[,1]=as.numeric(consec2[,1])
  
  bouts=rle(temp$visit)
  bouts=data.frame(cbind(bouts[[1]],bouts[[2]]))
  bouts$id=unique(temp$ID)
  bouts$species=unique(temp$Species)
  bouts2=bouts[which(bouts[,2]=="In"),]
  bouts2[,1]=as.numeric(bouts2[,1])
  consec2$index=rownames(consec2)
  bouts2$index=rownames(bouts2)
  revisits=rbind(bouts2,consec2)
  revisits=data.frame(revisits)
  revisits$index=as.numeric(revisits$index)
  revisits=revisits[with(revisits, order(index)),]
  revisits$index2=NA
  revisits$index3=NA
  revisits$index2[1]=1
  revisits$index3[1]=revisits$X1[1]
  revisits2=c()
  for(j in 2:nrow(revisits)){
    revisits$index2[j]= revisits$index3[j-1]+1
    revisits$index3[j]= revisits$index2[j]+(revisits$X1[j]-1)
  }
  revisits$pkuid = NA
  revisits$entry_time = NA
  revisits$exit_time = NA
  
  
  
  tempvis = c()
  count = 1
  for(k in 1:nrow(revisits)){
    tempvisit = revisits[k,]

    if(revisits$X2[k]=="Out"){
      tempvis[count] = list(tempvisit)
      count = count+1
      next
    }else{
      Patch = unique(temp$pkuid[seq(from = revisits$index2[k], to = revisits$index3[k], by = 1)])
      for(z in 1:length(Patch)){
        tempvisit$pkuid=Patch[z]
        
        tempvisit$entry_time = as.character(min(temp$timestamp[seq(from = revisits$index2[k], to = revisits$index3[k], by = 1)][which(temp[seq(from = revisits$index2[k], to = revisits$index3[k], by = 1),]$pkuid==Patch[z])], na.rm = TRUE))
        
        tempvisit$exit_time = as.character(max(temp$timestamp[seq(from = revisits$index2[k], to = revisits$index3[k], by = 1)][which(temp[seq(from = revisits$index2[k], to = revisits$index3[k], by = 1),]$pkuid==Patch[z])], na.rm = TRUE))
        
        tempvis[count] = list(tempvisit)
        count = count+1
      }

      

    }
    
  }
  revisits2 = do.call(rbind,tempvis)
  rm(tempvis)
  revisits2$entry_time=as.POSIXct(revisits2$entry_time, format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="America/Panama")
  
  revisits2$exit_time=as.POSIXct(revisits2$exit_time, format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="America/Panama")
  revisits$visit_duration=as.numeric(difftime(revisits$exit_time,revisits$entry_time,units="mins"))
  
  revisits2=revisits2[which(revisits2$X2=="In"),]
  revisits2 = revisits2[
    with(revisits2, order(id, entry_time)),
  ]
  
  revisits2$Between_visit_duration = NA
  revisits2$Visit_no=seq(from = 1, to =nrow(revisits2), by = 1)

  uniquepatches = unique(revisits2$pkuid)
  
  allvisits = c()
  for(l in 1:length(uniquepatches)){
    temppatches = revisits2[which(revisits2$pkuid==uniquepatches[l]),]
    
    for(ll in 2:nrow(temppatches)){
      if(nrow(temppatches)>1){
        temppatches$Between_visit_duration[1]=NA
        temppatches$Between_visit_duration[ll] = as.numeric(difftime(temppatches$entry_time[ll],temppatches$exit_time[ll-1],units="mins"))
      }else{
        temppatches$Between_visit_duration[1]=NA
      }
      
    }
    temppatches$revisit_no=seq(from = 1, to = nrow(temppatches), by = 1)
    allvisits[l] = list(temppatches)
  }
  revisits2 = data.table::rbindlist(allvisits)
  rm(allvisits)
  revisits2 = revisits2[
    with(revisits2, order(id, entry_time)),
  ]
  
  for(m in 2:nrow(revisits2)){
    if( revisits2$Between_visit_duration[m] <12 & revisits2$pkuid[m]==revisits2$pkuid[m-1]){
      revisits2$Visit_no[m] = revisits2$Visit_no[m-1]

    } else{
      revisits2$Visit_no[m] = revisits2$Visit_no[m-1]+1
    }
  }

  
  allvisits = unique(revisits2$Visit_no)
  
  revisits3 = c()
  for(n in 1:length(allvisits)){
    tempvisit = revisits2[which(revisits2$Visit_no==allvisits[n]),]
    tempvisit = tempvisit[
      with(tempvisit, order(id, entry_time)),
    ]
    ID = unique(tempvisit$id)
    Species = unique(tempvisit$species) 
    PKUID = unique(tempvisit$pkuid)
    entry_time = min(tempvisit$entry_time, na.rm = TRUE)
    exit_time = max(tempvisit$exit_time, na.rm = TRUE)
    visit_duration=as.numeric(difftime(exit_time,entry_time,units="mins"))
    entry_time = as.character(entry_time)
    exit_time = as.character(exit_time)
    Between_visit_duration  = tempvisit$Between_visit_duration[1]
    Visit_no = allvisits[n]
    Patch_centroid_x=centroids[which(centroids$pkuid==PKUID),1]
    Patch_centroid_y=centroids[which(centroids$pkuid==PKUID),2]
    
    res = data.frame(cbind(ID,Species,PKUID,Patch_centroid_x,Patch_centroid_y,entry_time,exit_time,visit_duration,Between_visit_duration,Visit_no))
    if(length(res)<10){
      print(n)
      break
    }
    revisits3[n]= list(res)
    
  }
  revisits3 = do.call(rbind,revisits3)
  revisits3$visit_duration = as.numeric(revisits3$visit_duration)
  revisits3$Between_visit_duration = as.numeric(revisits3$Between_visit_duration)
  revisits3$Visit_no = as.numeric(revisits3$Visit_no)
  revisits3$entry_time=as.POSIXct(revisits3$entry_time, format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="America/Panama")
  
  revisits3$exit_time=as.POSIXct(revisits3$exit_time, format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="America/Panama")
  
  #revisits3=revisits3[which(revisits3$visit_duration>=1),]
  revisits3$Visit_no=seq(from = 1, to =nrow(revisits3), by = 1)
  

  
  uniquepatches = unique(revisits3$PKUID)
  
  temppatches=split(revisits3,as.factor(revisits3$PKUID))
  allvisits = c()
  temppatches2=c()
  toremove=c()
  toremovecounter=1
  for(p in 1:length(temppatches)){
    allvisits <-temppatches[[p]]
    allvisits = allvisits[
      with(allvisits, order(ID, entry_time)),
    ]
    #temppatches <- revisits3[which(revisits3$PKUID==uniquepatches[p]),];
    indexrows = nrow(allvisits)
    #print(indexrows);
    for(q in 2:indexrows){
      
      if(indexrows>1){
        
        allvisits$Between_visit_duration[1]=NA;
        allvisits$Between_visit_duration[q] = as.numeric(difftime(allvisits$entry_time[q],allvisits$exit_time[q-1],units="mins"));
        
      }else{
        
        allvisits$Between_visit_duration[1]=NA
        
      }

    }
    
    if(indexrows>1){
      allvisits$revisit_no=seq(from = 1, to = nrow(allvisits), by = 1)
      for(r in 2:nrow(allvisits)){
        if( is.na(allvisits$Between_visit_duration[r])==FALSE & allvisits$Between_visit_duration[r] <12 & allvisits$PKUID[r]==allvisits$PKUID[r-1]){
          allvisits$revisit_no[r] = allvisits$revisit_no[r-1]
          
        } else{
          allvisits$revisit_no[r] = allvisits$revisit_no[r-1]+1
        }
      }
    } else{
      allvisits$revisit_no=NA
    }
    
    allvisits2 = unique(allvisits$revisit_no)
    allrev=c()
 
      #allvisits$revisit_no=seq(from = 1, to = nrow(allvisits), by = 1)
      for(r in 1:length(allvisits2)){
        if(indexrows>1){
        temprevisit = allvisits[which(allvisits$revisit_no==allvisits2[r]),]
        temprevisit = temprevisit[
          with(temprevisit, order(ID, entry_time)),
        ]
        ID = unique(temprevisit$ID)
        Species = unique(temprevisit$Species) 
        PKUID = unique(temprevisit$PKUID)
        entry_time = min(temprevisit$entry_time, na.rm = TRUE)
        exit_time = max(temprevisit$exit_time, na.rm = TRUE)
        visit_duration=as.numeric(difftime(exit_time,entry_time,units="mins"))
        entry_time = as.character(entry_time)
        exit_time = as.character(exit_time)
        Between_visit_duration  = temprevisit$Between_visit_duration[1]
        Visit_no = allvisits2[r]
        Patch_centroid_x=centroids[which(centroids$pkuid==PKUID),1]
        Patch_centroid_y=centroids[which(centroids$pkuid==PKUID),2]
        
        res = data.frame(cbind(ID,Species,PKUID,Patch_centroid_x,Patch_centroid_y,entry_time,exit_time,visit_duration,Between_visit_duration,Visit_no))
        if(length(res)<10){
          print(r)
          break
        }
        allrev[r]= list(res)
        
      }     else{
        res = data.frame(cbind(allvisits$ID,allvisits$Species,allvisits$PKUID,allvisits$Patch_centroid_x,allvisits$Patch_centroid_y,as.character(allvisits$entry_time),as.character(allvisits$exit_time),allvisits$visit_duration,allvisits$Between_visit_duration,allvisits$Visit_no))
        colnames(res)=c("ID","Species","PKUID","Patch_centroid_x","Patch_centroid_y","entry_time","exit_time","visit_duration","Between_visit_duration","Visit_no")
        allrev[r]= list(res)
      }
      }
    
    


    
    allvisits = do.call(rbind,allrev)
    allvisits$visit_duration = as.numeric(allvisits$visit_duration)
    allvisits$Between_visit_duration = as.numeric(allvisits$Between_visit_duration)
    allvisits$Visit_no = as.numeric(allvisits$Visit_no)
    allvisits$entry_time=as.POSIXct(allvisits$entry_time, format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="America/Panama")
    
    allvisits$exit_time=as.POSIXct(allvisits$exit_time, format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="America/Panama")
    
    allvisits = allvisits[
      with(allvisits, order(ID, entry_time)),
    ]
    
    #allvisits=allvisits[which(allvisits$visit_duration>=1),]
    #allvisits$Visit_no=seq(from = 1, to =nrow(allvisits), by = 1)
    
    indexrows = nrow(allvisits)
    #print(indexrows);
    for(q in 2:indexrows){
      
      if(indexrows>1){
        
        allvisits$Between_visit_duration[1]=NA;
        allvisits$Between_visit_duration[q] = as.numeric(difftime(allvisits$entry_time[q],allvisits$exit_time[q-1],units="mins"));
        
      }else{
        
        allvisits$Between_visit_duration[1]=NA
        
      }
      
    }
    
    allvisits$revisit_no=seq(from = 1, to = nrow(allvisits), by = 1)
    
    if(length(which(allvisits$visit_duration>=1))==0){
      toremove[toremovecounter]=p
      toremovecounter=toremovecounter+1
      temppatches2[p] = list(allvisits)
      
    } else{
      allvisits=allvisits[which(allvisits$visit_duration>=1),]
      allvisits$revisit_no=seq(from = 1, to = nrow(allvisits), by = 1)
      temppatches2[p] = list(allvisits)
    }

  }
  if(!is.null(toremove)){
    temppatches2=temppatches2[-toremove]
  }
  revisits3 = data.table::rbindlist(temppatches2)
  rm(toremove)
  rm(toremovecounter)
  rm(temppatches)
  rm(temppatches2)
  rm(allvisits)


  
  revisits3$visit_duration = as.numeric(revisits3$visit_duration)
  revisits3$Between_visit_duration = as.numeric(revisits3$Between_visit_duration)
  revisits3$Visit_no = as.numeric(revisits3$Visit_no)
  revisits3$entry_time=as.POSIXct(revisits3$entry_time, format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="America/Panama")
  
  revisits3$exit_time=as.POSIXct(revisits3$exit_time, format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="America/Panama")
  
  #revisits3=revisits3[which(revisits3$visit_duration>=1),]
  revisits3$Visit_no=seq(from = 1, to =nrow(revisits3), by = 1)
  
  revisits3 = revisits3[
    with(revisits3, order(ID, entry_time)),
  ]
  
  revisits3$Visit_no=seq(from = 1, to =nrow(revisits3), by = 1)
  
  
  revisits3$Time_since_previous_dip_visit = NA
  #revisits3$Between_visit_duration[1] = NA
  for(o in 2:nrow(revisits3)){
    revisits3$Time_since_previous_dip_visit[o] = as.numeric(difftime(revisits3$entry_time[o],revisits3$exit_time[o-1],units="mins"))
    
  }

  
  visitdata[i] = list(revisits3)
}


rm(datasplit2)
visitdata2=data.table::rbindlist(visitdata)
rm(visitdata)
gc()

#visitdata2=do.call(rbind,visitdata)
#visitdata2$visit_duration=visitdata2$visit_duration/60
#visitdata2$Between_visit_duration=visitdata2$Between_visit_duration/60

visitdata2$revisit_ID=paste(visitdata2$PKUID,visitdata2$Visit_no,sep="_")
visitdata2$ID=as.factor(visitdata2$ID)



write.csv(visitdata2, file="revisits_12min_threshold 2023_4_hour.csv")
gc()

vedba=read.csv("FFT veDBA.csv")
#long_visits2=read.csv("revisits_12min_threshold (2).csv")

gc()

vedba$timestamp<-as.POSIXct((vedba$timestamp), format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="UTC")
vedba$localtime<-lubridate::with_tz(vedba$timestamp, tzone = "America/Panama")
vedba$ID<- paste(vedba$individual.local.identifier, vedba$tag.local.identifier, sep = " ")

visitdata2$total_patch_vedba=NA
visitdata2$mean_patch_vedba=NA
visitdata2$vedba_above_35=NA
visitdata2$vedba_below_35=NA
visitdata2$median_patch_vedba=NA

individuals = unique(as.character(visitdata2$ID))


for(i in 1:44){
  print(i)
  indexes = which(as.character(visitdata2$ID)== individuals[i])
  temp = visitdata2[indexes,]
  pb = txtProgressBar(min = 0, max = nrow(temp), style = 3 ) 
  
  for(j in 1:nrow(temp)){

    VEDBAS = vedba[which(vedba$ID==individuals[i]),]
    VEDBAS = VEDBAS$mean_vedba[which(VEDBAS$localtime>=temp$entry_time[j] & VEDBAS$localtime < temp$exit_time[j])]
    
    visitdata2$total_patch_vedba[indexes[j]] = sum(VEDBAS, na.rm = TRUE)
    visitdata2$mean_patch_vedba[indexes[j]] = mean(VEDBAS, na.rm = TRUE)
    visitdata2$median_patch_vedba[indexes[j]] = median(VEDBAS, na.rm = TRUE)
    visitdata2$vedba_above_35[indexes[j]] = length(VEDBAS[which(VEDBAS>=35)])/length(VEDBAS)
    visitdata2$vedba_below_35[indexes[j]] = length(VEDBAS[which(VEDBAS<35)])/length(VEDBAS)
    
    setTxtProgressBar(pb, j)
    
  }
  close(pb) 
}

gc()

#visitdata2$trip_vedba=NA
#visitdata2$trip_vedba_mean=NA



visitdata2$VID=paste(visitdata2$ID,visitdata2$PKUID,sep="_")
visitdata2$VID=as.factor(visitdata2$VID)

visitdata2$date=lubridate::date(visitdata2$entry_time)


visitdata2_split=split(visitdata2,visitdata2$VID)


gc()

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



for(i in 1:length(visitdata2_split)){
  visitdata2_split[[i]]$day=NA
  if(unique(visitdata2_split[[i]]$Species)=="Potos flavus"){
    if(min(unique(lubridate::year(visitdata2_split[[i]]$entry_time)),na.rm=TRUE)<2017){
      Days=nights_2015
      Days2=nights_20152
    }else{
      Days=nights_2017
      Days2=nights_20172
    }
  }else{
    if(min(unique(lubridate::year(visitdata2_split[[i]]$entry_time)),na.rm=TRUE)<2017){
      Days=days_2015
      Days2=days_20152
    }else{
      Days=days_2017
      Days2=days_20172
    }
  }
  for(j in 2:length(Days)){
    visitdata2_split[[i]]$day[which(visitdata2_split[[i]]$entry_time>=Days[j-1] & visitdata2_split[[i]]$entry_time<Days[j])]=Days2[j-1]
  }
  
}

visitdata2=data.table::rbindlist(visitdata2_split)
rm(visitdata2_split)

visitdata2 = visitdata2[
  with(visitdata2, order(ID, entry_time)),
]

visitdata2$year = lubridate::year(visitdata2$entry_time)
visitdata2$study_period=NA
visitdata2$study_period[which(visitdata2$year<2017)]="First"
visitdata2$study_period[which(visitdata2$year>2016)]="Second"


firstday2015=lubridate::date(min(visitdata2$entry_time[which(visitdata2$year==2015)],na.rm=TRUE))-lubridate::days(5)
lastday2015=lubridate::date(max(visitdata2$entry_time[which(visitdata2$year==2016)],na.rm=TRUE))

firstday2017=lubridate::date(min(visitdata2$entry_time[which(visitdata2$year==2017& visitdata2$ID!="Merk 4665")],na.rm=TRUE))-lubridate::days(4)
lastday2017=lubridate::date(max(visitdata2$entry_time[which(visitdata2$year==2018& visitdata2$ID!="Merk 4665")],na.rm=TRUE))

weeks2015=seq(from=firstday2015, to=lastday2015,by="1 day")
weeks2017=seq(from=firstday2017, to=lastday2017,by="1 day")

weeks2015_2=sort(rep(seq(1:20),7))[1:length(weeks2015)]
weeks2017_2=sort(rep(seq(1:28),7))[1:length(weeks2017)]

weeks=cbind(c(as.character(weeks2015),as.character(weeks2017)),c(weeks2015_2,weeks2017_2))
for(i in 1:nrow(weeks)){
  visitdata2$week[which(as.character(lubridate::date(visitdata2$entry_time))==weeks[i,1])]=weeks[i,2]
  
}

write.csv(visitdata2, file="revisits_12min_threshold with vedba2023_4_hour.csv")
gc()

visitdata2$trip_vedba=NA
visitdata2$trip_vedba_mean=NA
visitdata2$trip_vedba_median=NA
visitdata2$trip_distance = NA

load("~/onetofour.RData")
load("~/fivetoeight.RData")
load("~/ninetotwelve.RData")
load("~/thirteentosexteen.RData")
load("~/seventeentoeighteen.RData")
load("~/nonteentotwentytwo.RData")
load("~/twentythree.RData")
load("~/twentyfour.RData")
load("~/twentyfive.RData")
load("~/twentysix.RData")
load("~/twentyseventothirty.RData")
load("~/thirtyonetothirtyfour.RData")
load("~/thirtyninetofourtytwo.RData")
load("~/fourty.RData")
load("~/fourtyone.RData")
load("~/fourtytwo.RData")
load("~/fourtythreetofourtyfour.RData")

thirtyninetofourtytwo=thirtyninetofourtytwo[-which(thirtyninetofourtytwo$ID == "Tony Stark 4659"),]

compileddata = list(onetofour,fivetoeight,ninetotwelve,
                    thirteentosexteen,seventeentoeighteen,
                    nonteentotwentytwo,twentythree, twentyfour, twentyfive, 
                    twentysix,twentyseventothirty,
                    thirtyonetothirtyfour,thirtyninetofourtytwo,fourty, fourtyone,
                    fourtytwo, fourtythreetofourtyfour
)

gc()

rm(list = c("onetofour","fivetoeight","ninetotwelve",
            "thirteentosexteen","seventeentoeighteen",
            "nonteentotwentytwo","twentythree", "twentyfour", "twentyfive", 
            "twentysix","twentyseventothirty",
            "thirtyonetothirtyfour","thirtyninetofourtytwo","fourty", "fourtyone",
            "fourtytwo", "fourtythreetofourtyfour"
))

gc()
compileddata = do.call(rbind, compileddata)

compileddata$loopID = paste(compileddata$ID, compileddata$Study_day, sep = "_")

compileddata$loopID = as.factor(compileddata$loopID)

compileddata$timestamp = lubridate::with_tz(compileddata$timestamp, tzone = "America/Panama")
gc()

datasplit=split(compileddata,compileddata$loopID)
for(i in 1:length(datasplit)){
  mintime = datasplit[[i]]$timestamp[1]+lubridate::hours(4)
  datasplit[[i]]=datasplit[[i]][which(datasplit[[i]]$timestamp<=mintime),]
}
compileddata=data.table::rbindlist(datasplit)
rm(datasplit)
gc()


visitdata2_split=split(visitdata2,visitdata2$ID)
tempvdba=c()
for(i in 1:44){
  print(i)
  temp = visitdata2_split[[i]]
  pb = txtProgressBar(min = 0, max = nrow(temp), style = 3 )
  for(j in 2:nrow(temp)){
    
    VEDBAS = vedba[which(vedba$ID==unique(temp$ID)),]
    VEDBAS = VEDBAS$mean_vedba[which(VEDBAS$localtime>=temp$exit_time[j-1] & VEDBAS$localtime < temp$entry_time[j])]
    
    tempGPS=compileddata[which(compileddata$ID==temp$ID[j] & compileddata$Study_day==temp$day[j]),]
    tempGPS=tempGPS[which(tempGPS$timestamp>temp$exit_time[j-1] & tempGPS$timestamp<temp$entry_time[j]),]

    temp$trip_vedba[j] = sum(VEDBAS, na.rm = TRUE)
    temp$trip_vedba_mean[j] = mean(VEDBAS, na.rm = TRUE)
    temp$trip_vedba_median[j] = median(VEDBAS, na.rm = TRUE)
    
    trip_distance = sum(Mod(diff(tempGPS$V1)), na.rm = TRUE)
    temp$trip_distance[j]=trip_distance
    
    setTxtProgressBar(pb, j)
    
  }
  close(pb)
  tempvdba[i]=list(temp)
}


visitdata2=data.table::rbindlist(tempvdba)
rm(visitdata2_split)
rm(tempvdba)
gc()

visitdata2 = visitdata2[
  with(visitdata2, order(ID, entry_time)),
]

write.csv(visitdata2, file="revisits_12min_threshold with vedba2023_4_hour.csv")
