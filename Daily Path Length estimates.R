#Load interpolated data (data run in parallel in chunks)
load("~/workspaces/onetofour.RData")
load("~/workspaces/fivetoeight.RData")
load("~/workspaces/ninetotwelve.RData")
load("~/workspaces/thirteentosexteen.RData")
load("~/workspaces/seventeentoeighteen.RData")
load("~/workspaces/nonteentotwentytwo.RData")
load("~/workspaces/twentythree.RData")
load("~/workspaces/twentyfour.RData")
load("~/workspaces/twentyfive.RData")
load("~/workspaces/twentysix.RData")
load("~/workspaces/twentyseventothirty.RData")
load("~/workspaces/thirtyonetothirtyfour.RData")
load("~/workspaces/thirtyninetofourtytwo.RData")
load("~/workspaces/fourty.RData")
load("~/workspaces/fourtyone.RData")
load("~/workspaces/fourtytwo.RData")
load("~/workspaces/fourtythreetofourtyfour.RData")

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

gc()

#Add study day labels 
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



splitdata = split(compileddata, compileddata$ID)

for(i in 1:length(splitdata)){
  splitdata[[i]]$Study_day=NA
  if(unique(splitdata[[i]]$Species)=="Potos flavus"){
    if(min(unique(lubridate::year(splitdata[[i]]$timestamp)),na.rm=TRUE)<2017){
      Days=nights_2015
      Days2=nights_20152
    }else{
      Days=nights_2017
      Days2=nights_20172
    }
  }else{
    if(min(unique(lubridate::year(splitdata[[i]]$timestamp)),na.rm=TRUE)<2017){
      Days=days_2015
      Days2=days_20152
    }else{
      Days=days_2017
      Days2=days_20172
    }
  }
  for(j in 2:length(Days)){
    splitdata[[i]]$Study_day[which(splitdata[[i]]$timestamp>=Days[j-1] & splitdata[[i]]$timestamp<Days[j])]=Days2[j-1]
  }
  
}
gc()
compileddata = data.table::rbindlist(splitdata)
compileddata$loopID = paste(compileddata$ID, compileddata$Study_day, sep = "_")

splitdata = split(compileddata, compileddata$loopID)

gc()
##Calculate path length if each day for interpolated data
interpolatedpathlength = c()
for(i in 1:length(splitdata)){
  temp = splitdata[[i]]
  ID = unique(temp$ID)
  Study_day = unique(temp$Study_day)
  Species = unique(temp$Species)
  pathlength= sum(Mod(diff(temp$V1)), na.rm = TRUE)
  interpolatedpathlength[i] = list(data.frame(cbind(ID, Species, Study_day, pathlength)))
  
}

interpolatedpathlength = data.table::rbindlist(interpolatedpathlength)


#Calculate path length of observed data without interpolations
data=read.csv("FFT_Cleaned_Final.csv")
data$timestamp<-as.POSIXct(as.character(data$timestamp), format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="UTC")

Gridded = split(data,as.factor(data$individual.local.identifier))

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



data=data.table::rbindlist(Gridded)
data=data[-which(data$individual.local.identifier=="Atlas 4673"),]
data=data[-which(data$individual.local.identifier=="Merk 4665"),]
data=data[-which(data$individual.local.identifier=="Serge 4670"),]
data=data[-which(data$individual.local.identifier=="Sahti 4693"),]
data=data[-which(data$individual.local.identifier=="Chloe 4052"),]

data$loopID = paste(trimws(as.character(data$individual.local.identifier)), trimws(as.character(data$day)), sep = "_")
data$loopID = as.factor(data$loopID)   
library(sp)
data2 <- SpatialPointsDataFrame(coords = data[,c(3,4)], data = data,
                                proj4string=CRS("+proj=longlat +datum=WGS84"))
data2 <- spTransform(data2, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
data2=data.frame(data2)
data2$Z = data2$location.long.2 +1i*data2$location.lat.2

splitdata2 = split(data2, data2$loopID)

originalpathlength = c()
for(i in 1:length(splitdata2)){
  #temp = splitdata2[[i]]
  ID = unique(splitdata2[[i]]$individual.local.identifier)
  Study_day = unique(splitdata2[[i]]$day)
  date = as.character(min(unique(lubridate::date(splitdata2[[i]]$timestamp))))
  Species = unique(splitdata2[[i]]$individual.taxon.canonical.name)
  pathlength= sum(Mod(diff(splitdata2[[i]]$Z)), na.rm = TRUE)
  originalpathlength[i] = list(data.frame(cbind(ID, Species, Study_day,date, pathlength)))
  
}

originalpathlength = do.call(rbind,originalpathlength)

pathlengthinfo = merge(originalpathlength, interpolatedpathlength, by.x = c("ID", "Species","Study_day"), by.y = c("ID", "Species","Study_day"), all.x = TRUE)

pathlengthinfo$pathlength.x = as.numeric(pathlengthinfo$pathlength.x)
pathlengthinfo$pathlength.y = as.numeric(pathlengthinfo$pathlength.y)

pathlengthinfo = pathlengthinfo[-which(pathlengthinfo$pathlength.y>10000),]

library(ggplot2)
plot1 = ggplot(pathlengthinfo, aes(pathlength.x, pathlength.y, color = Species)) + geom_point() + coord_equal()
plot1+facet_wrap(~factor(ID))

#Load and add precomputed daily path length using the CTMM method
ctmmspeed = read.csv("C:/Users/salavi/Documents/DailyPathLengths CTMM FFT 2023.csv")
ctmmspeed$date = lubridate::as_datetime(ctmmspeed$date)
ctmmspeed=ctmmspeed[,-1]

pathlengthinfo2 = merge(pathlengthinfo, ctmmspeed, by.x = c("ID", "Species","Study_day"), by.y = c("id", "species","day"), all.x = TRUE)


pathlengthinfo2=pathlengthinfo2[,c(-11)]
colnames(pathlengthinfo2)=c("ID","Species","Study_day","date","Straight_line_DJL","Interpolated_DJL","ctmm.estimate","ctmm.min","ctmm.max","tau")

plot2 = ggplot(pathlengthinfo2, aes(Straight_line_DJL, ctmm.estimate, color = Species)) + geom_point() + coord_equal()
plot2+facet_wrap(~factor(ID))

plot3 = ggplot(pathlengthinfo2, aes(Interpolated_DJL, ctmm.estimate, color = Species)) + geom_point() + coord_equal()
plot3+facet_wrap(~factor(ID))


ggplot(pathlengthinfo2, aes(Straight_line_DJL, ctmm.estimate, color = Species)) + geom_point() +
  geom_errorbar(aes(ymin=ctmm.min, ymax=ctmm.max), width=.2,
              position=position_dodge(.9))+ coord_equal()
library(brms)
library(cmdstanr)

write.csv(pathlengthinfo3, file = "FFT Daily Path Length2.csv")

#Look at crude uncertainties of path lengths
pathlengthinfo2$uncertainty_ctmm = pathlengthinfo2$ctmm.max-pathlengthinfo2$ctmm.min
pathlengthinfo2$uncertainty_SLD=NA
pathlengthinfo2$uncertainty_interp=NA

#Check if ctmm error bars include the estimates from the other two approaches
pathlengthinfo2[which(pathlengthinfo2$Straight_line_DJL>pathlengthinfo2$ctmm.min & pathlengthinfo2$Straight_line_DJL< pathlengthinfo2$ctmm.max & !is.infinite(pathlengthinfo2$ctmm.min) & !is.infinite(pathlengthinfo2$ctmm.max)),]$uncertainty_SLD = "In"
pathlengthinfo2[-which(pathlengthinfo2$Straight_line_DJL>pathlengthinfo2$ctmm.min & pathlengthinfo2$Straight_line_DJL< pathlengthinfo2$ctmm.max & !is.infinite(pathlengthinfo2$ctmm.min) & !is.infinite(pathlengthinfo2$ctmm.max)),]$uncertainty_SLD = "Out"

pathlengthinfo2[which(pathlengthinfo2$Interpolated_DJL>pathlengthinfo2$ctmm.min & pathlengthinfo2$Interpolated_DJL< pathlengthinfo2$ctmm.max & !is.infinite(pathlengthinfo2$ctmm.min) & !is.infinite(pathlengthinfo2$ctmm.max)) ,]$uncertainty_interp = "In"

pathlengthinfo2[-which(pathlengthinfo2$Interpolated_DJL>pathlengthinfo2$ctmm.min & pathlengthinfo2$Interpolated_DJL< pathlengthinfo2$ctmm.max & !is.infinite(pathlengthinfo2$ctmm.min) & !is.infinite(pathlengthinfo2$ctmm.max)) ,]$uncertainty_interp = "Out"

pathlengthinfo2[which(is.infinite(pathlengthinfo2$ctmm.min) | is.infinite(pathlengthinfo2$ctmm.max)),]$uncertainty_SLD = NA
pathlengthinfo2[which(is.infinite(pathlengthinfo2$ctmm.min) | is.infinite(pathlengthinfo2$ctmm.max)),]$uncertainty_interp = NA

pathlengthinfo3=pathlengthinfo2[complete.cases(pathlengthinfo2$Interpolated_DJL),]
pathlengthinfo3=pathlengthinfo3[complete.cases(pathlengthinfo3$ID),]

pathlengthinfo3$uncertainty_ctmm2 = pathlengthinfo3$uncertainty_ctmm/2


ggplot(pathlengthinfo3, aes(Straight_line_DJL, ctmm.estimate, color = uncertainty_SLD)) + geom_point()+ facet_wrap(~ID)

ggplot(pathlengthinfo3, aes(Interpolated_DJL, ctmm.estimate, color = uncertainty_interp)) + geom_point()+ facet_wrap(~ID)

length(which(pathlengthinfo3$uncertainty_ctmm2>=500))
length(which(pathlengthinfo3$uncertainty_ctmm2<500))/nrow(pathlengthinfo3)
length(which(is.infinite(pathlengthinfo3$ctmm.min) | is.infinite(pathlengthinfo3$ctmm.max)))
length(which(pathlengthinfo3$uncertainty_SLD=="In"))/nrow(pathlengthinfo3)
length(which(pathlengthinfo3$uncertainty_SLD=="Out"))
length(which(pathlengthinfo3$uncertainty_interp=="In"))/nrow(pathlengthinfo3)
length(which(pathlengthinfo3$uncertainty_interp=="Out"))

pathlengthinfo3$Straight_line_DJL[1]-pathlengthinfo3$Interpolated_DJL[1]
pathlengthinfo3$Straight_line_DJL[2]-pathlengthinfo3$Interpolated_DJL[2]

#Check actual relationship between SLD and interpolated estimates
options(mc.cores = parallel::detectCores()) 

distancemodel=brm(pathlength.y ~ pathlength.x, 
                    data= data3,
                    family = categorical(link = "logit"),
                    iter = 4000,
                    prior = c(
                      prior(normal(0, 1.5), class = Intercept, dpar = muM5),
                      prior(normal(0, 1.5), class = Intercept, dpar = muM10),
                      prior(normal(0, 1.5), class = Intercept, dpar = muM15),
                      prior(normal(0, 1.5), class = Intercept, dpar = muM20),
                      prior(normal(0, 1.5), class = Intercept, dpar = muM30)
                    ),
                    save_pars = save_pars(all = TRUE),
                    backend = "cmdstanr",
                    threads = threading(8),
                    control = list(max_treedepth = 10, adapt_delta = .999))
