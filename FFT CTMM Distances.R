library(sp)
library(rgdal)


data=read.csv("FFT_Cleaned_Final.csv")
data$timestamp<-as.POSIXct(as.character(data$timestamp), format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="UTC")
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

count=1
fixes=c()
fixesperday=c()
for(i in 1:length(Gridded)){
  temp = Gridded[[i]]
  temp = split(temp, as.factor(as.character(temp$day)))
  for (j in 1:length(temp)){
    fixes = nrow(temp[[j]])
    studyday = unique(temp[[j]]$day)
    ID = unique(temp[[j]]$individual.local.identifier)
    Species = unique(temp[[j]]$individual.taxon.canonical.name)
    tobind = data.frame(cbind(fixes, studyday, ID, Species))
    fixesperday[count] = list(tobind)
    count = count+1
  }
  
}
fixesperday = do.call(rbind, fixesperday)
fixesperday$fixes = as.numeric(fixesperday$fixes)
fixesperday$studyday = as.numeric(fixesperday$studyday)
fixesperday=fixesperday[-which(fixesperday$ID=="Merk 4665"),]


data=do.call(rbind,Gridded)
data=data[-which(data$individual.local.identifier=="Atlas 4673"),]
data=data[-which(data$individual.local.identifier=="Merk 4665"),]
data=data[-which(data$individual.local.identifier=="Serge 4670"),]
data=data[-which(data$individual.local.identifier=="Sahti 4693"),]
data=data[-which(data$individual.local.identifier=="Chloe 4052"),]

###reconstruct data to continuous time
###Model fitting can be slow, may consider running on cluster in future
library(lubridate)
library(ctmm)
library(sp)
tokeep=names(data)
data2=as.telemetry(data, keep = tokeep) 
dt=c(240) #Set sampling rate for variogram
control <- list(method="pNewton",cores=-1)
PROTO <- ctmm(error=TRUE,circle=FALSE)


variograms=lapply(1:44,function(i) variogram(data2[[i]],dt=dt)) #calculate variogram for each animal 
#GUESS=lapply(1:45,function(i) ctmm.guess(data2[[i]], CTMM=PROTO,variogram=variograms[[i]],interactive=FALSE)) ##Guess initial parameter values for model fitting
library(crayon)
res=c()  ##I create this vector that I will use in my loop
speeds=c() ##I will store things here later
count=1 
for(i in 1:44) {
  days=unique(data2[[i]]$day) 
  Species=unique(data2[[i]]$individual.taxon.canonical.name) 
  for(j in 1:length(days)){
    #select data for day#
    data.subset= data2[[i]][which(data2[[i]]$day==days[j]),]
    if(nrow(data.subset)<135){
      message("Not enough data for day ",days[j]," : ",date, " for ", Species," ", data2[[i]]@info$identity)
      message("Only ",nrow(data.subset)," fixes")
      next
    } else {
      date=as.character(as.Date(data.subset$timestamp[1]))
      id=data2[[i]]@info$identity
      #calculating duration of sampling period#
      samp.time=diff(c(data.subset$t[1],data.subset$t[nrow(data.subset)]))
      guess=ctmm.guess(data2[[i]],CTMM=PROTO,variogram=variograms[[i]],interactive=FALSE)
      cat(yellow("Fitting model","\n"))
      
      #op <- options(warn=2)
      
      fits=try(ctmm.select(data.subset,CTMM=guess,method="pHREML",trace = TRUE, control=list(method="pNewton",cores=-1)))
      if(class(fits)=="try-error"){
        message("Failed to fit model")
        x=data.frame(cbind(days[j],NA,NA,NA,NA,date,id,as.character(Species)))
        names(x)=c("day","ctmm.estimate","ctmm.min","ctmm.max","tau","date","id","species")
        res[count]=list(x)
        count = count + 1
        next

      }
      cat(yellow("Done fitting","\n"))
      #op <- options(warn=0)
      message("estimating distance traveled on day ",j," :",date, " for ", Species," ", data2[[i]]@info$identity)
      ctmm_speed=try(speed(object=data.subset,CTMM=fits,units=FALSE,robust=TRUE))
      if(class(ctmm_speed)=="try-error"){
        ctmm_speed=try(speed(object=data.subset,CTMM=fits,units=FALSE,robust=TRUE,fast=FALSE))
        if(class(ctmm_speed)=="try-error"){
          message("Failed to estimate speed")
          x=data.frame(cbind(days[j],NA,NA,NA,length(fits$tau),date,id,as.character(Species)))
          names(x)=c("day","ctmm.estimate","ctmm.min","ctmm.max","tau","date","id","species")
          res[count]=list(x)
          count = count + 1
          next
        }
      } else{
        ctmm_dist=ctmm_speed$CI*samp.time
        print(ctmm_dist[2])
        rownames(ctmm_dist)="distance(meters)"
        x=data.frame(cbind(days[j],ctmm_dist[2],ctmm_dist[1],ctmm_dist[3],length(fits$tau),date,id,as.character(Species)))
        names(x)=c("day","ctmm.estimate","ctmm.min","ctmm.max","tau","date","id","species")
        res[count]=list(x)
        count = count + 1
      }


    }
    
    
  }
  #res=as.data.frame(do.call(rbind,res))
  #   speeds[i]=list(res)
}

speeds=data.table::rbindlist(res)
speeds=as.data.frame(speeds)
write.csv(speeds,file="DailyPathLengths CTMM FFT 2023.csv")
