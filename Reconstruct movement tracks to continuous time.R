library(sp)
library(rgdal)


data=read.csv("FFT_Cleaned_Final.csv")
data$timestamp<-as.POSIXct(as.character(data$timestamp), format="%Y-%m-%d %H:%M:%OS",origin="01-01-1900", tz="UTC")

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

Gridded = data

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

###reconstruct data to continuous time
###Model fitting can be slow, may consider running on cluster 
library(lubridate)
library(ctmm)
library(sp)
tokeep=names(data)
data2=as.telemetry(data, keep = tokeep) 
dt=c(240) #Set sampling rate for variogram
control <- list(method="pNewton",cores=-1)
PROTO <- ctmm(error=TRUE,circle=FALSE)


variograms=lapply(1:45,function(i) variogram(data2[[i]],dt=dt)) #calculate variogram for each animal 
GUESS=lapply(1:45,function(i) ctmm.guess(data2[[i]], CTMM=PROTO,variogram=variograms[[i]],interactive=FALSE)) ##Guess initial parameter values for model fitting
library(crayon)
sims=c()
SIMS=c()
count = 1
for(i in 1:45) {
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
          
          fits=ctmm.fit(data.subset,CTMM=guess,method="pHREML",control=list(method="pNewton",cores=-1))
          
          cat(yellow("Done fitting","\n"))
          
          
          for(k in 1:10){
            print(paste(id,days[j],"sim",k,sep="_"))
            SIM <- simulate(data=data.subset,object=fits, res=1,complete=TRUE,dt=1) #Simulate from movement model at 1 second intervals
            sim=data.frame(SIM)
            index1=which(names(sim)=="longitude")
            index2=which(names(sim)=="latitude")
            sim <- SpatialPointsDataFrame(coords = sim[,c(index1,index2)], data = sim,
                                          proj4string=CRS("+proj=longlat +datum=WGS84"))
            sim <- spTransform(sim, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
            sim=data.frame(sim)
            sims[k]=list(sim)
          }
          
          sims10=do.call(rbind,sims[c(1:10)]) ##Pull all the simulations for an individual together
          sims10=as.data.frame(sims10)
          sims10=sims10[order(sims10$timestamp),] ##sorts by timestamp
          DT <- data.table::data.table(sims10)
          DT$z=DT$longitude.1+1i*DT$latitude.1 ## make coordinates complex numbers for convenience 
          averagesim10=DT[, mean(z,na.rm=T), by = timestamp] #Average all simulations for an individual, summarizes data with mean of each timestamp
          averagesim10=as.data.frame(averagesim10)
          averagesim10$ID=data2[[i]]@info$identity ##annotates all locations with the indiviual being analyzed, for later refenerence when individual are combined
          averagesim10$X=Re(averagesim10$V1) #gets x (real) component of complex mean of locations
          averagesim10$Y=Im(averagesim10$V1) #gets y (imaginary) component of complex mean of locations
          averagesim10$Species = Species
          averagesim10$Study_day = days[j]
          SIMS[count]=list(averagesim10)
          
          
          count = count + 1
        }

    
  }
  #res=as.data.frame(do.call(rbind,res))
  #   speeds[i]=list(res)
}

data=do.call(rbind,SIMS)
