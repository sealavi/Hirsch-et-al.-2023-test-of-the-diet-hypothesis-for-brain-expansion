library(sp)
library(rgdal)
library(lubridate)
library(sp)
library(ctmm)
library(sf)
library(rgeos)

load("~/RanSimFFT workspace.RData") #load precomputed OUF models
data=read.csv("C:/Users/salavi/Documents/FFT_Cleaned_Final.csv")
data$timestamp=as.POSIXct(as.character(data$timestamp),format="%Y-%m-%d %H:%M:%S",origin="01-01-1900",tz="UTC")

###creat study day labels
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

##Potos has to be treated different because they are nocturnal
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
rm(Gridded)

##Fit continuous time movement models
data2=as.telemetry(data)  
dt=c(240) #Set sampling rate for variogram
control <- list(method="pNewton",cores=-1)
PROTO <- ctmm(error=TRUE,circle=FALSE)
PROTO2 <- ctmm(error=TRUE,mean="periodic",period=c(24 %#% "hours",1 %#% "month"),circle=FALSE)



variograms=lapply(1:49,function(i) variogram(data2[[i]],dt=dt)) #Estimate variograms for each animal



BM=ctmm(tau = Inf) #Parameters for brownian movement
IID=ctmm(tau = NULL) #Parameters for IID (totally unrealistic probably wont use)
#OUF models computed previously for analysis in previous paper. Loaded earlier
GUESS=lapply(1:49,function(i) ctmm.guess(data2[[i]], CTMM=PROTO,variogram=variograms[[i]],interactive=FALSE)) ##Guess starting values for parameters

GUESS2=lapply(1:49,function(i) ctmm.guess(data2[[i]], CTMM=PROTO2,variogram=variograms[[i]],interactive=FALSE))

#Fit IID models
FITS_IID=lapply(1:49,function(i) {
  print(i)
  ctmm.fit(data2[[i]],IID,control=control)})

#Fit Brownian motion models
FITS_BM=lapply(1:49,function(i) {
  print(i)
  ctmm.fit(data2[[i]],BM,control=control)})

gc()

#Create unique identifier for each animal on each day
data$loopID = paste(data$individual.local.identifier, data$day, sep = "_")

datasplit = split(data, as.factor(data$individual.local.identifier))

#create vector of times to use in simulations so they match the day lengths of observed data
times = c()
times2 = c()
for(i in 1:length(datasplit)){
  times3=c()
  temptimes = c()
  tempdata = split(datasplit[[i]], as.factor(datasplit[[i]]$loopID))
  
  for(j in 1:length(tempdata)){
    timestamps = seq(from = lubridate::floor_date(min(tempdata[[j]]$timestamp, na.rm = TRUE), unit = "mins"), to = lubridate::floor_date(max(tempdata[[j]]$timestamp, na.rm = TRUE ), unit = "mins"), by = "240 sec")
    t = seq(1:length(timestamps))
    
    temptimes[j] = list(t)
    times3[j] = list(as.character(timestamps))
  }
  temptimes=unlist(temptimes)
  times[i] = list(unlist(times3))
  times2[i] = length(temptimes)
}
  
#Load island and tree crown shapefiles  
BCI = readOGR("C:/Users/salavi/Documents/BCI_outline.shp")
BCI = sp::spTransform(BCI, sp::CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))

crowns= readOGR("C:/Users/salavi/Downloads/documents-export-2022-03-23/BCI_Dipteryx_Master_final.shp")
crowns = spTransform(crowns, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
centroids=geosphere::centroid(crowns) #calculate centroid fo crowns 
centroids=data.frame(centroids)
centroids$pkuid=crowns$pkuid

##Create buffers around crowns to vary distance at which a tree can be encountered
crowns10=gBuffer(crowns, width = 10, byid = TRUE)
crowns50=gBuffer(crowns, width = 50, byid = TRUE)
crowns100=gBuffer(crowns, width = 100, byid = TRUE)
crowns150=gBuffer(crowns, width = 150, byid = TRUE)
crowns200=gBuffer(crowns, width = 200, byid = TRUE)


## Simulate random movement using models 
simulations_IID = c()
simulations_BM = c()
simulations_OUF = c()

count = 1

for(i in 1:49){
  print(i)
  temp = datasplit[[i]]
  
  #calculate midpoint of homerange 
  midpoint=temp
  midpoint <- sp::SpatialPointsDataFrame(coords = midpoint[,c(3,4)], data = midpoint,
                                        proj4string=sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  midpoint <- spTransform(midpoint, sp::CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
  midpoint=data.frame(midpoint)
  midpoint$Z=midpoint$location.long.2+1i*midpoint$location.lat.2
  midpoint=mean(midpoint$Z, na.rm = TRUE)
  days = unique(na.omit(temp$day))
 # pb = txtProgressBar(min = 0, max = length(days), style = 3 )
  # for(j in 1:length(days)){
    # temp2 = temp[which(temp$day==days[j]),]
    timestamps = as.POSIXct(times[[i]],format="%Y-%m-%d %H:%M:%S",origin="01-01-1900",tz="UTC") #Appropriate vector of times to simulate along

    #Simulate IID model
    t = seq(1:times2[i])
    SIM_IID = simulate(FITS_IID[[i]],t=t, dt = 240)
    SIM_IID = data.frame(SIM_IID)
    SIM_IID$timestamp=timestamps
    #SIM_IID$day=temp$day
    SIM_IID$ID=unique(temp$individual.local.identifier)
    SIM_IID$Species=unique(temp$individual.taxon.canonical.name)
    SIM_IID$z=SIM_IID$x+1i*SIM_IID$y
    SIM_IID$z_utm=SIM_IID$z+midpoint #Ensure data are within actual homerange
    SIM_IID$utm_E=Re(SIM_IID$z_utm)
    SIM_IID$utm_N=Im(SIM_IID$z_utm)
    
    #Simulate Brownian movement model
    SIM_BM = simulate(FITS_BM[[i]],t=t, dt = 240)
    SIM_BM = data.frame(SIM_BM)
    SIM_BM$timestamp=timestamps
    #SIM_BM$day=temp$day
    SIM_BM$ID=unique(temp$individual.local.identifier)
    SIM_BM$Species=unique(temp$individual.taxon.canonical.name)
    SIM_BM$z=SIM_BM$x+1i*SIM_BM$y
    SIM_BM$z_utm=SIM_BM$z+midpoint #Ensure data are within actual homerange
    SIM_BM$utm_E=Re(SIM_BM$z_utm)
    SIM_BM$utm_N=Im(SIM_BM$z_utm)
    
    timestamps2=seq(from = min(timestamps), to = max(timestamps), by = "10 sec")
    t2 = seq(1:length(timestamps2))
    
    #Simulate Ornstein Eulenbeck F movement
    SIM_OUF = simulate(OUF_FITS[[i]],t=t2, dt = 10)
    SIM_OUF = data.frame(SIM_OUF)
    SIM_OUF$timestamp=timestamps2
    SIM_OUF = SIM_OUF[which(SIM_OUF$timestamp%in%timestamps==TRUE),]
    #SIM_OUF$day=temp$day
    SIM_OUF$ID=unique(temp$individual.local.identifier)
    SIM_OUF$Species=unique(temp$individual.taxon.canonical.name)
    SIM_OUF$z=SIM_OUF$x+1i*SIM_OUF$y
    SIM_OUF$z_utm=SIM_OUF$z+midpoint #Ensure data are within actual homerange
    SIM_OUF$utm_E=Re(SIM_OUF$z_utm)
    SIM_OUF$utm_N=Im(SIM_OUF$z_utm)


    simulations_IID[i]= list(SIM_IID)
    simulations_BM[i] = list(SIM_BM)
    simulations_OUF[i] = list(SIM_OUF)
    
    #count = count+1
  #   setTxtProgressBar(pb, j)
  # }
  # close(pb) 
  
  }
gc()

##Create dataframe of tree encounters (overlap with crown or not) for each buffer level
IID=c()
BM=c()
OUF=c()
count=1
for(i in 1:49){
  temp = simulations_IID[[i]]
  temp <- sp::SpatialPointsDataFrame(coords = temp[,c(9,10)], data = temp,
                                                     proj4string=sp::CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
  SIM_IIDpoly=over(temp,crowns)
  SIM_IIDpoly10=over(temp,crowns10)
  SIM_IIDpoly50=over(temp,crowns50)
  SIM_IIDpoly100=over(temp,crowns100)
  SIM_IIDpoly150=over(temp,crowns150)
  SIM_IIDpoly200=over(temp,crowns200)
  
  

  temp = data.frame(temp)
  temp10=cbind(temp,SIM_IIDpoly10)
  temp50=cbind(temp,SIM_IIDpoly50)
  temp100=cbind(temp,SIM_IIDpoly100)
  temp150=cbind(temp,SIM_IIDpoly150)
  temp200=cbind(temp,SIM_IIDpoly200)
  temp=cbind(temp,SIM_IIDpoly)
  temp$buffer = 0
  temp10$buffer = 10
  temp50$buffer = 50
  temp100$buffer = 100
  temp150$buffer = 150
  temp200$buffer = 200
  
  tempBM = simulations_BM[[i]]
  tempBM <- sp::SpatialPointsDataFrame(coords = tempBM[,c(9,10)], data = tempBM,
                                     proj4string=sp::CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
  SIM_BMpoly=over(tempBM,crowns)
  SIM_BMpoly10=over(tempBM,crowns10)
  SIM_BMpoly50=over(tempBM,crowns50)
  SIM_BMpoly100=over(tempBM,crowns100)
  SIM_BMpoly150=over(tempBM,crowns150)
  SIM_BMpoly200=over(tempBM,crowns200)
  
  
  
  tempBM = data.frame(tempBM)
  tempBM10=cbind(tempBM,SIM_BMpoly10)
  tempBM50=cbind(tempBM,SIM_BMpoly50)
  tempBM100=cbind(tempBM,SIM_BMpoly100)
  tempBM150=cbind(tempBM,SIM_BMpoly150)
  tempBM200=cbind(tempBM,SIM_BMpoly200)
  tempBM=cbind(tempBM,SIM_BMpoly)
  tempBM$buffer = 0
  tempBM10$buffer = 10
  tempBM50$buffer = 50
  tempBM100$buffer = 100
  tempBM150$buffer = 150
  tempBM200$buffer = 200
  
  
  tempOUF = simulations_OUF[[i]]
  tempOUF <- sp::SpatialPointsDataFrame(coords = tempOUF[,c(11,12)], data = tempOUF,
                                       proj4string=sp::CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
  SIM_OUFpoly=over(tempOUF,crowns)
  SIM_OUFpoly10=over(tempOUF,crowns10)
  SIM_OUFpoly50=over(tempOUF,crowns50)
  SIM_OUFpoly100=over(tempOUF,crowns100)
  SIM_OUFpoly150=over(tempOUF,crowns150)
  SIM_OUFpoly200=over(tempOUF,crowns200)
  
  
  
  tempOUF = data.frame(tempOUF)
  tempOUF10=cbind(tempOUF,SIM_OUFpoly10)
  tempOUF50=cbind(tempOUF,SIM_OUFpoly50)
  tempOUF100=cbind(tempOUF,SIM_OUFpoly100)
  tempOUF150=cbind(tempOUF,SIM_OUFpoly150)
  tempOUF200=cbind(tempOUF,SIM_OUFpoly200)
  tempOUF=cbind(tempOUF,SIM_OUFpoly)
  tempOUF$buffer = 0
  tempOUF10$buffer = 10
  tempOUF50$buffer = 50
  tempOUF100$buffer = 100
  tempOUF150$buffer = 150
  tempOUF200$buffer = 200
  
  
  IID[count]=list(temp)
  BM[count]=list(tempBM)
  OUF[count]=list(tempOUF)
  
  count=count+1
  IID[count]=list(temp10)
  BM[count]=list(tempBM10)
  OUF[count]=list(tempOUF10)
  
  count=count+1
  IID[count]=list(temp50)
  BM[count]=list(tempBM50)
  OUF[count]=list(tempOUF50)
  
  count=count+1
  IID[count]=list(temp100)
  BM[count]=list(tempBM100)
  OUF[count]=list(tempOUF100)
  
  count=count+1
  IID[count]=list(temp150)
  BM[count]=list(tempBM150)
  OUF[count]=list(tempOUF150)
  
  count=count+1
  IID[count]=list(temp200)
  BM[count]=list(tempBM200)
  OUF[count]=list(tempOUF200)
  
  count=count+1
}


gc()


#Add study day labels to simulated data

for(i in 1:length(IID)){
  IID[[i]]$day=NA
  if(unique(IID[[i]]$ID)=="Potos flavus"){
    if(min(unique(lubridate::year(IID[[i]]$timestamp)),na.rm=TRUE)<2017){
      Days=nights_2015
      Days2=nights_20152
    }else{
      Days=nights_2017
      Days2=nights_20172
    }
  }else{
    if(min(unique(lubridate::year(IID[[i]]$timestamp)),na.rm=TRUE)<2017){
      Days=days_2015
      Days2=days_20152
    }else{
      Days=days_2017
      Days2=days_20172
    }
  }
  for(j in 2:length(Days)){
    IID[[i]]$day[which(IID[[i]]$timestamp>=Days[j-1] & IID[[i]]$timestamp<Days[j])]=Days2[j-1]
  }
  
}

for(i in 1:length(BM)){
  BM[[i]]$day=NA
  if(unique(BM[[i]]$ID)=="Potos flavus"){
    if(min(unique(lubridate::year(BM[[i]]$timestamp)),na.rm=TRUE)<2017){
      Days=nights_2015
      Days2=nights_20152
    }else{
      Days=nights_2017
      Days2=nights_20172
    }
  }else{
    if(min(unique(lubridate::year(BM[[i]]$timestamp)),na.rm=TRUE)<2017){
      Days=days_2015
      Days2=days_20152
    }else{
      Days=days_2017
      Days2=days_20172
    }
  }
  for(j in 2:length(Days)){
    BM[[i]]$day[which(BM[[i]]$timestamp>=Days[j-1] & BM[[i]]$timestamp<Days[j])]=Days2[j-1]
  }
  
}

for(i in 1:length(OUF)){
  OUF[[i]]$day=NA
  if(unique(OUF[[i]]$ID)=="Potos flavus"){
    if(min(unique(lubridate::year(OUF[[i]]$timestamp)),na.rm=TRUE)<2017){
      Days=nights_2015
      Days2=nights_20152
    }else{
      Days=nights_2017
      Days2=nights_20172
    }
  }else{
    if(min(unique(lubridate::year(OUF[[i]]$timestamp)),na.rm=TRUE)<2017){
      Days=days_2015
      Days2=days_20152
    }else{
      Days=days_2017
      Days2=days_20172
    }
  }
  for(j in 2:length(Days)){
    OUF[[i]]$day[which(OUF[[i]]$timestamp>=Days[j-1] & OUF[[i]]$timestamp<Days[j])]=Days2[j-1]
  }
  
}

IID=data.table::rbindlist(IID)
BM=data.table::rbindlist(BM)
OUF=data.table::rbindlist(OUF)

#data=data.table::rbindlist(Gridded)

#save(simulations_IID, file = "//10.126.19.90/salavi/simulations_IID.RData")
# rm(simulations_IID)
# gc()

# 
# simulations_IID=data.table::rbindlist(simulations_IID)
# simulations_BM=data.table::rbindlist(simulations_BM)
# simulations_OUF=data.table::rbindlist(simulations_OUF)

IID$timestamp = lubridate::with_tz(IID$timestamp, tzone = "America/Panama")
BM$timestamp = lubridate::with_tz(BM$timestamp, tzone = "America/Panama")
OUF$timestamp = lubridate::with_tz(OUF$timestamp, tzone = "America/Panama")

gc()

IID$visit=NA
BM$visit=NA
OUF$visit=NA

gc()

IID$LoopID = paste(IID$ID, IID$day, IID$buffer, sep = "_")
BM$LoopID = paste(BM$ID, BM$day, BM$buffer, sep = "_")
OUF$LoopID = paste(OUF$ID, OUF$day, OUF$buffer, sep = "_")

#simulations_IID$visit=NA

datasplit2=split(IID,as.factor(IID$LoopID))
gc()

##Calculate diptryx encounter rates per day for IID data
count=1
IID_visitdata=c()
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
  
  if(length(temp$visit[which(temp$visit=="In")])==0){
    ID = unique(temp$ID)
    Species = unique(temp$Species) 
    Day = unique(temp$day)
    Model = "IID"
    buffer = unique(temp$buffer)
    Visits = 0
    
    res = data.frame(cbind(ID,Species,Day,Model,buffer,Visits))
    IID_visitdata[i]= list(res)
    next
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
  if(nrow(revisits)>1){
    for(j in 2:nrow(revisits)){
      revisits$index2[j]= revisits$index3[j-1]+1
      revisits$index3[j]= revisits$index2[j]+(revisits$X1[j]-1)
    } 
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
  revisits2$Between_visit_duration = NA
  revisits2$Visit_no=seq(from = 1, to =nrow(revisits2), by = 1)
  
  ID = unique(tempvisit$id)
  Species = unique(tempvisit$species) 
  Day = unique(temp$day)
  Model = "IID"
  buffer = unique(temp$buffer)
  Visits = nrow(revisits2)

  res = data.frame(cbind(ID,Species,Day,Model,buffer,Visits))
  IID_visitdata[i]= list(res)
  }

IID_visitdata = data.table::rbindlist(IID_visitdata)




datasplit2=split(BM,as.factor(BM$LoopID))
gc()

##Calculate diptryx encounter rates per day for Brownian motion data
#Could have made this a function but just modifying the same loop as before

count=1
BM_visitdata=c()
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
  
  if(length(temp$visit[which(temp$visit=="In")])==0){
    ID = unique(temp$ID)
    Species = unique(temp$Species) 
    Day = unique(temp$day)
    Model = "BM"
    buffer = unique(temp$buffer)
    Visits = 0
    
    res = data.frame(cbind(ID,Species,Day,Model,buffer,Visits))
    BM_visitdata[i]= list(res)
    next
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
  if(nrow(revisits)>1){
    for(j in 2:nrow(revisits)){
      revisits$index2[j]= revisits$index3[j-1]+1
      revisits$index3[j]= revisits$index2[j]+(revisits$X1[j]-1)
    } 
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
  revisits2$Between_visit_duration = NA
  revisits2$Visit_no=seq(from = 1, to =nrow(revisits2), by = 1)
  
  ID = unique(tempvisit$id)
  Species = unique(tempvisit$species) 
  Day = unique(temp$day)
  Model = "BM"
  buffer = unique(temp$buffer)
  Visits = nrow(revisits2)
  
  res = data.frame(cbind(ID,Species,Day,Model,buffer,Visits))
  BM_visitdata[i]= list(res)
}

BM_visitdata = data.table::rbindlist(BM_visitdata)






datasplit2=split(OUF,as.factor(OUF$LoopID))
gc()

##Calculate diptryx encounter rates per day for OUF data
#Could have made this a function but just modifying the same loop as before

count=1
OUF_visitdata=c()
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
  
  if(length(temp$visit[which(temp$visit=="In")])==0){
    ID = unique(temp$ID)
    Species = unique(temp$Species) 
    Day = unique(temp$day)
    Model = "OUF"
    Visits = 0
    buffer = unique(temp$buffer)
    
    res = data.frame(cbind(ID,Species,Day,Model,buffer,Visits))
    OUF_visitdata[i]= list(res)
    next
  }
  
    if(length(temp$visit[which(temp$visit=="In")])==0){
    ID = unique(temp$ID)
    Species = unique(temp$Species) 
    Day = unique(temp$day)
    Model = "OUF"
    Visits = 0
    buffer = unique(temp$buffer)
    
    res = data.frame(cbind(ID,Species,Day,Model,buffer,Visits))
    OUF_visitdata[i]= list(res)
    next
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
  
  if(nrow(revisits)==1){
    ID = unique(temp$ID)
    Species = unique(temp$Species) 
    Day = unique(temp$day)
    Model = "OUF"
    Visits = 1
    buffer = unique(temp$buffer)
    
    res = data.frame(cbind(ID,Species,Day,Model,buffer,Visits))
    OUF_visitdata[i]= list(res)
    next
  }
  
  revisits2=c()
  if(nrow(revisits)>1){
    for(j in 2:nrow(revisits)){
      revisits$index2[j]= revisits$index3[j-1]+1
      revisits$index3[j]= revisits$index2[j]+(revisits$X1[j]-1)
    } 
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
  revisits2$Between_visit_duration = NA
  revisits2$Visit_no=seq(from = 1, to =nrow(revisits2), by = 1)
  
  ID = unique(tempvisit$id)
  Species = unique(tempvisit$species) 
  Day = unique(temp$day)
  Model = "OUF"
  Visits = nrow(revisits2)
  buffer = unique(temp$buffer)
  
  res = data.frame(cbind(ID,Species,Day,Model,buffer,Visits))
  OUF_visitdata[i]= list(res)
}

OUF_visitdata = data.table::rbindlist(OUF_visitdata)


data <- sp::SpatialPointsDataFrame(coords = data[,c(3,4)], data = data,
                                       proj4string=sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
data <- spTransform(data, sp::CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
datapoly=over(data,crowns)
datapoly10=over(data,crowns10)
datapoly50=over(data,crowns50)
datapoly100=over(data,crowns100)
datapoly150=over(data,crowns150)
datapoly200=over(data,crowns200)


#Calculate crown encounter rates for actual data accross buffers
data = data.frame(data)
data310=cbind(data,datapoly10)
data350=cbind(data,datapoly50)
data3100=cbind(data,datapoly100)
data3150=cbind(data,datapoly150)
data3200=cbind(data,datapoly200)
data3=cbind(data,datapoly)
data3$buffer = 0
data310$buffer = 10
data350$buffer = 50
data3100$buffer = 100
data3150$buffer = 150
data3200$buffer = 200

obs=c()
count=1
obs[count]=list(data3)

count=count+1
obs[count]=list(data310)

count=count+1
obs[count]=list(data350)

count=count+1
obs[count]=list(data3100)

count=count+1
obs[count]=list(data3150)

count=count+1
obs[count]=list(data3200)

obs = data.table::rbindlist(obs)

obs$LoopID = paste(obs$individual.local.identifier, obs$day, obs$buffer, sep = "_")
obs$visit=NA

datasplit2=split(obs,as.factor(obs$LoopID))
gc()

count=1
data_visitdata=c()
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
  
  if(length(temp$visit[which(temp$visit=="In")])==0){
    ID = unique(temp$individual.local.identifier)
    Species = unique(temp$individual.taxon.canonical.name) 
    Day = unique(temp$day)
    Model = "Observed"
    Visits = 0
    
    buffer = unique(temp$buffer)
    
    res = data.frame(cbind(ID,Species,Day,Model,buffer,Visits))
    data_visitdata[i]= list(res)
    next
  }
  
  if(length(temp$visit[which(temp$visit=="In")])==0){
    ID = unique(temp$individual.local.identifier)
    Species = unique(temp$individual.taxon.canonical.name) 
    Day = unique(temp$day)
    Model = "Observed"
    Visits = 0
    buffer = unique(temp$buffer)
    
    res = data.frame(cbind(ID,Species,Day,Model,buffer,Visits))
    data_visitdata[i]= list(res)
    next
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
  
  if(nrow(revisits)==1){
    ID = unique(temp$individual.local.identifier)
    Species = unique(temp$individual.taxon.canonical.name) 
    Day = unique(temp$day)
    Model = "Observed"
    Visits = 1
    buffer = unique(temp$buffer)
    
    res = data.frame(cbind(ID,Species,Day,Model,buffer,Visits))
    data_visitdata[i]= list(res)
    next
  }
  
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
  revisits2$Between_visit_duration = NA
  revisits2$Visit_no=seq(from = 1, to =nrow(revisits2), by = 1)
  
  ID = unique(temp$individual.local.identifier)
  Species = unique(temp$individual.taxon.canonical.name) 
  Day = unique(temp$day)
  Model = "Observed"
  Visits = nrow(revisits2)
  buffer = unique(temp$buffer)
  
  res = data.frame(cbind(ID,Species,Day,Model,buffer,Visits))
  data_visitdata[i]= list(res)
}



data_visitdata = data.table::rbindlist(data_visitdata)

#Merge simulated and real data
Alldata=rbind(data_visitdata,IID_visitdata,BM_visitdata,OUF_visitdata)
Alldata=Alldata[-which(Alldata$Species=="Pecari tajacu"),]
Alldata=Alldata[-which(Alldata$ID=="Atlas 4673"),]
Alldata=Alldata[-which(Alldata$ID=="Serge 4670"),]
Alldata$Visits=as.numeric(Alldata$Visits)
Alldata$Visits2=Alldata$Visits+1
Alldata$Model=factor(Alldata$Model, levels = c("Observed","IID","BM","OUF"))
Alldata$Species=as.factor(Alldata$Species)
Alldata$ID=as.factor(Alldata$ID)
Alldata$buffer=factor(Alldata$buffer, levels = c("0","10","50","100", "150","200"))
Alldata$buffer2=as.numeric(as.character(Alldata$buffer))
Alldata$Species_model=paste(Alldata$Species, Alldata$Model, sep = "_")
Alldata$Species_model=as.factor(Alldata$Species_model)
Alldata2=Alldata[which(as.character(Alldata$buffer)=="0"),]
Alldata3=Alldata2[-which(as.character(Alldata2$Model)=="IID"),]

library(brms)
library(cmdstanr)
library(emmeans)
library(ggplot2)

theme_set(theme_classic())

#Visualize encounter rates accross models and buffers
a=ggplot(data=Alldata)+geom_boxplot(aes(x = buffer, y = Visits, color = Species))+ facet_wrap(~Model, scales = "free")

ggplot(data=Alldata)+geom_violin(aes(x = buffer, y = Visits, color = Species))+ facet_wrap(~Model, scales = "free")

b=ggplot(data=Alldata)+geom_boxplot(aes(x = buffer, y = Visits, color = Model))+ facet_wrap(~Species, scales = "free")

ggplot(data=Alldata)+geom_violin(aes(x = buffer, y = Visits, color = Model))+ facet_wrap(~Species, scales = "free")

ggpubr::ggarrange(a,b,common.legend = FALSE, legend = "bottom")

#Statistical comparison of treatments
options(mc.cores = parallel::detectCores()) 

model = brm(bf(Visits2 ~ Species + Model + Species:Model +(1|ID)),
            data = Alldata3,
            family = Gamma(link = "log"),
            backend = "cmdstanr",
            threads = threading(2),
            refresh = 5,
            save_pars = save_pars(all = TRUE),
            control = list(max_treedepth = 10, adapt_delta = .8))
summary(model)
pp_check(model, ndraws=100)
Modelplot=plot(conditional_effects(model), plot = F)[[1]]
Modelplot

EMM=emmeans(model , specs = pairwise ~ Species:Model)
plot(EMM)


par(mfrow=c(2,2))
BMsplit=split(BM, as.factor(BM$ID))
IIDsplit=split(IID, as.factor(IID$ID))
OUFsplit=split(OUF, as.factor(OUF$ID))
datasplit=split(data, as.factor(data$individual.local.identifier))
plot(BMsplit[[6]]$z_utm, asp = 1, type = "l", main = "Bob BM")
#plot(IIDsplit[[6]]$z_utm, asp = 1, type = "l", main = "Bob IID")
plot(OUFsplit[[6]]$z_utm, asp = 1, type = "l", main = "Bob OUF")
plot(datasplit[[6]]$location.long.2,datasplit[[6]]$location.lat.2, asp = 1, type = "l", main = "Bob Observed")
