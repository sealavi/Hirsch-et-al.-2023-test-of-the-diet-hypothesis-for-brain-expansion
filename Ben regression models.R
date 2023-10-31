### daily test data
DailyFFT <- read.csv("Daily_stats_for_Ben_2023_updated.csv")
EQ <- read.csv("EQ.csv")
Trackback <- read.csv("Trackback.csv")
fruit <- read.csv("Ro_fruitfall.csv")
DailyFFT2h=read.csv("Daily_stats_for_Ben_2023_2_hour_updated.csv")
DailyFFT4h=read.csv("Daily_stats_for_Ben_2023_4_hour_updated.csv")

Ro_fruitfall=read.csv("Ro_fruitfall.csv")


library(dplyr)
library(tidyverse) #for data wrangling
library(ggplot2)
library(lubridate)
library(emmeans)   #for estimating marginal means
library(lme4)

### Calculate efficiency metrics
theme_set(theme_classic())

DailyFFT <- DailyFFT %>% mutate(SLeff = (active_minutes*60/Straight_line_DJL))
DailyFFT <- DailyFFT %>% mutate(Intereff = (active_minutes*60/Interpolated_DJL))
DailyFFT <- DailyFFT %>% mutate(raweff = (visit_duration_sum/Straight_line_DJL))
DailyFFT <- DailyFFT %>% mutate(numbereff = (Number_Dip_visits/Straight_line_DJL))
DailyFFT <- DailyFFT %>% mutate(numbereff12 = (Number_Dip_visits_above_12/Straight_line_DJL))
DailyFFT <- DailyFFT %>% mutate(numbereff_interp = (Number_Dip_visits/Interpolated_DJL))
DailyFFT <- DailyFFT %>% mutate(numbereff12_interp = (Number_Dip_visits_above_12/Interpolated_DJL))
DailyFFT <-left_join(DailyFFT,Trackback,by="Species")
DailyFFT <-left_join(DailyFFT,EQ,by="Species")
DailyFFT <-left_join(DailyFFT,Ro_fruitfall,by="Week_year")
DailyFFT <- DailyFFT %>% mutate(TBeff = ((active_minutes*60*Trackback)/Straight_line_DJL))
DailyFFT <- DailyFFT %>% mutate(TBeff_interp = ((active_minutes*60*Trackback)/Interpolated_DJL))
DailyFFT <- DailyFFT %>% filter(Species != "")
DailyFFT <- DailyFFT %>% mutate(DiptDensity = ((dip_in_range/Home_range_estimate)))
DailyFFT <- DailyFFT %>% filter(ID != "Sahti 4693")
DailyFFT <- DailyFFT %>% filter(ID != "Atlas 4673")
DailyFFT <- DailyFFT %>% filter(ID != "Chloe 4052")
DailyFFT <- DailyFFT %>% filter(week != "NA")


### Append fruit data to data frame
DailyFFT$fruit_count=NA
DailyFFT$fruit_mean=NA
DailyFFT$fruit_sd=NA
for(i in 1:nrow(fruit)){
  DailyFFT$fruit_count[which(DailyFFT$Species==fruit$Species[i] & DailyFFT$Week_year==fruit$Week_year[i])]=fruit$count[i]
  DailyFFT$fruit_mean[which(DailyFFT$Species==fruit$Species[i] & DailyFFT$Week_year==fruit$Week_year[i])]=fruit$mean[i]
  DailyFFT$fruit_sd[which(DailyFFT$Species==fruit$Species[i] & DailyFFT$Week_year==fruit$Week_year[i])]=fruit$sd[i]
  
}


## weeks with fruitfall data 4-13
DailyFFT2 <- DailyFFT %>% filter(week > 3)
DailyFFT2 <- DailyFFT2 %>% filter(week < 14)


options(mc.cores = parallel::detectCores()) ##parallelize the chains

library(lubridate)
DailyFFT2$date=as.POSIXct(as.character(DailyFFT2$date), format="%Y-%m-%d",origin="01-01-1900", tz="UTC")
DailyFFT2$date2=as.numeric(DailyFFT2$date)
DailyFFT2$ID=as.factor(DailyFFT2$ID)
DailyFFT2$Species=as.factor(DailyFFT2$Species)

#save(DailyFFT2, file = "DailyFFT2.RData")
#write.csv(DailyFFT2,'DailyFFT2.csv')

## 2hr dataset
DailyFFT2H <- DailyFFT2h %>% mutate(SLeff = (active_minutes*60/Straight_line_DJL_2H))
DailyFFT2H <- DailyFFT2H %>% mutate(Intereff = (active_minutes*60/Interpolated_DJL_2H))
DailyFFT2H <- DailyFFT2H %>% mutate(raweff = (visit_duration_sum/Straight_line_DJL_2H))
DailyFFT2H <- DailyFFT2H %>% mutate(numbereff = (Number_Dip_visits/Straight_line_DJL_2H))
DailyFFT2H <- DailyFFT2H %>% mutate(numbereff12 = (Number_Dip_visits_above_12/Straight_line_DJL_2H))
DailyFFT2H <- DailyFFT2H %>% mutate(numbereff_interp = (Number_Dip_visits/Interpolated_DJL_2H))
DailyFFT2H <- DailyFFT2H %>% mutate(numbereff12_interp = (Number_Dip_visits_above_12/Interpolated_DJL_2H))
DailyFFT2H <-left_join(DailyFFT2H,Trackback,by="Species")
DailyFFT2H <-left_join(DailyFFT2H,EQ,by="Species")
DailyFFT2H <-left_join(DailyFFT2H,Ro_fruitfall,by="Week_year")
DailyFFT2H <- DailyFFT2H %>% mutate(TBeff = ((active_minutes*60*Trackback)/Straight_line_DJL_2H))
DailyFFT2H <- DailyFFT2H %>% mutate(TBeff_interp = ((active_minutes*60*Trackback)/Interpolated_DJL_2H))
DailyFFT2H <- DailyFFT2H %>% filter(Species != "")
DailyFFT2H <- DailyFFT2H %>% mutate(DiptDensity = ((dip_in_range/Home_range_estimate)))
DailyFFT2H <- DailyFFT2H %>% filter(ID != "Sahti 4693")
DailyFFT2H <- DailyFFT2H %>% filter(ID != "Atlas 4673")
DailyFFT2H <- DailyFFT2H %>% filter(ID != "Chloe 4052")
DailyFFT2H <- DailyFFT2H %>% filter(week != "NA")



DailyFFT2H$fruit_count=NA
DailyFFT2H$fruit_mean=NA
DailyFFT2H$fruit_sd=NA
for(i in 1:nrow(fruit)){
  DailyFFT2H$fruit_count[which(DailyFFT2H$Species==fruit$Species[i] & DailyFFT2H$Week_year==fruit$Week_year[i])]=fruit$count[i]
  DailyFFT2H$fruit_mean[which(DailyFFT2H$Species==fruit$Species[i] & DailyFFT2H$Week_year==fruit$Week_year[i])]=fruit$mean[i]
  DailyFFT2H$fruit_sd[which(DailyFFT2H$Species==fruit$Species[i] & DailyFFT2H$Week_year==fruit$Week_year[i])]=fruit$sd[i] }

## weeks with fruitfall data 4-13
DailyFFT2H <- DailyFFT2H %>% filter(week > 3)
DailyFFT2H <- DailyFFT2H %>% filter(week < 14)

DailyFFT2H$date=as.POSIXct(as.character(DailyFFT2H$date), format="%Y-%m-%d",origin="01-01-1900", tz="UTC")
DailyFFT2H$date2=as.numeric(DailyFFT2H$date)
DailyFFT2H$ID=as.factor(DailyFFT2H$ID)
DailyFFT2H$Species=as.factor(DailyFFT2H$Species)

#save(DailyFFT2H, file = "DailyFFT2H.RData")
#write.csv(DailyFFT2H,'DailyFFT2H.csv')

## 4hr dataset
DailyFFT4H  <- DailyFFT4h  %>% mutate(SLeff = (active_minutes*60/Straight_line_DJL_4H))
DailyFFT4H  <- DailyFFT4H  %>% mutate(Intereff = (active_minutes*60/Interpolated_DJL_4H))
DailyFFT4H  <- DailyFFT4H  %>% mutate(raweff = (visit_duration_sum/Straight_line_DJL_4H))
DailyFFT4H  <- DailyFFT4H  %>% mutate(numbereff = (Number_Dip_visits/Straight_line_DJL_4H))
DailyFFT4H  <- DailyFFT4H  %>% mutate(numbereff12 = (Number_Dip_visits_above_12/Straight_line_DJL_4H))
DailyFFT4H  <- DailyFFT4H  %>% mutate(numbereff_interp = (Number_Dip_visits/Interpolated_DJL_4H))
DailyFFT4H  <- DailyFFT4H  %>% mutate(numbereff12_interp = (Number_Dip_visits_above_12/Interpolated_DJL_4H))
DailyFFT4H  <-left_join(DailyFFT4H ,Trackback,by="Species")
DailyFFT4H  <-left_join(DailyFFT4H ,EQ,by="Species")
DailyFFT4H  <-left_join(DailyFFT4H ,Ro_fruitfall,by="Week_year")
DailyFFT4H  <- DailyFFT4H  %>% mutate(TBeff = ((active_minutes*60*Trackback)/Straight_line_DJL_4H))
DailyFFT4H  <- DailyFFT4H  %>% mutate(TBeff_interp = ((active_minutes*60*Trackback)/Interpolated_DJL_4H))
DailyFFT4H  <- DailyFFT4H  %>% filter(Species != "")
DailyFFT4H  <- DailyFFT4H  %>% mutate(DiptDensity = ((dip_in_range/Home_range_estimate)))
DailyFFT4H  <- DailyFFT4H  %>% filter(ID != "Sahti 4693")
DailyFFT4H  <- DailyFFT4H  %>% filter(ID != "Atlas 4673")
DailyFFT4H  <- DailyFFT4H  %>% filter(ID != "Chloe 4052")
DailyFFT4H  <- DailyFFT4H  %>% filter(week != "NA")


DailyFFT4H$fruit_count=NA
DailyFFT4H$fruit_mean=NA
DailyFFT4H$fruit_sd=NA
for(i in 1:nrow(fruit)){
  DailyFFT4H$fruit_count[which(DailyFFT4H$Species==fruit$Species[i] & DailyFFT4H$Week_year==fruit$Week_year[i])]=fruit$count[i]
  DailyFFT4H$fruit_mean[which(DailyFFT4H$Species==fruit$Species[i] & DailyFFT4H$Week_year==fruit$Week_year[i])]=fruit$mean[i]
  DailyFFT4H$fruit_sd[which(DailyFFT4H$Species==fruit$Species[i] & DailyFFT4H$Week_year==fruit$Week_year[i])]=fruit$sd[i] }

## weeks with fruitfall data 4-13
DailyFFT4H <- DailyFFT4H %>% filter(week > 3)
DailyFFT4H <- DailyFFT4H %>% filter(week < 14)

DailyFFT4H$date=as.POSIXct(as.character(DailyFFT4H$date), format="%Y-%m-%d",origin="01-01-1900", tz="UTC")
DailyFFT4H$date2=as.numeric(DailyFFT4H$date)
DailyFFT4H$ID=as.factor(DailyFFT4H$ID)
DailyFFT4H$Species=as.factor(DailyFFT4H$Species)

#save(DailyFFT4H, file = "DailyFFT4H.RData")
#write.csv(DailyFFT4H,'DailyFFT4H.csv')

#remove NAs and convert durations to proper units
DailyFFT2  <- DailyFFT2  %>% filter(TBeff != "NA")
DailyFFT2  <- DailyFFT2  %>% mutate(DipAct = (active_minutes*60))

DailyFFT2H  <- DailyFFT2H  %>% filter(TBeff != "NA")
DailyFFT2H  <- DailyFFT2H  %>% mutate(DipAct = (active_minutes*60))

DailyFFT4H  <- DailyFFT4H  %>% filter(TBeff != "NA")
DailyFFT4H  <- DailyFFT4H  %>% mutate(DipAct = (active_minutes*60))

##Remove data where two taus couldn't be estimated
DailyFFT2=DailyFFT2[which(DailyFFT2$tau==2),]
DailyFFT2H=DailyFFT2H[which(DailyFFT2H$tau==2),]
DailyFFT4H=DailyFFT4H[which(DailyFFT4H$tau==2),]

#remove INF values
DailyFFT2=DailyFFT2[-which(is.infinite(DailyFFT2$ctmm.estimate)==TRUE),]
DailyFFT2H=DailyFFT2H[-which(is.infinite(DailyFFT2H$ctmm.estimate)==TRUE),]
DailyFFT4H=DailyFFT4H[-which(is.infinite(DailyFFT4H$ctmm.estimate)==TRUE),]

##sanity checking
everyday=unique(paste(DailyFFT2$ID,DailyFFT2$Day,sep=" ,day "))
everyday_2h=unique(paste(DailyFFT2H$ID,DailyFFT2H$Day,sep=" ,day "))
everyday_4h=unique(paste(DailyFFT4H$ID,DailyFFT4H$Day,sep=" ,day "))

remove2H=setdiff(everyday_2h, everyday)
remove4H=setdiff(everyday_4h, everyday)
remove2Hb=setdiff(everyday, everyday_2h)
remove4Hb=setdiff(everyday, everyday_4h)

#Ensure that datasets are comparable
for(i in 1:length(remove2H)){
  if(length(which(paste(DailyFFT2H$ID,DailyFFT2H$Day,sep=" ,day ")==remove2H[i]))==0){next}
  DailyFFT2H=DailyFFT2H[-which(paste(DailyFFT2H$ID,DailyFFT2H$Day,sep=" ,day ")==remove2H[i]),]
}
for(i in 1:length(remove4H)){
  if(length(which(paste(DailyFFT4H$ID,DailyFFT4H$Day,sep=" ,day ")==remove4H[i]))==0){next}
  DailyFFT4H=DailyFFT4H[-which(paste(DailyFFT4H$ID,DailyFFT4H$Day,sep=" ,day ")==remove4H[i]),]
}
for(i in 1:length(remove2Hb)){
  if(length(which(paste(DailyFFT2$ID,DailyFFT2$Day,sep=" ,day ")==remove2Hb[i]))==0){next}
  DailyFFT2=DailyFFT2[-which(paste(DailyFFT2$ID,DailyFFT2$Day,sep=" ,day ")==remove2Hb[i]),]
}
for(i in 1:length(remove4Hb)){
  if(length(which(paste(DailyFFT2$ID,DailyFFT2$Day,sep=" ,day ")==remove4Hb[i]))==0){next}
  DailyFFT2=DailyFFT2[-which(paste(DailyFFT2$ID,DailyFFT2$Day,sep=" ,day ")==remove4Hb[i]),]
}

everyday=unique(paste(DailyFFT2$ID,DailyFFT2$Day,sep=" ,day "))
everyday_2h=unique(paste(DailyFFT2H$ID,DailyFFT2H$Day,sep=" ,day "))
everyday_4h=unique(paste(DailyFFT4H$ID,DailyFFT4H$Day,sep=" ,day "))

remove2H=setdiff(everyday_2h, everyday)
remove4H=setdiff(everyday_4h, everyday)
remove2Hb=setdiff(everyday, everyday_2h)
remove4Hb=setdiff(everyday, everyday_4h)


#Add EQ and taxonomic data
DailyFFT2$Clade=NA
DailyFFT2H$Clade=NA
DailyFFT4H$Clade=NA

DailyFFT2$EQ=NA
DailyFFT2H$EQ=NA
DailyFFT4H$EQ=NA


DailyFFT2$Clade[which(DailyFFT2$Species==as.character("Ateles geoffroyi"))]="Primate"
DailyFFT2$Clade[which(DailyFFT2$Species==as.character("Cebus capucinus"))]="Primate"
DailyFFT2$Clade[which(DailyFFT2$Species==as.character("Nasua narica"))]="Procyonid"
DailyFFT2$Clade[which(DailyFFT2$Species==as.character("Potos flavus"))]="Procyonid"

DailyFFT2H$Clade[which(DailyFFT2H$Species==as.character("Ateles geoffroyi"))]="Primate"
DailyFFT2H$Clade[which(DailyFFT2H$Species==as.character("Cebus capucinus"))]="Primate"
DailyFFT2H$Clade[which(DailyFFT2H$Species==as.character("Nasua narica"))]="Procyonid"
DailyFFT2H$Clade[which(DailyFFT2H$Species==as.character("Potos flavus"))]="Procyonid"

DailyFFT4H$Clade[which(DailyFFT4H$Species==as.character("Ateles geoffroyi"))]="Primate"
DailyFFT4H$Clade[which(DailyFFT4H$Species==as.character("Cebus capucinus"))]="Primate"
DailyFFT4H$Clade[which(DailyFFT4H$Species==as.character("Nasua narica"))]="Procyonid"
DailyFFT4H$Clade[which(DailyFFT4H$Species==as.character("Potos flavus"))]="Procyonid"

DailyFFT2$EQ[which(DailyFFT2$Species==as.character("Ateles geoffroyi"))]=EQ$EQ[which(EQ$Species=="Ateles geoffroyi")]
DailyFFT2$EQ[which(DailyFFT2$Species==as.character("Cebus capucinus"))]=EQ$EQ[which(EQ$Species=="Cebus capucinus")]
DailyFFT2$EQ[which(DailyFFT2$Species==as.character("Nasua narica"))]=EQ$EQ[which(EQ$Species=="Nasua narica")]
DailyFFT2$EQ[which(DailyFFT2$Species==as.character("Potos flavus"))]=EQ$EQ[which(EQ$Species=="Potos flavus")]

DailyFFT2H$EQ[which(DailyFFT2H$Species==as.character("Ateles geoffroyi"))]=EQ$EQ[which(EQ$Species=="Ateles geoffroyi")]
DailyFFT2H$EQ[which(DailyFFT2H$Species==as.character("Cebus capucinus"))]=EQ$EQ[which(EQ$Species=="Cebus capucinus")]
DailyFFT2H$EQ[which(DailyFFT2H$Species==as.character("Nasua narica"))]=EQ$EQ[which(EQ$Species=="Nasua narica")]
DailyFFT2H$EQ[which(DailyFFT2H$Species==as.character("Potos flavus"))]=EQ$EQ[which(EQ$Species=="Potos flavus")]

DailyFFT4H$EQ[which(DailyFFT4H$Species==as.character("Ateles geoffroyi"))]=EQ$EQ[which(EQ$Species=="Ateles geoffroyi")]
DailyFFT4H$EQ[which(DailyFFT4H$Species==as.character("Cebus capucinus"))]=EQ$EQ[which(EQ$Species=="Cebus capucinus")]
DailyFFT4H$EQ[which(DailyFFT4H$Species==as.character("Nasua narica"))]=EQ$EQ[which(EQ$Species=="Nasua narica")]
DailyFFT4H$EQ[which(DailyFFT4H$Species==as.character("Potos flavus"))]=EQ$EQ[which(EQ$Species=="Potos flavus")]


##Build regression models
library(brms)
library(cmdstanr)
options(mc.cores = parallel::detectCores()) ##parallelize the chains

### New models
## 0- Dipt Density, Species, and species*sex interaction
## 1- Dipt Density, Species
## 2- Dipt Density, Clade 
## 3- Dipt Density, EQ

## Model class names:
## All day 
## Models for route efficiency = raweff = A models
## Models for track back corrected route efficiency = TBeff = B models
## 2hr
## Models for route efficiency = raweff = 2A models
## Models for track back corrected route efficiency = TBeff = 2B models
## 4hr
## Models for route efficiency = raweff = 4A models
## Models for track back corrected route efficiency = TBeff = 4B models

## class A - all day raweff models

#Scale continuous covariates
DailyFFT2$DiptDensity_scaled=(DailyFFT2$DiptDensity-mean(DailyFFT2$DiptDensity,na.rm = TRUE))/sd(DailyFFT2$DiptDensity,na.rm = TRUE)
DailyFFT2$Ro_fruitfall_scaled=(DailyFFT2$DiptDensity-mean(DailyFFT2$Ro_fruitfall,na.rm = TRUE))/sd(DailyFFT2$Ro_fruitfall,na.rm = TRUE)
DailyFFT2$week_scaled=(DailyFFT2$week-mean(DailyFFT2$week,na.rm = TRUE))/sd(DailyFFT2$week,na.rm = TRUE)
DailyFFT2$EQ_scaled=(DailyFFT2$EQ-mean(DailyFFT2$EQ,na.rm = TRUE))/sd(DailyFFT2$EQ,na.rm = TRUE)

DailyFFT2H$DiptDensity_scaled=(DailyFFT2H$DiptDensity-mean(DailyFFT2H$DiptDensity,na.rm = TRUE))/sd(DailyFFT2H$DiptDensity,na.rm = TRUE)
DailyFFT2H$Ro_fruitfall_scaled=(DailyFFT2H$DiptDensity-mean(DailyFFT2H$Ro_fruitfall,na.rm = TRUE))/sd(DailyFFT2H$Ro_fruitfall,na.rm = TRUE)
DailyFFT2H$week_scaled=(DailyFFT2H$week-mean(DailyFFT2H$week,na.rm = TRUE))/sd(DailyFFT2H$week,na.rm = TRUE)
DailyFFT2H$EQ_scaled=(DailyFFT2H$EQ-mean(DailyFFT2H$EQ,na.rm = TRUE))/sd(DailyFFT2H$EQ,na.rm = TRUE)

DailyFFT4H$DiptDensity_scaled=(DailyFFT4H$DiptDensity-mean(DailyFFT4H$DiptDensity,na.rm = TRUE))/sd(DailyFFT4H$DiptDensity,na.rm = TRUE)
DailyFFT4H$Ro_fruitfall_scaled=(DailyFFT4H$DiptDensity-mean(DailyFFT4H$Ro_fruitfall,na.rm = TRUE))/sd(DailyFFT4H$Ro_fruitfall,na.rm = TRUE)
DailyFFT4H$week_scaled=(DailyFFT4H$week-mean(DailyFFT4H$week,na.rm = TRUE))/sd(DailyFFT4H$week,na.rm = TRUE)
DailyFFT4H$EQ_scaled=(DailyFFT4H$EQ-mean(DailyFFT4H$EQ,na.rm = TRUE))/sd(DailyFFT4H$EQ,na.rm = TRUE)

source("https://raw.githubusercontent.com/sealavi/GAM-first-derivative-functions/main/gam%20first%20derivitave%20plot%20function_general.R")

## class A - all day raweff models
nonlinear_modelA0=brm(bf((Intereff)~s(DiptDensity_scaled)+Species+Sex:Species+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2,
                      family = Gamma(link="log"),
                      init=0,
                      iter = 6000,
                      chains=2,
                      prior=c(prior(normal(0,1.45),class="Intercept"),
                              prior(normal(0,1.45),class="b"),
                              prior(gamma(.01,.01),class="shape"),
                              prior(normal(0,1.45),class="sd")),
                      save_pars = save_pars(all = TRUE),
                      backend = "cmdstanr",
                      threads = threading(2),
                      control=list(adapt_delta=0.9999, max_treedepth = 10))

save(nonlinear_modelA0, file = "nonlinear_modelA0.RData")
summary(nonlinear_modelA0)

#precompute metrics for model comparison and evaluation 
nonlinear_modelA0 = add_criterion(nonlinear_modelA0, c("waic","bayes_R2","loo", "loo_R2"), 
                                  control=list(adapt_delta=0.9999, max_treedepth = 10), 
                       backend = "cmdstanr",
                       reloo = TRUE)
#estimate derivative of splines
deriv_plot(model=nonlinear_modelA0,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelA0_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_modelA0,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelA0_DiptDensity_scaled_plot2")
#Pairwise comparisons 
EMM=emmeans(nonlinear_modelA0 , specs = pairwise ~ Species)
EMM
plot(EMM)
EMM2=emmeans(nonlinear_modelA0 , specs = pairwise ~ Sex:Species)
EMM2
plot(EMM2)
save(nonlinear_modelA0, file = "nonlinear_modelA0.RData")

loo(nonlinear_modelA0)
loo_R2(nonlinear_modelA0)

nonlinear_modelA1=brm(bf((Intereff)~s(DiptDensity_scaled)+Species+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2,
                      family = Gamma(link="log"),
                      init=0,
                      iter = 3000,
                      chains=2,
                      prior=c(prior(normal(0,1.45),class="Intercept"),
                              prior(normal(0,1.45),class="b"),
                              prior(gamma(.01,.01),class="shape"),
                              prior(normal(0,1.45),class="sd")),
                      save_pars = save_pars(all = TRUE),
                      backend = "cmdstanr",
                      threads = threading(2),
                      control=list(adapt_delta=0.9999, max_treedepth = 10))

save(nonlinear_modelA1, file="nonlinear_modelA1.RData")
summary(nonlinear_modelA1)
plot(conditional_effects(nonlinear_modelA1,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_modelA1)
nonlinear_modelA1 = add_criterion(nonlinear_modelA1, c("waic","bayes_R2","loo", "loo_R2"), 
                                  control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                  backend = "cmdstanr",
                                  reloo = TRUE)
loo_R2(nonlinear_modelA1)

deriv_plot(model=nonlinear_modelA1,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelA1_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_modelA1,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelA1_DiptDensity_scaled_plot2")
EMMa=emmeans(nonlinear_modelA1 , specs = pairwise ~ Species)
EMMa
EMMplotA =plot(EMMa)
save(nonlinear_modelA1, file = "nonlinear_modelA1.RData")


nonlinear_modelA2=brm(bf((Intereff)~s(DiptDensity_scaled)+Clade+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2,
                      family = Gamma(link="log"),
                      init=0,iter = 5000,
                      chains=2,
                      prior=c(prior(normal(0,1.45),class="Intercept"),
                              prior(normal(0,1.45),class="b"),
                              prior(gamma(.01,.01),class="shape"),
                              prior(normal(0,1.45),class="sd")),
                      save_pars = save_pars(all = TRUE),
                      backend = "cmdstanr",
                      threads = threading(2),
                      control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_modelA2, file = "nonlinear_modelA2.RData")
summary(nonlinear_modelA2)
plot(conditional_effects(nonlinear_modelA2,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_modelA2)
nonlinear_modelA2 = add_criterion(nonlinear_modelA2, c("waic","bayes_R2","loo", "loo_R2"), 
                                  control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                  backend = "cmdstanr",
                                  reloo = TRUE)
loo_R2(nonlinear_modelA2)
deriv_plot(model=nonlinear_modelA2,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelA2_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_modelA2,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelA2_DiptDensity_scaled_plot2")
EMMb=emmeans(nonlinear_modelA2 , specs = pairwise ~ Clade)
EMMb
plot(EMMb)
save(nonlinear_modelA2, file = "nonlinear_modelA2.RData")

nonlinear_modelA3=brm(bf((Intereff)~s(DiptDensity_scaled)+s(EQ, k=4)+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2,
                      family = Gamma(link="log"),
                      init=0,
                      iter = 3000,
                      chains=2,
                      prior=c(prior(normal(0,1.45),class="Intercept"),
                              prior(normal(0,1.45),class="b"),
                              prior(gamma(.01,.01),class="shape"),
                              prior(normal(0,1.45),class="sd")),
                      save_pars = save_pars(all = TRUE),
                      backend = "cmdstanr",
                      threads = threading(2),
                      control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_modelA3, file = "nonlinear_modelA3.RData")
summary(nonlinear_modelA3)
plot(conditional_effects(nonlinear_modelA3,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_modelA3)
nonlinear_modelA3 = add_criterion(nonlinear_modelA3, c("waic","bayes_R2","loo", "loo_R2"), 
                                  control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                  backend = "cmdstanr",
                                  reloo = TRUE)
loo_R2(nonlinear_modelA3)
deriv_plot(model=nonlinear_modelA3,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelA3_DiptDensity_scaled_plot")
deriv_plot(model=nonlinear_modelA3,term="s(EQ, k=4)",main="EQ",eps=0.001,output="nonlinear_modelA3_EQ_plot")
save(nonlinear_modelA3, file = "nonlinear_modelA3.RData")
deriv_plot2(model=nonlinear_modelA3,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelA3_DiptDensity_scaled_plot2")
deriv_plot2(model=nonlinear_modelA3,term="s(EQ, k=4)",main="EQ",eps=0.001,output="nonlinear_modelA3_EQ_plot2")
save(nonlinear_modelA3, file = "nonlinear_modelA3.RData")



## Class B - all day TBeff models

nonlinear_modelB0=brm(bf((TBeff_interp)~s(DiptDensity_scaled)+Species+Sex:Species+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2,
                      family = Gamma(link="log"),
                      init=0,
                      iter = 8000,
                      chains=2,
                      prior=c(prior(normal(0,1.45),class="Intercept"),
                              prior(normal(0,1.45),class="b"),
                              prior(gamma(.01,.01),class="shape"),
                              prior(normal(0,1.45),class="sd")),
                      save_pars = save_pars(all = TRUE),
                      backend = "cmdstanr",
                      threads = threading(2),
                      control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_modelB0, file = "nonlinear_modelB0.RData")
summary(nonlinear_modelB0)
plot(conditional_effects(nonlinear_modelB0,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_modelB0)
nonlinear_modelB0 = add_criterion(nonlinear_modelB0, c("waic","bayes_R2","loo", "loo_R2"), 
                                  control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                  backend = "cmdstanr",
                                  reloo = TRUE)
loo_R2(nonlinear_modelB0)
deriv_plot(model=nonlinear_modelB0,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelB0_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_modelB0,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelB0_DiptDensity_scaled_plot2")
EMMd=emmeans(nonlinear_modelB0 , specs = pairwise ~ Species)
EMMd
plot(EMMd)
save(nonlinear_modelB0, file = "nonlinear_modelB0.RData")



nonlinear_modelB1=brm(bf((TBeff_interp)~s(DiptDensity_scaled)+Species+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2,
                      family = Gamma(link="log"),
                      init=0,
                      iter = 3000,
                      chains=2,
                      prior=c(prior(normal(0,1.45),class="Intercept"),
                              prior(normal(0,1.45),class="b"),
                              prior(gamma(.01,.01),class="shape"),
                              prior(normal(0,1.45),class="sd")),
                      save_pars = save_pars(all = TRUE),
                      backend = "cmdstanr",
                      threads = threading(2),
                      control=list(adapt_delta=0.9999, max_treedepth = 10))

save(nonlinear_modelB1, file = "nonlinear_modelB1.RData")
summary(nonlinear_modelB1)
plot(conditional_effects(nonlinear_modelB1,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_modelB1)
nonlinear_modelB1 = add_criterion(nonlinear_modelB1, c("waic","bayes_R2","loo", "loo_R2"), 
                                  control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                  backend = "cmdstanr",
                                  reloo = TRUE)
loo_R2(nonlinear_modelB1)
deriv_plot(model=nonlinear_modelB1,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelB1_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_modelB1,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelB1_DiptDensity_scaled_plot2")
EMMd=emmeans(nonlinear_modelB1 , specs = pairwise ~ Species)
EMMd
EMMplotB = plot(EMMd)
save(nonlinear_modelB1, file = "nonlinear_modelB1.RData")



nonlinear_modelB2=brm(bf((TBeff_interp)~s(DiptDensity_scaled)+Clade+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2,
                      family = Gamma(link="log"),
                      init=0,iter = 5000,
                      chains=2,
                      prior=c(prior(normal(0,1.45),class="Intercept"),
                              prior(normal(0,1.45),class="b"),
                              prior(gamma(.01,.01),class="shape"),
                              prior(normal(0,1.45),class="sd")),
                      save_pars = save_pars(all = TRUE),
                      backend = "cmdstanr",
                      threads = threading(2),
                      control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_modelB2, file = "nonlinear_modelB2.RData")
summary(nonlinear_modelB2)
plot(conditional_effects(nonlinear_modelB2,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_modelB2)
nonlinear_modelB2 = add_criterion(nonlinear_modelB2, c("waic","bayes_R2","loo", "loo_R2"), 
                                  control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                  backend = "cmdstanr",
                                  reloo = TRUE)
loo_R2(nonlinear_modelB2)
deriv_plot(model=nonlinear_modelB2,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelB2_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_modelB2,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelB2_DiptDensity_scaled_plot2")
EMMe=emmeans(nonlinear_modelB2 , specs = pairwise ~ Clade)
EMMe
plot(EMMe)
save(nonlinear_modelB2, file = "nonlinear_modelB2.RData")

nonlinear_modelB3=brm(bf((TBeff_interp)~s(DiptDensity_scaled)+s(EQ, k=4)+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2,
                      family = Gamma(link="log"),
                      init=0,
                      iter = 3000,
                      chains=2,
                      prior=c(prior(normal(0,1.45),class="Intercept"),
                              prior(normal(0,1.45),class="b"),
                              prior(gamma(.01,.01),class="shape"),
                              prior(normal(0,1.45),class="sd")),
                      save_pars = save_pars(all = TRUE),
                      backend = "cmdstanr",
                      threads = threading(2),
                      control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_modelB3, file = "nonlinear_modelB3.RData")
summary(nonlinear_modelB3)
plot(conditional_effects(nonlinear_modelB3,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_modelB3)
nonlinear_modelB3 = add_criterion(nonlinear_modelB3, c("waic","bayes_R2","loo", "loo_R2"), 
                                  control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                  backend = "cmdstanr",
                                  reloo = TRUE)
loo_R2(nonlinear_modelB3)
deriv_plot(model=nonlinear_modelB3,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelB3_DiptDensity_scaled_plot")
deriv_plot(model=nonlinear_modelB3,term="s(EQ, k=4)",main="EQ",eps=0.001,output="nonlinear_modelB3_EQ_plot")
deriv_plot2(model=nonlinear_modelB3,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_modelB3_DiptDensity_scaled_plot2")
deriv_plot2(model=nonlinear_modelB3,term="s(EQ, k=4)",main="EQ",eps=0.001,output="nonlinear_modelB3_EQ_plot2")
save(nonlinear_modelB3, file = "nonlinear_modelB3.RData")



## class 2A - 2hr raweff models
nonlinear_model2A0=brm(bf((Intereff+.00001)~s(DiptDensity_scaled)+Species+Sex:Species+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 6000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model2A0, file = "nonlinear_model2A0.RData")
summary(nonlinear_model2A0)
plot(conditional_effects(nonlinear_model2A0,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model2A0)
nonlinear_model2A0 = add_criterion(nonlinear_model2A0, c("waic","bayes_R2","loo", "loo_R2"), 
                                  control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                  backend = "cmdstanr",
                                  reloo = TRUE)
loo_R2(nonlinear_model2A0)
deriv_plot(model=nonlinear_model2A0,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2A0_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_model2A0,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2A0_DiptDensity_scaled_plot2")
EMMf=emmeans(nonlinear_model2A0 , specs = pairwise ~ Species)
EMMf
plot(EMMf)
EMMf2=emmeans(nonlinear_model2A0 , specs = pairwise ~ Sex:Species)
EMMf2
plot(EMMf2)
save(nonlinear_model2A0, file = "nonlinear_model2A0.RData")

nonlinear_model2A1=brm(bf((Intereff+.00001)~s(DiptDensity_scaled)+Species+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 3000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model2A1, file = "nonlinear_model2A1.RData")
summary(nonlinear_model2A1)
plot(conditional_effects(nonlinear_model2A1,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model2A1)
nonlinear_model2A1 = add_criterion(nonlinear_model2A1, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model2A1)
deriv_plot(model=nonlinear_model2A1,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2A1_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_model2A1,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2A1_DiptDensity_scaled_plot2")
EMMg=emmeans(nonlinear_model2A1 , specs = pairwise ~ Species)
EMMg
EMMplotE = plot(EMMg)
save(nonlinear_model2A1, file = "nonlinear_model2A1.RData")

nonlinear_model2A2=brm(bf((Intereff+.00001)~s(DiptDensity_scaled)+Clade+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 5000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model2A2, file = "nonlinear_model2A2.RData")
summary(nonlinear_model2A2)
plot(conditional_effects(nonlinear_model2A2,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model2A2)
nonlinear_model2A2 = add_criterion(nonlinear_model2A2, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model2A2)
deriv_plot(model=nonlinear_model2A2,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2A2_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_model2A2,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2A2_DiptDensity_scaled_plot2")
EMMh=emmeans(nonlinear_model2A2 , specs = pairwise ~ Clade)
EMMh
plot(EMMh)
save(nonlinear_model2A2, file = "nonlinear_model2A2.RData")

nonlinear_model2A3=brm(bf((Intereff+.00001)~s(DiptDensity_scaled)+s(EQ, k=4)+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 3000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model2A3, file = "nonlinear_model2A3.RData")
summary(nonlinear_model2A3)
plot(conditional_effects(nonlinear_model2A3,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model2A3)
nonlinear_model2A3 = add_criterion(nonlinear_model2A3, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model2A3)
deriv_plot(model=nonlinear_model2A3,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2A3_DiptDensity_scaled_plot")
deriv_plot(model=nonlinear_model2A3,term="s(EQ, k=4)",main="EQ",eps=0.001,output="nnonlinear_model2A3_EQ_plot")
deriv_plot2(model=nonlinear_model2A3,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2A3_DiptDensity_scaled_plot2")
deriv_plot2(model=nonlinear_model2A3,term="s(EQ, k=4)",main="EQ",eps=0.001,output="nnonlinear_model2A3_EQ_plot2")
save(nonlinear_model2A3, file = "nonlinear_model2A3.RData")


## Class 2B - 2hr TBeff models

nonlinear_model2B0=brm(bf((TBeff_interp+.00001)~s(DiptDensity_scaled)+Species+Sex:Species+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 3000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model2B0, file = "nonlinear_model2B0.RData")
summary(nonlinear_model2B0)
plot(conditional_effects(nonlinear_model2B0,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model2B0)
nonlinear_model2B0 = add_criterion(nonlinear_model2B0, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model2B0)
deriv_plot(model=nonlinear_model2B0,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2B0_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_model2B0,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2B0_DiptDensity_scaled_plot2")
EMMi=emmeans(nonlinear_model2B0 , specs = pairwise ~ Species)
EMMi
plot(EMMi)
EMMi2=emmeans(nonlinear_model2B0 , specs = pairwise ~ Sex:Species)
EMMi2
EMMplotF= plot(EMMi2)
save(nonlinear_model2B0, file = "nonlinear_model2B0.RData")



nonlinear_model2B1=brm(bf((TBeff_interp+.00001)~s(DiptDensity_scaled)+Species+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 3000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model2B1, file = "nonlinear_model2B1.RData")
summary(nonlinear_model2B1)
plot(conditional_effects(nonlinear_model2B1,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model2B1)
nonlinear_model2B1 = add_criterion(nonlinear_model2B1, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model2B1)
deriv_plot(model=nonlinear_model2B1,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2B1_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_model2B1,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2B1_DiptDensity_scaled_plot2")
EMMj=emmeans(nonlinear_model2B1 , specs = pairwise ~ Species)
EMMj
plot(EMMj)
save(nonlinear_model2B1, file = "nonlinear_model2B1.RData")

nonlinear_model2B2=brm(bf((TBeff_interp+.00001)~s(DiptDensity_scaled)+Clade+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 5000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model2B2, file = "nonlinear_model2B2.RData")
summary(nonlinear_model2B2)
plot(conditional_effects(nonlinear_model2B2,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model2B2)
nonlinear_model2B2 = add_criterion(nonlinear_model2B2, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model2B2)
deriv_plot(model=nonlinear_model2B2,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2B2_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_model2B2,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2B2_DiptDensity_scaled_plot2")
EMMk=emmeans(nonlinear_model2B2 , specs = pairwise ~ Clade)
EMMk
plot(EMMk)
save(nonlinear_model2B2, file = "nonlinear_model2B2.RData")

nonlinear_model2B3=brm(bf((TBeff_interp+.00001)~s(DiptDensity_scaled)+s(EQ, k=4)+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT2H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 3000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model2B3, file = "nonlinear_model2B3.RData")
summary(nonlinear_model2B3)
plot(conditional_effects(nonlinear_model2B3,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model2B3)
nonlinear_model2B3 = add_criterion(nonlinear_model2B3, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model2B3)
deriv_plot(model=nonlinear_model2B3,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2B3_DiptDensity_scaled_plot")
deriv_plot(model=nonlinear_model2B3,term="s(EQ, k=4)",main="EQ",eps=0.001,output="nonlinear_model2B3_EQ_plot")
deriv_plot2(model=nonlinear_model2B3,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model2B3_DiptDensity_scaled_plot2")
deriv_plot2(model=nonlinear_model2B3,term="s(EQ, k=4)",main="EQ",eps=0.001,output="nonlinear_model2B3_EQ_plot2")
save(nonlinear_model2B3, file = "nonlinear_model2B3.RData")



## class 4A - 4hr raweff models
nonlinear_model3A0=brm(bf((Intereff+.0001)~s(DiptDensity_scaled)+Species+Sex:Species+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT4H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 6000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model3A0, file = "nonlinear_model3A0.RData")
summary(nonlinear_model3A0)
plot(conditional_effects(nonlinear_model3A0,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model3A0)
nonlinear_model3A0 = add_criterion(nonlinear_model3A0, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model3A0)
deriv_plot(model=nonlinear_model3A0,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model3A0_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_model3A0,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model3A0_DiptDensity_scaled_plot2")
EMMl=emmeans(nonlinear_model3A0 , specs = pairwise ~ Species)
EMMl
plot(EMMl)
EMMl2=emmeans(nonlinear_model3A0 , specs = pairwise ~ Sex:Species)
EMMl2
plot(EMMl2)
save(nonlinear_model3A0, file = "nonlinear_model3A0.RData")

nonlinear_model3A1=brm(bf((Intereff+.0001)~s(DiptDensity_scaled)+Species+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT4H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 3000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model3A1, file = "nonlinear_model3A1.RData")
summary(nonlinear_model3A1)
plot(conditional_effects(nonlinear_model3A1,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model3A1)
nonlinear_model3A1 = add_criterion(nonlinear_model3A1, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model3A1)
deriv_plot(model=nonlinear_model3A1,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model3A1_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_model3A1,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model3A1_DiptDensity_scaled_plot2")
EMMm=emmeans(nonlinear_model3A1 , specs = pairwise ~ Species)
EMMm
EMMplotC= plot(EMMm)
save(nonlinear_model3A1, file = "nonlinear_model3A1.RData")

nonlinear_model3A2=brm(bf((Intereff+.0001)~s(DiptDensity_scaled)+Clade+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT4H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 5000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model3A2, file = "nonlinear_model3A2.RData")
summary(nonlinear_model3A2)
plot(conditional_effects(nonlinear_model3A2,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model3A2)
nonlinear_model3A2 = add_criterion(nonlinear_model3A2, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model3A2)
deriv_plot(model=nonlinear_model3A2,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model3A2_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_model3A2,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model3A2_DiptDensity_scaled_plot2")
EMMn=emmeans(nonlinear_model3A2 , specs = pairwise ~ Clade)
EMMn
plot(EMMn)
save(nonlinear_model3A2, file = "nonlinear_model3A2.RData")


nonlinear_model3A3=brm(bf((Intereff+.0001)~s(DiptDensity_scaled)+s(EQ, k=4)+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT4H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 3000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model3A3, file = "nonlinear_model3A3.RData")
summary(nonlinear_model3A3)
plot(conditional_effects(nonlinear_model3A3,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model3A3)
nonlinear_model3A3 = add_criterion(nonlinear_model3A3, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model3A3)
deriv_plot(model=nonlinear_model3A3,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model3A3_DiptDensity_scaled_plot")
deriv_plot(model=nonlinear_model3A3,term="s(EQ, k=4)",main="EQ",eps=0.001,output="nonlinear_model3A3_EQ_plot")
save(nonlinear_model3A3, file = "nonlinear_model3A3.RData")
deriv_plot2(model=nonlinear_model3A3,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model3A3_DiptDensity_scaled_plot2")
deriv_plot2(model=nonlinear_model3A3,term="s(EQ, k=4)",main="EQ",eps=0.001,output="nonlinear_model3A3_EQ_plot2")
save(nonlinear_model3A3, file = "nonlinear_model3A3.RData")



## Class 4B - 4hr TBeff models

nonlinear_model4B0=brm(bf((TBeff_interp+.0001)~s(DiptDensity_scaled)+Species+Sex:Species+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT4H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 8000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model4B0, file = "nonlinear_model4B0.RData")
summary(nonlinear_model4B0)
plot(conditional_effects(nonlinear_model4B0,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model4B0)
nonlinear_model4B0 = add_criterion(nonlinear_model4B0, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model4B0)
deriv_plot(model=nonlinear_model4B0,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model4B0_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_model4B0,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model4B0_DiptDensity_scaled_plot2")
EMMo=emmeans(nonlinear_model4B0 , specs = pairwise ~ Species)
EMMo
plot(EMMo)
EMMo2=emmeans(nonlinear_model4B0 , specs = pairwise ~ Sex:Species)
EMMo2
plot(EMMo2)
save(nonlinear_model4B0, file = "nonlinear_model4B0.RData")

nonlinear_model4B1=brm(bf((TBeff_interp+.0001)~s(DiptDensity_scaled)+Species+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT4H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 3000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model4B1, file = "nonlinear_model4B1.RData")
summary(nonlinear_model4B1)
plot(conditional_effects(nonlinear_model4B1,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model4B1)
nonlinear_model4B1 = add_criterion(nonlinear_model4B1, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model4B1)
deriv_plot(model=nonlinear_model4B1,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model4B1_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_model4B1,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model4B1_DiptDensity_scaled_plot2")
EMMp=emmeans(nonlinear_model4B1 , specs = pairwise ~ Species)
EMMp
EMMplotD= plot(EMMp)
save(nonlinear_model4B1, file = "nonlinear_model4B1.RData")

nonlinear_model4B2=brm(bf((TBeff_interp+.0001)~s(DiptDensity_scaled)+Clade+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT4H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 5000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model4B2, file = "nonlinear_model4B2.RData")
summary(nonlinear_model4B2)
plot(conditional_effects(nonlinear_model4B2,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model4B2)
nonlinear_model4B2 = add_criterion(nonlinear_model4B2, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model4B2)
deriv_plot(model=nonlinear_model4B2,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model4B2_DiptDensity_scaled_plot")
deriv_plot2(model=nonlinear_model4B2,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model4B2_DiptDensity_scaled_plot2")
EMMq=emmeans(nonlinear_model4B2 , specs = pairwise ~ Clade)
EMMq
plot(EMMq)
save(nonlinear_model4B2, file = "nonlinear_model4B2.RData")

nonlinear_model4B3=brm(bf((TBeff_interp+.0001)~s(DiptDensity_scaled)+s(EQ, k=4)+(1 + DiptDensity_scaled | ID), decomp = "QR"), data=DailyFFT4H,
                       family = Gamma(link="log"),
                       init=0,
                       iter = 3000,
                       chains=2,
                       prior=c(prior(normal(0,1.45),class="Intercept"),
                               prior(normal(0,1.45),class="b"),
                               prior(gamma(.01,.01),class="shape"),
                               prior(normal(0,1.45),class="sd")),
                       save_pars = save_pars(all = TRUE),
                       backend = "cmdstanr",
                       threads = threading(2),
                       control=list(adapt_delta=0.9999, max_treedepth = 10))
save(nonlinear_model4B3, file = "nonlinear_model4B3.RData")
summary(nonlinear_model4B3)
plot(conditional_effects(nonlinear_model4B3,spaghetti=TRUE)) ##plots the regression lines/conditional effects
pp_check(nonlinear_model4B3)
nonlinear_model4B3 = add_criterion(nonlinear_model4B3, c("waic","bayes_R2","loo", "loo_R2"), 
                                   control=list(adapt_delta=0.9999, max_treedepth = 10), 
                                   backend = "cmdstanr",
                                   reloo = TRUE)
loo_R2(nonlinear_model4B3)
deriv_plot(model=nonlinear_model4B3,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model4B3_DiptDensity_scaled_plot")
deriv_plot(model=nonlinear_model4B3,term="s(EQ, k=4)",main="EQ",eps=0.001,output="nonlinear_model4B3_EQ_plot")
deriv_plot2(model=nonlinear_model4B3,term="s(DiptDensity_scaled)",main="DiptDensity_scaled",eps=0.001,output="nonlinear_model4B3_DiptDensity_scaled_plot2")
deriv_plot2(model=nonlinear_model4B3,term="s(EQ, k=4)",main="EQ",eps=0.001,output="nonlinear_model4B3_EQ_plot2")
save(nonlinear_model4B3, file = "nonlinear_model4B3.RData")


##Commence with model comparison using LOOCV

loo_compare(nonlinear_modelA0,nonlinear_modelA1,nonlinear_modelA2,nonlinear_modelA3,criterion = "loo") #nonlinear_modelA1 willing model
loo_compare(nonlinear_model2A0,nonlinear_model2A1,nonlinear_model2A2,nonlinear_model2A3,criterion = "loo") #nonlinear_model2A0 winning model
loo_compare(nonlinear_model3A0,nonlinear_model3A1,nonlinear_model3A2,nonlinear_model3A3,criterion = "loo") #nonlinear_model3A1 winning model
loo_compare(nonlinear_modelB0,nonlinear_modelB1,nonlinear_modelB2,nonlinear_modelB3,criterion = "loo") #nonlinear_modelB1 winning model
loo_compare(nonlinear_model2B0,nonlinear_model2B1,nonlinear_model2B2,nonlinear_model2B3,criterion = "loo") #nonlinear_model2B1 winning model
loo_compare(nonlinear_model4B0,nonlinear_model4B1,nonlinear_model4B2,nonlinear_model4B3,criterion = "loo") #nonlinear_model4B1 winning model

loo(nonlinear_modelA1)

##Edit derivative plots
nonlinear_modelA1_DiptDensity_scaled_plot2 = nonlinear_modelA1_DiptDensity_scaled_plot2+
  ggtitle("A")+ ylab("First Derivative") + xlab("Dipteryx Density") + theme(legend.position = "none")
nonlinear_modelB1_DiptDensity_scaled_plot2 = nonlinear_modelB1_DiptDensity_scaled_plot2+
  ggtitle("B")+ ylab("First Derivative") + xlab("Dipteryx Density")+ theme(legend.position = "none")
nonlinear_model3A1_DiptDensity_scaled_plot2 = nonlinear_model3A1_DiptDensity_scaled_plot2+
  ggtitle("C")+ ylab("First Derivative") + xlab("Dipteryx Density")+ theme(legend.position = "none")
nonlinear_model4B1_DiptDensity_scaled_plot2 = nonlinear_model4B1_DiptDensity_scaled_plot2+
  ggtitle("D")+ ylab("First Derivative") + xlab("Dipteryx Density")+ theme(legend.position = "none")
nonlinear_model2A1_DiptDensity_scaled_plot2 = nonlinear_model2A1_DiptDensity_scaled_plot2+
  ggtitle("E")+ ylab("First Derivative") + xlab("Dipteryx Density")+ theme(legend.position = "none")
nonlinear_model2B0_DiptDensity_scaled_plot2 = nonlinear_model2B0_DiptDensity_scaled_plot2+
  ggtitle("F")+ ylab("First Derivative") + xlab("Dipteryx Density")+ theme(legend.position = "none")

#Plot in multipanel 
ggpubr::ggarrange(nonlinear_modelA1_DiptDensity_scaled_plot2, nonlinear_modelB1_DiptDensity_scaled_plot2, nonlinear_model3A1_DiptDensity_scaled_plot2, nonlinear_model4B1_DiptDensity_scaled_plot2, nonlinear_model2A1_DiptDensity_scaled_plot2, nonlinear_model2B0_DiptDensity_scaled_plot2,ncol=2, nrow=3)

#edit emmeans plots
EMMplotA = EMMplotA + ggtitle("A")+ ylab("") + xlab("Efficiency (EMM)") + theme(legend.position = "none")
EMMplotB = EMMplotB + ggtitle("B")+ ylab("") + xlab("Efficiency (EMM)") + theme(legend.position = "none")
EMMplotC = EMMplotC + ggtitle("C")+ ylab("") + xlab("Efficiency (EMM)") + theme(legend.position = "none")
EMMplotD = EMMplotD + ggtitle("D")+ ylab("") + xlab("Efficiency (EMM)") + theme(legend.position = "none")
EMMplotE = EMMplotE + ggtitle("E")+ ylab("") + xlab("Efficiency (EMM)") + theme(legend.position = "none")
EMMplotF = EMMplotF + ggtitle("F")+ ylab("") + xlab("Efficiency (EMM)") + theme(legend.position = "none")

#Plot as multipanel 
ggpubr::ggarrange(EMMplotA, EMMplotB, EMMplotC, EMMplotD, EMMplotE, EMMplotF,ncol=2, nrow=3)


