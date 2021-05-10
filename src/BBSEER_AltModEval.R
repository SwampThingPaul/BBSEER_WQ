## 
## BBSEER
## WQ alternative evaluation
## 
##
## Code was compiled by Paul Julian
## contact info: paul.julian@floridadep.gov

## BAD ## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
# devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);
library(plyr)
library(reshape)

#GIS Libraries
# library(rgdal)
# library(rgeos)
# library(tmap)
# library(raster)

library(Hmisc)
library(dunn.test)

library(flextable)
library(magrittr)

wd="C:/Julian_LaCie/_GitHub/BBSEER_WQ"

paths=paste0(wd,c("/Plots/","/Exports/","/Data/","/src/","/GIS"))
# Folder.Maker(paths);#One and done. Creates folders in working DIR
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]

# -------------------------------------------------------------------------
dates=date.fun(c("1999-05-01","2020-04-30"))

sites=data.frame(SITE=c("S20F","S20G","S21A","S21","S123","S22","S26","S25A","S25B","S27","S28","S29"),
                 DBKEY=c(91441,91442,91445,91447,91367,91449,91469,91457,91456,91470,91471,91472),
                 WQSite=c("MW01","MI01","PR01","BL02","CD01A","SP01","MR01","MR01","MR01","LR01","BS01","SK01"),
                 Canal=c("C103","G95","C102","C1","C100","C2","C6","C5","C4","C7","C8","C9"))

storet_win=data.frame(SITE=c("S20F","S20G","S21A","S21","S123","S22","S26","S25A","S25B","S27","S28","S29"),
                      WQSite=c("MW01","MI01","PR01","BL02","CD01A","SP01","MR01","MR01","MR01","LR01","BS01","SK01"))

wq.nnc=data.frame(WQSite=c("MW01","MI01","PR01","BL02","CD01A","SP01","MR01","LR01","BS01","SK01"),
                  NNC_seg=paste0("ENRH",c(rep(6,4),3,3,9,9,5,5)))

# Water Quality -----------------------------------------------------------
## https://floridadep.gov/sites/default/files/62-160_help-document_0.pdf
dat.qual=data.frame(QUALIFIER=c(NA,"!","A","D","E","F","I","R","T","U","*","?","B","H","J","K","L","M","N","O","Q","V","Y","Z"),
                    FATALYN=c("N",rep("N",9),rep("Y",14)))

## STORET data (downloaded 2020-12-11;)
## https://floridadep.gov/dear/watershed-services-program/content/winstoret
storet.dat1=read.table(paste0(data.path,"STORET/Water_Quality_Results_dwn20210208.txt"),sep="|",header=T)
storet.dat2=read.table(paste0(data.path,"STORET/Water_Quality_Results_dwn20201211_2.txt"),sep="|",header=T)
storet.dat=rbind(storet.dat1,storet.dat2)
storet.dat$Date.EST=date.fun(storet.dat$Act.Date,form="%m/%d/%Y")
storet.dat=subset(storet.dat,Org.Name%in%c("Dade Environmental Resource Management (Florida)", "South Florida Water Management District"))

range(storet.dat$Date.EST)
unique(storet.dat$Medium)
unique(storet.dat$Matrix)
unique(storet.dat$Act.Type)
unique(storet.dat$Act.Category)
# QA/QC qualifiers 
quals=as.character(unique(storet.dat$VQ))
spl=strsplit(quals,split="")
quals=data.frame(VQ=quals,q1=sapply(spl,"[",1),q2=sapply(spl,"[",2),q3=sapply(spl,"[",3))
quals$Fatal=with(quals,ifelse(q1%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q2%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q3%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER,"Y","N"))

sort(unique(storet.dat$Result.Value))
subset(storet.dat,Result.Value=="*Non-detect")
subset(storet.dat,Result.Value==0)

storet.dat$Result.Value=with(storet.dat,ifelse(Result.Value=="*Non-detect",MDL*-1,as.numeric(as.character(Result.Value))))
subset(storet.dat,is.na(Result.Value))
storet.dat$HalfMDL=with(storet.dat,ifelse(Act.Type=="Field Msr/Obs",Result.Value,
                                          ifelse(abs(Result.Value)<=MDL,MDL/2,abs(Result.Value))))
sort(unique(storet.dat$HalfMDL))
storet.dat=merge(storet.dat,quals[,c("VQ","Fatal")],"VQ",all.x=T)
storet.dat=subset(storet.dat,Fatal=="N")
unique(as.character(storet.dat$Result.Units))
dput(unique(as.character(storet.dat$Characteristic)))

ddply(storet.dat,c("Characteristic","Result.Units"),summarise,n.val=N.obs(Station.ID))
storet.parameters=data.frame(Characteristic=c("Turbidity", "Total Suspended Solids (TSS)", "Nitrogen, ammonia (NH3) as NH3", 
                                              "Phosphorus as PO4", "Nitrogen, Nitrite (NO2) + Nitrate (NO3) as N", 
                                              "Phosphorus, phosphate (PO4) as PO4","Chlorophyll a, free of pheophytin", 
                                              "Apparent Color","Nitrogen, Kjeldahl","Phosphorus, orthophosphate as PO4","pH", "Dissolved oxygen saturation", 
                                              "Dissolved oxygen (DO)", "Salinity", "Specific conductance","Temperature, water", "Nitrogen, ammonia as N","Phosphorus, orthophosphate as P","Phosphorus","Phosphorus as P","Chlorophyll a, corrected for pheophytin"),
                             param=c("Turb","TSS","NH4","TP","NOx","TP","Chla_c","Color","TKN","SRP","pH","DOSat","DO","Sal","SPC","Temp","NH4","SRP","TP","TP","Chla_c"))
unique(as.character(storet.dat$Characteristic))[unique(as.character(storet.dat$Characteristic))%in%storet.parameters$Characteristic==FALSE]

nrow(storet.dat)
storet.dat=merge(storet.dat,storet.parameters,"Characteristic")

storet.xtab=cast(storet.dat,Station.ID+Date.EST~param,value="HalfMDL",mean)
storet.xtab$WY=WY(storet.xtab$Date.EST)
storet.xtab$TN=NA
storet.xtab$TN=with(storet.xtab,TN_Combine(NOx,TKN,TN))
storet.xtab$DIN=with(storet.xtab,NH4+NOx)

storet.xtab$TPReversal=with(storet.xtab,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>(TP*1.3),1,0)));# Reversals identified as 1 reversals consistent with TP rule evaluation
storet.xtab$TNReversal=with(storet.xtab,ifelse(is.na(DIN)==T|is.na(TN)==T,0,ifelse(DIN>(TN*1.3),1,0)));

sum(storet.xtab$TNReversal,na.rm=T)
sum(storet.xtab$TPReversal,na.rm=T)
subset(storet.xtab,TPReversal==1)

par(family="serif",oma=c(1,1,1,1),mar=c(4,4,1,1))
layout(matrix(1:2,1,2,byrow=F))
plot(TN~DIN,storet.xtab,ylab="TN (mg L\u207B\u00b9)",xlab="DIN (mg L\u207B\u00b9)",pch=21,bg=ifelse(TNReversal==1,"dodgerblue1",NA),col=adjustcolor("grey",0.8));abline(0,1,col="dodgerblue1")
plot(TP~SRP,storet.xtab,ylab="TP (mg L\u207B\u00b9)",xlab="SRP (mg L\u207B\u00b9)",pch=21,bg=ifelse(TPReversal==1,"red",NA),col=adjustcolor("grey",0.8));abline(0,1,col="red")

storet.xtab$TP=with(storet.xtab,ifelse(TPReversal==1,NA,TP))
storet.xtab$SRP=with(storet.xtab,ifelse(TPReversal==1,NA,SRP))
storet.xtab$TN=with(storet.xtab,ifelse(TNReversal==1,NA,TN))
storet.xtab$DIN=with(storet.xtab,ifelse(TNReversal==1,NA,DIN))

plot(TP~Date.EST,subset(storet.xtab,Station.ID=="SK01"))
plot(TN~Date.EST,subset(storet.xtab,Station.ID=="SK01"))
plot(NH4~Date.EST,subset(storet.xtab,Station.ID=="SK01"))
plot(TN~Date.EST,subset(storet.xtab,Station.ID=="MW01"))
dev.off()

## WIN data (downloaded 2020-12-01)
win.dat=read.table(paste0(data.path,"WIN/WIN_WAVES_GUEST_20201211072815_36488.txt"),skip=7,sep="|",header=T)
win.dat$Activity.Start.Date.Time=date.fun(as.character(win.dat$Activity.Start.Date.Time),form="%m/%d/%Y %H:%M:%S")
win.dat$Date.EST=date.fun(win.dat$Activity.Start.Date.Time)
win.dat$Station.ID=win.dat$Monitoring.Location.ID
unique(win.dat$Organization.ID)

# QA/QC qualifiers 
quals=as.character(unique(win.dat$Value.Qualifier))
spl=strsplit(quals,split="")
quals=data.frame(Value.Qualifier=quals,q1=sapply(spl,"[",1),q2=sapply(spl,"[",2),q3=sapply(spl,"[",3),q4=sapply(spl,"[",4))
quals$Fatal=with(quals,ifelse(q1%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q2%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q3%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q4%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER,"Y","N"))

unique(win.dat$Matrix)
win.dat.clean=subset(win.dat,Matrix=="AQUEOUS-Surface Water")
win.dat.clean$HalfMDL=with(win.dat.clean,ifelse(Sample.Collection.Type=="Field Testing-Discrete",DEP.Result.Value.Number,
                                                ifelse(DEP.Result.Value.Number<=DEP.MDL,DEP.MDL/2,DEP.Result.Value.Number)))
win.dat.clean=merge(win.dat.clean,quals[,c("Value.Qualifier","Fatal")],"Value.Qualifier",all.x=T)
win.dat.clean=subset(win.dat.clean,Fatal=="N")

dput(unique(win.dat.clean$DEP.Analyte.Name))
win.parameters=data.frame(DEP.Analyte.Name=c("Ammonia (N)", 
                                             "Chlorophyll a- uncorrected", 
                                             "Chlorophyll a, free of pheophytin","Chlorophyll a- corrected", "Dissolved Oxygen","Dissolved Oxygen Saturation", "Nitrate-Nitrite (N)", 
                                             "Nitrogen- Total Kjeldahl", "Orthophosphate (P)", "Phosphorus- Total", 
                                             "Specific Conductance", "Temperature, Water","pH","Salinity","Turbidity","Residues- Nonfilterable (TSS)","Color- Apparent"),
                          param=c("NH4","Chla","Chla_c","Chla_c","DO","DOSat","NOx","TKN","SRP","TP","SPC","Temp","pH","Sal","Turb","TSS","Color"))
unique(win.dat.clean$DEP.Analyte.Name)%in%win.parameters$DEP.Analyte.Name

win.dat.clean=merge(win.dat.clean,win.parameters,"DEP.Analyte.Name")

win.xtab=cast(win.dat.clean,Station.ID+Date.EST~param,value="HalfMDL",mean)
win.xtab$WY=WY(win.xtab$Date.EST)
win.xtab$TN=NA
win.xtab$TN=with(win.xtab,TN_Combine(NOx,TKN,TN))
win.xtab$DIN=with(win.xtab,NH4+NOx)

win.xtab$TPReversal=with(win.xtab,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>(TP*1.3),1,0)));# Reversals identified as 1 reversals consistent with TP rule evaluation
win.xtab$TNReversal=with(win.xtab,ifelse(is.na(DIN)==T|is.na(TN)==T,0,ifelse(DIN>(TN*1.3),1,0)));

sum(win.xtab$TNReversal,na.rm=T)
sum(win.xtab$TPReversal,na.rm=T)
subset(win.xtab,TPReversal==1)

par(family="serif",oma=c(1,1,1,1),mar=c(4,4,1,1))
layout(matrix(1:2,1,2,byrow=F))
plot(TN~DIN,win.xtab,ylab="TN (mg L\u207B\u00b9)",xlab="DIN (mg L\u207B\u00b9)",pch=21,bg=ifelse(TNReversal==1,"dodgerblue1",NA),col=adjustcolor("grey",0.8));abline(0,1,col="dodgerblue1")
plot(TP~SRP,win.xtab,ylab="TP (mg L\u207B\u00b9)",xlab="SRP (mg L\u207B\u00b9)",pch=21,bg=ifelse(TPReversal==1,"red",NA),col=adjustcolor("grey",0.8));abline(0,1,col="red")

win.xtab$TP=with(win.xtab,ifelse(TPReversal==1,NA,TP))
win.xtab$SRP=with(win.xtab,ifelse(TPReversal==1,NA,SRP))
win.xtab$TN=with(win.xtab,ifelse(TNReversal==1,NA,TN))
win.xtab$DIN=with(win.xtab,ifelse(TNReversal==1,NA,DIN))

plot(TP~Date.EST,subset(win.xtab,Station.ID=="SK01"))
plot(TN~Date.EST,subset(win.xtab,Station.ID=="SK01"))

## 
names(win.xtab)[names(win.xtab)%in%names(storet.xtab)==F]
storet.xtab$Chla=NA
storet.xtab$Chla_c=NA
storet.xtab$Color=NA

storet.xtab=storet.xtab[,names(win.xtab)]
dat.xtab=rbind(storet.xtab,win.xtab)
dat.xtab$season=FL.Hydroseason(dat.xtab$Date.EST)

plot(TP~Date.EST,subset(dat.xtab,Station.ID=="MW01"))
plot(TN~Date.EST,subset(dat.xtab,Station.ID=="MW01"))

dat.xtab=merge(dat.xtab,storet_win,by.x="Station.ID",by.y="WQSite",all.x=T)

## DBHydro
dbhydro.sites=subset(sites,is.na(WQSite)==F)[,c("SITE","WQSite")]
params=data.frame(Test.Number=c(18,21,80,20,25,23,61,179,7,16,8,9,98,12,13),
                  param=c("NOx","TKN","TN","NH4","TP","SRP","Chla","Chla","Temp","TSS","DO","SPC","Sal","Turb","Color"))

wmd.dat=DBHYDRO_WQ(dates[1],dates[2],dbhydro.sites$WQSite,params$Test.Number)
wmd.dat=merge(wmd.dat,params,"Test.Number")
unique(wmd.dat$Collection.Method)
wmd.dat=subset(wmd.dat,Collection.Method=="G")

wmd.dat.xtab=cast(wmd.dat,Station.ID+Date.EST~param,value="HalfMDL",mean)
wmd.dat.xtab$DIN=with(wmd.dat.xtab,NH4+NOx)
wmd.dat.xtab$TN=with(wmd.dat.xtab, TN_Combine(NOx,TKN,TN))
wmd.dat.xtab$DOSat=with(wmd.dat.xtab,DO_PerSat(Temp,DO,Sal))
wmd.dat.xtab$WY=WY(wmd.dat.xtab$Date.EST)

wmd.dat.xtab$TPReversal=with(wmd.dat.xtab,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>(TP*1.3),1,0)));# Reversals identified as 1 reversals consistent with TP rule evaluation
wmd.dat.xtab$TNReversal=with(wmd.dat.xtab,ifelse(is.na(DIN)==T|is.na(TN)==T,0,ifelse(DIN>(TN*1.3),1,0)));

# sum(wmd.dat.xtab$TNReversal,na.rm=T)
# sum(wmd.dat.xtab$TPReversal,na.rm=T)
# subset(wmd.dat.xtab,TPReversal==1)

#par(family="serif",oma=c(1,1,1,1),mar=c(4,4,1,1))
#layout(matrix(1:2,1,2,byrow=F))
plot(TN~DIN,wmd.dat.xtab,ylab="TN (mg L\u207B\u00b9)",xlab="DIN (mg L\u207B\u00b9)",pch=21,bg=ifelse(TNReversal==1,"dodgerblue1",NA),col=adjustcolor("grey",0.8));abline(0,1,col="dodgerblue1")
plot(TP~SRP,wmd.dat.xtab,ylab="TP (mg L\u207B\u00b9)",xlab="SRP (mg L\u207B\u00b9)",pch=21,bg=ifelse(TPReversal==1,"red",NA),col=adjustcolor("grey",0.8));abline(0,1,col="red")

names(dat.xtab)[names(dat.xtab)%in%names(wmd.dat.xtab)==F]

wmd.dat.xtab$Chla_c=as.numeric(NA)
wmd.dat.xtab$pH=as.numeric(NA)
wmd.dat.xtab$season=FL.Hydroseason(wmd.dat.xtab$Date.EST)
wmd.dat.xtab=merge(wmd.dat.xtab,dbhydro.sites,by.x="Station.ID",by.y="WQSite",all.x=T)

wmd.dat.xtab=wmd.dat.xtab[,names(dat.xtab)]
str(wmd.dat.xtab)
str(dat.xtab)
dat.xtab=rbind(dat.xtab,wmd.dat.xtab)

head(dat.xtab)
IDvars=c("Station.ID","Date.EST","WY","season","SITE")
paramsvars=c("Chla", "Chla_c", "Color", "DO","DOSat", "NH4", "NOx", "pH", "Sal", "SPC", "SRP", "Temp", "TKN","TP", "TSS", "Turb", "WY", "TN", "DIN")
dat.xtab.melt=melt(dat.xtab[,c(IDvars,paramsvars)],id.vars = IDvars)
dat.xtab.melt=subset(dat.xtab.melt,is.na(value)==F)
dat.xtab=cast(dat.xtab.melt,Station.ID+Date.EST+WY+season+SITE~variable,value="value",mean)
dat.xtab$TP=dat.xtab$TP*1000
dat.xtab$month=as.numeric(format(dat.xtab$Date.EST,"%m"))


# POR WQ ------------------------------------------------------------------
POR.month.wq=ddply(dat.xtab,c("Station.ID","month"),summarise,mean.TP=mean(TP,na.rm=T),sd.TP=sd(TP,na.rm=T),Q10=quantile(TP,probs=0.10,na.rm=T),Q90=quantile(TP,probs=0.9,na.rm=T),min.val=min(TP,na.rm=T),max.val=max(TP,na.rm=T))

POR.WY.AM.wq=ddply(dat.xtab,c("Station.ID","WY"),summarise,mean.TP=mean(TP,na.rm=T))
POR.year.wq=ddply(dat.xtab,c("Station.ID"),summarise,sd.Yr.TP=sd(TP,na.rm=T))
            

test.S29.jan=subset(POR.month.wq,SITE=="S29"&month==1)

set.seed(123)
nsim=500
#tmp=rlnorm(nsim,log(test.S29.jan$mean.TP),log(test.S29.jan$sd.TP))
tmp=with(test.S29.jan,runif(nsim,mean.TP-sd.TP,mean.TP+sd.TP))
tmp

plot(tmp)
abline(h=test.S29.jan$mean.TP,col="red")
mean(tmp)

tmp2=with(test.S29.jan,rnorm(nsim,mean.TP,sd.TP))

plot(tmp2)
abline(h=test.S29.jan$mean.TP,col="red")
mean(tmp2)


# Discharge ---------------------------------------------------------------
q.dat=data.frame()
for(i in 1:nrow(sites)){
  tmp=DBHYDRO_daily(dates[1],dates[2],sites$DBKEY[i])
  tmp$DBKEY=as.character(sites$DBKEY[i])
  q.dat=rbind(q.dat,tmp)
  print(i)
}

q.dat$WY=WY(q.dat$Date)
q.dat$Date.EST=date.fun(q.dat$Date)
q.dat=merge(q.dat,sites,"DBKEY")

range(q.dat$Data.Value,na.rm=T)
q.dat$Data.Value=with(q.dat,ifelse(Data.Value<0,0,Data.Value));# positive flow only

q.dat.tflow=ddply(q.dat,c("WQSite","Date.EST","WY"),summarise,TFlow.cfs=sum(Data.Value,na.rm=T))

flow.wq=merge(q.dat.tflow,dat.xtab,by.x=c("WQSite","Date.EST","WY"),by.y=c("Station.ID","Date.EST","WY"),all.x=T)
flow.wq=merge(flow.wq,wq.nnc,"WQSite")

plot(TN~Date.EST,subset(flow.wq,WQSite=="MW01"))

flow.wq$TP.int=with(flow.wq,ave(TP,WQSite,FUN=function(x) dat.interp(x)))
flow.wq$TP.load=with(flow.wq,Load.Calc.kg(TFlow.cfs,TP.int/1000))
flow.wq$TN.int=with(flow.wq,ave(TN,WQSite,FUN=function(x) dat.interp(x)))
flow.wq$TN.load=with(flow.wq,Load.Calc.kg(TFlow.cfs,TN.int))

flow.wq$month=as.numeric(format(flow.wq$Date.EST,"%m"))
flow.wq$CY=as.numeric(format(flow.wq$Date.EST,"%Y"))

TLoad.obs=ddply(flow.wq,c("WY","WQSite"),summarise,TFlow=sum(cfs.to.acftd(TFlow.cfs),na.rm=T),TPLoad.kg=sum(TP.load,na.rm=T),TNLoad.kg=sum(TN.load,na.rm=T))
TLoad.obs$TP.FWM.ugL=with(TLoad.obs,(TPLoad.kg/(TFlow*1.233e6))*1e9)
TLoad.obs$TN.FWM.mgL=with(TLoad.obs,(TNLoad.kg/(TFlow*1.233e6))*1e6)


cv.vals=ddply(flow.wq,c("WY","WQSite"),summarise,TP.CV=cv.per(TP),TN.CV=cv.per(TN),Q.CV=cv.per(TFlow.cfs))
# cv.vals=ddply(flow.wq,c("WY","WQSite"),summarise,TP.CV=EnvStats::cv(TP,na.rm=T),TN.CV=EnvStats::cv(TN,na.rm=T),Q.CV=EnvStats::cv(TFlow.cfs,na.rm=T))
site.val=rev(c("MW01","MI01","PR01","BL02","CD01A","SP01","MR01","LR01","BS01","SK01"))
discharge.site=rev(c("S20F","S20G","S21A","S21","S123","S22","S26/S25A/S25B","S27","S28","S29"))

par(family="serif",mar=c(2,1.75,1,2),oma=c(2,2,0.75,2))
layout(matrix(1:10,5,2))

for(k in 1:length(site.val)){
  plot(TP.CV~WY,cv.vals,xlim=c(2005,2020),ylim=c(0.1,15),type="n",log="y")
  with(subset(cv.vals,WQSite==site.val[k]),lines(WY,TP.CV,col="red"))
  with(subset(cv.vals,WQSite==site.val[k]),lines(WY,TN.CV,col="blue"))
  with(subset(cv.vals,WQSite==site.val[k]),lines(WY,Q.CV,col="black"))
  
  mtext(side=3,site.val[k])
}


par(family="serif",mar=c(2,1.75,1,2),oma=c(2,2,0.75,2))
layout(matrix(1:10,5,2))
for(k in 1:length(site.val)){
  plot(TP.CV~Q.CV,cv.vals,xlim=c(0.5,15),type="n",log="x")
  with(subset(cv.vals,WQSite==site.val[k]),points(Q.CV,TP.CV,pch=21,bg="red"))
  with(subset(cv.vals,WQSite==site.val[k]),points(Q.CV,TN.CV,pch=21,bg="blue"))
  mtext(side=3,site.val[k])
}



# -------------------------------------------------------------------------
flow.wq=subset(flow.wq,WY%in%seq(2005,2020,1))
TLoad.obs=subset(TLoad.obs,WY%in%seq(2005,2020,1))

unique(flow.wq$WQSite)

site.val=rev(c("MW01","MI01","PR01","BL02","CD01A","SP01","MR01","LR01","BS01","SK01"))
discharge.site=rev(c("S20F","S20G","S21A","S21","S123","S22","S26/S25A/S25B","S27","S28","S29"))

monthmean.load.all=data.frame()
for(j in 1:length(site.val)){
tmp.Q=subset(flow.wq,WQSite==site.val[j])
tmp.wq=subset(POR.month.wq,Station.ID==site.val[j])

tmp.Q.wq=tmp.Q.wq[order(tmp.Q.wq$WQSite,tmp.Q.wq$Date.EST),]
tmp.Q.wq=merge(tmp.Q,tmp.wq[,c("month","mean.TP")],all.x=T,"month")
tmp.Q.wq$TPLoad=with(tmp.Q.wq,Load.Calc.kg(TFlow.cfs,mean.TP/1000))

tmp.WY.load=ddply(tmp.Q.wq,c("WY","WQSite"),summarise,TFlow=sum(cfs.to.acftd(TFlow.cfs),na.rm=T),TPLoad.kg=sum(TPLoad,na.rm=T))
monthmean.load.all=rbind(monthmean.load.all,tmp.WY.load)

}
monthmean.load.all$TP.FWM.ugL=with(monthmean.load.all,(TPLoad.kg/(TFlow*1.233e6))*1e9)

monthmean.load.all=monthmean.load.all[order(monthmean.load.all$WQSite,monthmean.load.all$WY),]

par(family="serif",mar=c(2,1.75,1,2),oma=c(2,2,0.75,2))
layout(matrix(1:10,5,2))

for(k in 1:length(site.val)){
  plot(TP.FWM.ugL~WY,subset(TLoad.obs,WQSite==site.val[k]),xlim=c(2005,2020),ylim=c(0,30),type="n")
  with(subset(TLoad.obs,WQSite==site.val[k]),lines(WY,TP.FWM.ugL,col="red"))
  with(subset(monthmean.load.all,WQSite==site.val[k]),lines(WY,TP.FWM.ugL,col="blue"))
  mtext(side=3,site.val[k])
}


# Monte Carlo -------------------------------------------------------------
set.seed(123)

MC.load.all=data.frame()
conc.sim=data.frame()
nsim=50
for(j in 1:length(site.val)){
pb <- txtProgressBar(min = 0, max = nsim, style = 3)
for(i in 1:nsim){
tmp.Q=subset(flow.wq,WQSite==site.val[j])
tmp.wq=subset(POR.month.wq,Station.ID==site.val[j])
# tmp.wq$TP.dist=apply(tmp.wq,1,FUN=function(x) rnorm(1,(tmp.wq$mean.TP),(tmp.wq$sd.TP))) # +rnorm(1,0,sd(tmp.wq$mean.TP))
# tmp.wq$TP.dist=apply(tmp.wq,1,FUN=function(x) runif(1,tmp.wq$mean.TP-tmp.wq$sd.TP,tmp.wq$mean.TP+tmp.wq$sd.TP))+rnorm(12,0,sd(tmp.wq$mean.TP))#+rnorm(1,0,sd(subset(POR.WY.AM.wq,Station.ID==site.val[j])$mean.TP))
# tmp.wq$TP.dist=apply(tmp.wq,1,FUN=function(x) runif(1,tmp.wq$Q10,tmp.wq$Q90))+rnorm(1,0,sd(tmp.wq$mean.TP))
conc.sim=rbind(data.frame(tmp.wq[,c("month","Station.ID","mean.TP","TP.dist")],iter=i),conc.sim)

tmp.Q.wq=merge(tmp.Q,tmp.wq[,c("month","TP.dist")],all.x=T,"month")
tmp.Q.wq$TPLoad.unif=with(tmp.Q.wq,Load.Calc.kg(TFlow.cfs,TP.dist/1000))

tmp.WY.load=ddply(tmp.Q.wq,c("WY","WQSite"),summarise,TFlow=sum(cfs.to.acftd(TFlow.cfs),na.rm=T),TPLoad.kg=sum(TPLoad.unif,na.rm=T))
tmp.WY.load$iter=i
MC.load.all=rbind(MC.load.all,tmp.WY.load)
setTxtProgressBar(pb, i)
}
print(j)
}

ddply(conc.sim,c("Station.ID","iter"),summarise,mean.val=mean(TP.dist))
MC.load.all$TP.FWM.ugL=with(MC.load.all,(TPLoad.kg/(TFlow*1.233e6))*1e9)

tmp.x1=subset(monthmean.load.all,WQSite==site.val[2])
tmp.y1=ddply(subset(MC.load.all,WQSite==site.val[2]),"WY",summarise,mean.val=mean(TP.FWM.ugL))

plot(tmp.y1$mean.val~tmp.x1$TP.FWM.ugL,ylab="Monte Carlo Sim FWM",xlab="Monthly POR Mean FWM")
plot(tmp.y1$mean.val~subset(TLoad.obs,WQSite==site.val[1])$TP.FWM.ugL,ylab="Monte Carlo Sim FWM",xlab="Obs FWM")


xlim.val=c(2005,2020);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0,8000);by.y=2000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
max.val=ddply(TLoad.obs,"WQSite",summarise,max.val=max(TPLoad.kg))
max.val[order(match(max.val$WQSite,site.val)),]
ylim.max.val=c(10000,2500,6000,18000,3500,3000,6000,3000,500,5000)

# png(filename=paste0(plot.path,"SimLoad_TP.png"),width=6,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,1.75,0.75,2),oma=c(3,3,0.75,2))
layout(matrix(c(1:10),5,2))

for(k in 1:length(site.val)){
ylim.val=c(0,ylim.max.val[k]);by.y=ylim.max.val[k]/2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TPLoad.kg~WY,MC.load.all,ylim=ylim.val,xlim=xlim.val,axes=F,ann=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")

with(subset(MC.load.all,WQSite==site.val[k]),points(WY,TPLoad.kg,pch=21,bg=NA,col=adjustcolor("grey",0.5)))
sum.dat=ddply(subset(MC.load.all,WQSite==site.val[k]),"WY",summarise,mean.val=mean(TPLoad.kg))
with(sum.dat,pt_line(WY,mean.val,1,"dodgerblue1",1,21,"dodgerblue1"))
with(subset(TLoad.obs,WQSite==site.val[k]),pt_line(WY,TPLoad.kg,1,"indianred1",1,21,"indianred1"))
if(k%in%c(5,10)){axis_fun(1,xmaj,xmin,xmaj,line=-0.5)}else{axis_fun(1,xmaj,xmin,NA)}
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,paste0("Strucutre: ",discharge.site[k]),cex=0.75)
if(k==6){
  legend("topright",c("Simulated Data","Mean Sim","Observed"),
         lty=NA,lwd=c(0.5,0.1,0.1),
         pch=c(21,21,21),pt.bg=c(NA,"dodgerblue1","indianred1"),
         col=c("grey","black","black"),
         pt.cex=1,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0,yjust=0.5)
}
}
mtext(side=1,line=1,outer=T,"Water Year")
mtext(side=2,line=1.5,outer=T,"TP Load (kg Yr\u207B\u00B9)")
dev.off()

par(family="serif",mar=c(2,1.75,1,2),oma=c(2,2,0.75,2))
layout(matrix(1:10,5,2))

for(k in 1:length(site.val)){
  plot(TP.FWM.ugL~WY,subset(TLoad.obs,WQSite==site.val[k]),xlim=c(2005,2020),ylim=c(0,30),type="n")
  with(subset(TLoad.obs,WQSite==site.val[k]),lines(WY,TP.FWM.ugL,col="red"))
  with(MC.load.all,points(WY,TP.FWM.ugL,pch=21,col=adjustcolor("grey",0.25)))
  test=ddply(subset(MC.load.all,WQSite==site.val[k]),"WY",summarise,mean.val=mean(TP.FWM.ugL))
  with(test,lines(WY,mean.val,col="blue"))
  with(subset(monthmean.load.all,WQSite==site.val[k]),lines(WY,TP.FWM.ugL,col="black"))
  
  mtext(side=3,site.val[k])
}
dev.off()

ddply(TLoad.obs,"WQSite",summarise,cv.val=cv.per(TP.FWM.ugL))
ddply(ddply(MC.load.all,c("WQSite","WY"),summarise,mean.val=mean(TP.FWM.ugL)),"WQSite",summarise,cv.val=cv.per(mean.val))
ddply(monthmean.load.all,"WQSite",summarise,cv.val=cv.per(TP.FWM.ugL))
