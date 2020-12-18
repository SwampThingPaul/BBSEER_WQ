## 
## BBSEER
## Structure WQ WQBEL-like eval
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
library(rgdal)
library(rgeos)
library(tmap)
library(raster)

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
GIS.path=paths[5]

gen.GIS="C:/Julian_LaCie/_GISData"
db.path=paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""); 

utm17=CRS("+proj=utm +zone=17 +datum=WGS84 +units=m")
wgs84=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

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

# GIS ---------------------------------------------------------------------
ogrListLayers(paste0(gen.GIS,"/SFER_GIS_Geodatabase.gdb"))

study.area=spTransform(readOGR(GIS.path,"BBSEER_PRJBND_09092020"),utm17)
estuary.nnc=spTransform(readOGR(GIS.path,"bbseer_nnc"),utm17)

# wmd.monitoring=spTransform(readOGR(paste0(gen.GIS,"/SFWMD_Monitoring_20200221"),"Environmental_Monitoring_Stations"),utm17)
# wq.mon=subset(wmd.monitoring,ACTIVITY_S=="Surface Water Grab")
# subset(wq.mon,SITE%in%wq.nnc$WQSite)@data
storet.sites=read.table(paste0(data.path,"STORET/Station_Results_dwn20201211_2.txt"),sep="|",header=T)
storet.sites=subset(storet.sites,Org.Name%in%c("Dade Environmental Resource Management (Florida)", "South Florida Water Management District"))
storet.sites$screen=with(storet.sites,ifelse(Station.ID%in%c("MI01","MW01")&County=="DADE",0,1))
storet.sites=subset(storet.sites,screen==1)

spl=sapply(strsplit(as.character(storet.sites$Latitude),split="\\s+"),as.numeric)
storet.sites$Latitude= spl[1,]+spl[2,]/60+spl[3,]/3600
spl=sapply(strsplit(as.character(storet.sites$Longitude),split="\\s+"),as.numeric)
storet.sites$Longitude= (spl[1,]+spl[2,]/60+spl[3,]/3600)*-1
storet.sites=SpatialPointsDataFrame(storet.sites[,c("Longitude","Latitude")],data=storet.sites,proj4string=wgs84)
storet.sites=spTransform(storet.sites,utm17)

shore=spTransform(readOGR(paste0(gen.GIS,"/SFER_GIS_Geodatabase.gdb"),"SFWMD_Shoreline"),utm17)
canals=spTransform(readOGR(paste0(gen.GIS,"/SFER_GIS_Geodatabase.gdb"),"SFWMD_Canals"),utm17)
# basins=spTransform(readOGR(paste0(gen.GIS,"/AHED_release/AHED_20171102.gdb"),"WATERSHED"),utm17)
roads=spTransform(readOGR(paste0(gen.GIS,"/FDOT"),"FDOT_Roads_SFWMDClip"),utm17)
evpa=spTransform(readOGR(paste0(gen.GIS,"/SFER_GIS_Geodatabase.gdb"),"EPA_Boundary"),utm17)
enp.shore=spTransform(readOGR(paste0(gen.GIS,"/SFER_GIS_Geodatabase.gdb"),"Shoreline_ENPClip"),utm17)
wcas=spTransform(readOGR(paste0(gen.GIS,"/SFER_GIS_Geodatabase.gdb"),"WCAs"),utm17)

# png(filename=paste0(plot.path,"NNC_site_map.png"),width=4.25,height=5.5,units="in",res=200,type="windows",bg="white")
# tiff(filename=paste0(plot.path,"NNC_site_map.tiff"),width=4.25,height=5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",oma=c(0.5,0.5,0.5,0.5),mar=c(0.1,0.1,0.1,0.1))
bbox.lims=bbox(subset(estuary.nnc,ESTUARY_RE%in%"ENRH"))
plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05)
plot(subset(estuary.nnc,ESTUARY_RE%in%"ENRH"),add=T,col=adjustcolor("white",0.5),border="grey",lwd=0.1)
plot(subset(estuary.nnc,ESTUARY_SE%in%paste0("ENRH",c(6,3,9,5))),add=T,col="indianred1",border="grey",lwd=0.1)
plot(enp.shore,add=T,col="white",border=NA)
plot(wcas,add=T,col="white",border=NA)
plot(evpa,add=T,border="red",lwd=1.25)
plot(canals,add=T,lwd=1,col="lightblue")
plot(roads,add=T,lwd=0.5,lty=1,col="lightgrey")
plot(storet.sites,add=T,pch=21,bg="dodgerblue1",lwd=0.1)
text(storet.sites,"Station.ID",halo=T,cex=0.75,pos=2)
plot(study.area,add=T,col=NA,lty=2)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=1,seg.len=4)
legend("bottomright",legend=c("EPA","BBSEER Study Area","ENR","ENR - nearshore","Mointoring Locations\n(DERM & SFWMD)"),
       pch=c(22,NA,22,22,21),
       lty=c(NA,2,NA,NA,NA),
       lwd=c(1,1,0.5,0.5,0.1),
       col=c("red","black","grey","grey","black"),
       pt.bg=c("white",NA,adjustcolor("white",0.5),"indianred1","dodgerblue1"),
       pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

box(lwd=1)
dev.off()

# Water Quality -----------------------------------------------------------
## https://floridadep.gov/sites/default/files/62-160_help-document_0.pdf
dat.qual=data.frame(QUALIFIER=c(NA,"!","A","D","E","F","I","R","T","U","*","?","B","H","J","K","L","M","N","O","Q","V","Y","Z"),
                    FATALYN=c("N",rep("N",9),rep("Y",14)))

## STORET data (downloaded 2020-12-11;)
## https://floridadep.gov/dear/watershed-services-program/content/winstoret
storet.dat=read.table(paste0(data.path,"STORET/Water_Quality_Results_dwn20201211_2.txt"),sep="|",header=T)
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
plot(TN~Date.EST,subset(storet.xtab,Station.ID=="MW01"))

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


plot(TP~Date.EST,subset(dat.xtab,Station.ID=="MR01"),xlim=date.fun(c("2000-05-01","2005-05-01")),log="y")

unique(subset(storet.dat,Station.ID=="MR01")$Org.Name)
unique(dat.xtab$Station.ID)

range(dat.xtab$TP,na.rm=T)
par(family="serif",mar=c(2,2,1,0.5),oma=c(2.5,2.5,0.25,0.25));
layout(matrix(1:12,3,4))

xlim.val=dates;xmaj=seq(xlim.val[1],xlim.val[2],"5 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
ylim.val=c(0.5,250);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")# by.y=50;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:nrow(wq.nnc)){
  plot(TP~Date.EST,dat.xtab,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")  
  with(subset(dat.xtab,Station.ID==as.character(wq.nnc$WQSite[i])),pt_line(Date.EST,TP,2,"dodgerblue1",1,21,"dodgerblue1"))
  axis_fun(1,xmaj,xmin,format(xmaj,"%m-%Y"),line=-0.5)
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,adj=0,as.character(wq.nnc$WQSite[i]))

  }

range(dat.xtab$TN,na.rm=T)
par(family="serif",mar=c(2,2,1,0.5),oma=c(2.5,2.5,0.25,0.25));
layout(matrix(1:12,3,4))

xlim.val=dates;xmaj=seq(xlim.val[1],xlim.val[2],"5 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
ylim.val=c(0.05,5);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")# by.y=50;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:nrow(wq.nnc)){
  plot(TN~Date.EST,dat.xtab,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")  
  with(subset(dat.xtab,Station.ID==as.character(wq.nnc$WQSite[i])),pt_line(Date.EST,TN,2,"indianred1",1,21,"indianred1"))
  axis_fun(1,xmaj,xmin,format(xmaj,"%m-%Y"),line=-0.5)
  axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
  mtext(side=3,adj=0,as.character(wq.nnc$WQSite[i]))
  
}

# Discharge ---------------------------------------------------------------
q.dat=data.frame()
for(i in 1:nrow(sites)){
  tmp=DBHYDRO_daily(dates[1],dates[2],sites$DBKEY[i])
  tmp$DBKEY=as.character(sites$DBKEY[i])
  q.dat=rbind(q.dat,tmp)
}

q.dat$WY=WY(q.dat$Date)
q.dat$Date.EST=date.fun(q.dat$Date)
q.dat=merge(q.dat,sites,"DBKEY")

range(q.dat$Data.Value,na.rm=T)
q.dat$Data.Value=with(q.dat,ifelse(Data.Value<0,0,Data.Value));# positive flow only

q.dat.tflow=ddply(q.dat,c("WQSite","Date.EST","WY"),summarise,TFlow.cfs=sum(Data.Value,na.rm=T))

flow.wq=merge(q.dat.tflow,dat.xtab,by.x=c("WQSite","Date.EST","WY"),by.y=c("Station.ID","Date.EST","WY"),all.y=T)
flow.wq=merge(flow.wq,wq.nnc,"WQSite")

plot(TN~Date.EST,subset(flow.wq,WQSite=="MW01"))


# -------------------------------------------------------------------------
# TP
flow.wq$TP.flow=with(flow.wq,ifelse(TFlow.cfs==0|is.na(TFlow.cfs)==T,NA,TP))
TP.samp.size=cast(flow.wq,WQSite+WY~season,value="TP.flow",fun.aggregate = function(x)N.obs(x))
TP.samp.size$TSamp=rowSums(TP.samp.size[,c("A_Wet","B_Dry")],na.rm=T)
TP.samp.size$sea.screen=with(TP.samp.size, ifelse(A_Wet>0&B_Dry>0&TSamp>=4,1,0))
subset(TP.samp.size,sea.screen==0)

# TN
flow.wq$TN.flow=with(flow.wq,ifelse(TFlow.cfs==0|is.na(TFlow.cfs)==T,NA,TN))
TN.samp.size=cast(flow.wq,WQSite+WY~season,value="TN.flow",fun.aggregate = function(x)N.obs(x))
TN.samp.size$TSamp=rowSums(TN.samp.size[,c("A_Wet","B_Dry")],na.rm=T)
TN.samp.size$sea.screen=with(TN.samp.size, ifelse(A_Wet>0&B_Dry>0&TSamp>=4,1,0))
subset(TN.samp.size,sea.screen==0)


TP.annual=ddply(merge(flow.wq,TP.samp.size,c("WQSite","WY")),c("WQSite","WY","sea.screen"),summarise,
      TP.GM.flow=exp(mean(log(TP.flow),na.rm=T)),
      TP.GM.all=exp(mean(log(TP),na.rm=T)),
      TP.FWM=wtd.mean(TP,TFlow.cfs,na.rm=T),N.val=N.obs(TP),N.flow.val=N.obs(ifelse(TFlow.cfs==0|is.na(TFlow.cfs)==T,NA,TP)))
TP.annual=merge(TP.annual,wq.nnc,"WQSite")
plot(TP.FWM~TP.GM.flow,TP.annual);abline(0,1)
plot(TP.FWM~TP.GM.all,TP.annual);abline(0,1)

TN.annual=ddply(merge(flow.wq,TN.samp.size,c("WQSite","WY")),c("WQSite","WY","sea.screen"),summarise,
                TN.GM.flow=exp(mean(log(TN.flow),na.rm=T)),
                TN.GM.all=exp(mean(log(TN),na.rm=T)),
                TN.FWM=wtd.mean(TN,TFlow.cfs,na.rm=T),N.val=N.obs(TN),N.flow.val=N.obs(ifelse(TFlow.cfs==0|is.na(TFlow.cfs)==T,NA,TN)))
TN.annual=merge(TN.annual,wq.nnc,"WQSite")


# NNC Rescale -------------------------------------------------------------

## ENRH6 Rescale
TP.NNC=0.007*1000
TN.NNC=0.48

## TP
H6.TP.annual=subset(TP.annual,NNC_seg=="ENRH6"&sea.screen==1)
boxplot(TP.GM.flow~WQSite,subset(H6.TP.annual,sea.screen==1))

with(subset(H6.TP.annual,sea.screen==1),dunn.test(TP.GM.flow,WQSite))
with(subset(H6.TP.annual,sea.screen==1),dunn.test(TP.FWM,WQSite))

plot(TP.FWM~TP.GM.flow,subset(H6.TP.annual,sea.screen==1))
abline(0,1)
mod=mblm::mblm(TP.FWM~TP.GM.flow,subset(H6.TP.annual,sea.screen==1),repeated = T)
mod=lm(TP.FWM~TP.GM.flow+1,subset(H6.TP.annual,sea.screen==1))
summary(mod)
abline(mod,col="red")

H6.TP.annual$FWM_GM=with(H6.TP.annual,TP.FWM/TP.GM.flow)
H6.TP.rescale=ddply(H6.TP.annual,"WQSite",summarise,mean.GM=mean(TP.GM.flow,na.rm=T),mean.FWM=mean(TP.FWM,na.rm=T),RescaleFac=TP.NNC/mean.GM)
H6.TP.annual=merge(H6.TP.annual,H6.TP.rescale,"WQSite")
H6.TP.annual$FWM.RS=with(H6.TP.annual,TP.FWM*RescaleFac)

H6.TP.sum=ddply(H6.TP.annual,"WQSite",summarise,mean.GM=mean(TP.GM.flow,na.rm=T),mean.FWM=mean(TP.FWM,na.rm=T),sd.lnRSFWM=sd(log(FWM.RS),na.rm=T),N.lnRSFWM=N.obs(FWM.RS))
H6.TP.sum$df=with(H6.TP.sum,N.lnRSFWM-1)
H6.TP.sum$pool.var=with(H6.TP.sum,(sd.lnRSFWM^2)*df)

k.val=N.obs(H6.TP.rescale$RescaleFac)
n.siteyrs=with(H6.TP.annual,N.obs(FWM.RS))
dof=n.siteyrs-k.val
mean.RS=mean(H6.TP.annual$FWM.RS)
SE.RS=SE(H6.TP.annual$FWM.RS)
mean.ln.RS=mean(log(H6.TP.annual$FWM.RS))
sd.ln.RS=sd(log(H6.TP.annual$FWM.RS))
pooled.sd.ln.RS=with(H6.TP.sum,sqrt(sum(pool.var))/sum(df))
prob=0.10
Tp=abs(qt(prob,dof));
ann.FWM.lim=exp(mean.ln.RS+pooled.sd.ln.RS*Tp);ann.FWM.lim
LT.FWM.lim=mean.RS+pooled.sd.ln.RS*(Tp/sqrt(n.siteyrs));LT.FWM.lim

result.tab.H6=data.frame(Parameter=c("Downstream NNC","Number of Sites; k","Number of Site Years; N","Degree of Freedom; df","Mean Rescaled FWM; LTFWM","SE of Rescaled FWM; SE LTFWM","Annual Ln Mean FWM; m","Annual Ln SD FWM; std","Pooled Ln SD FWM; s","Assumed Probability; p","Students-t; Tp","Long Term FWM Limit"))
result.tab.H6=cbind(result.tab.H6,data.frame(Value.TP=c(TP.NNC,k.val,n.siteyrs,dof,mean.RS,SE.RS,mean.ln.RS,sd.ln.RS,pooled.sd.ln.RS,prob,Tp,LT.FWM.lim)))

## TN
H6.TN.annual=subset(TN.annual,NNC_seg=="ENRH6"&sea.screen==1)

boxplot(TN.GM.flow~WQSite,subset(H6.TN.annual,sea.screen==1))
with(subset(H6.TN.annual,sea.screen==1),dunn.test(TN.GM.flow,WQSite))

boxplot(TN.FWM~WQSite,subset(H6.TN.annual,sea.screen==1))
with(subset(H6.TN.annual,sea.screen==1),dunn.test(TN.FWM,WQSite))

plot(TN.FWM~TN.GM.flow,subset(H6.TN.annual,sea.screen==1))
abline(0,1)
mod=mblm::mblm(TN.FWM~TN.GM.flow,H6.TN.annual,repeated = T)
mod=lm(TN.FWM~TN.GM.flow+1,subset(H6.TN.annual,sea.screen==1))
summary(mod)
abline(mod,col="red")

H6.TN.annual$FWM_GM=with(H6.TN.annual,TN.FWM/TN.GM.flow)
H6.TN.rescale=ddply(H6.TN.annual,"WQSite",summarise,mean.GM=mean(TN.GM.flow,na.rm=T),mean.FWM=mean(TN.FWM,na.rm=T),RescaleFac=TN.NNC/mean.GM)
H6.TN.annual=merge(H6.TN.annual,H6.TN.rescale,"WQSite")
H6.TN.annual$FWM.RS=with(H6.TN.annual,TN.FWM*RescaleFac)

H6.TN.sum=ddply(H6.TN.annual,"WQSite",summarise,mean.GM=mean(TN.GM.flow,na.rm=T),mean.FWM=mean(TN.FWM,na.rm=T),sd.lnRSFWM=sd(log(FWM.RS),na.rm=T),N.lnRSFWM=N.obs(FWM.RS))
H6.TN.sum$df=with(H6.TN.sum,N.lnRSFWM-1)
H6.TN.sum$pool.var=with(H6.TN.sum,(sd.lnRSFWM^2)*df)

k.val=N.obs(H6.TN.rescale$RescaleFac)
n.siteyrs=with(H6.TN.annual,N.obs(FWM.RS))
dof=n.siteyrs-k.val
mean.RS=mean(H6.TN.annual$FWM.RS)
SE.RS=SE(H6.TN.annual$FWM.RS)
mean.ln.RS=mean(log(H6.TN.annual$FWM.RS))
sd.ln.RS=sd(log(H6.TN.annual$FWM.RS))
pooled.sd.ln.RS=with(H6.TN.sum,sqrt(sum(pool.var))/sum(df))
prob=0.10
Tp=abs(qt(prob,dof));
ann.FWM.lim=exp(mean.ln.RS+pooled.sd.ln.RS*Tp);ann.FWM.lim
LT.FWM.lim=mean.RS+pooled.sd.ln.RS*(Tp/sqrt(n.siteyrs));LT.FWM.lim

result.tab.H6=cbind(result.tab.H6,data.frame(Value.TN=c(TN.NNC,k.val,n.siteyrs,dof,mean.RS,SE.RS,mean.ln.RS,sd.ln.RS,pooled.sd.ln.RS,prob,Tp,LT.FWM.lim)))
# write.csv(result.tab.H6,paste0(export.path,"H6NNC_rescale.csv"),row.names = F)

result.tab.H6%>%
  flextable()%>%
  colformat_num(i=1:4,j=2,digits=0)%>%
  colformat_num(i=c(5:7,11),j=2,digits=2)%>%
  colformat_num(i=8:9,j=2,digits=3)%>%
  colformat_num(i=10,j=2,digits=1)%>%
  colformat_num(i=12,j=2,digits=0)%>%
  colformat_num(i=1,j=3,digits=2)%>%
  colformat_num(i=2:4,j=3,digits=0)%>%
  colformat_num(i=c(5:7,11),j=3,digits=2)%>%
  colformat_num(i=8:9,j=3,digits=3)%>%
  colformat_num(i=10,j=3,digits=1)%>%
  colformat_num(i=12,j=3,digits=2)%>%
  width(width=c(2.25,1,1))%>%
  add_header(Value.TP="Total\nPhosphorus\n(\u03BCg L\u207B\u00B9)",
             Value.TN="Total\nNitrogen\n(mg L\u207B\u00B9)")%>%
  set_header_labels("Value.TP"="Value","Value.TN"="Value")%>%
  align(j=2:3,align="center",part="all")%>%
  font(fontname = "serif",part="all")%>%
  fontsize(size=10,part="body")%>%
  fontsize(size=12,part="header")%>%
  bold(part="header")


# png(filename=paste0(plot.path,"H6_FWMGM.png"),width=4.25,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1.25,2,1,0.25),oma=c(1.5,1.5,0.25,0.5));
layout(matrix(1:2,2,1))
xlim.val=c(0,15);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,35);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TP.FWM~TP.GM.flow,H6.TP.annual,xlim=xlim.val,ylim=ylim.val,ann=F,axes=F,type="n",xaxs="i",yaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(H6.TP.annual,points(TP.GM.flow,TP.FWM,pch=21,bg="grey",lwd=0.01))
with(H6.TP.rescale,points(mean.GM,mean.FWM,pch=22,bg="dodgerblue1",lwd=0.01))
abline(0,1,lty=2)
# abline(mblm::mblm(TP.FWM~TP.GM.flow,H6.TP.annual),col="indianred1")
axis_fun(1,xmaj,xmin,xmaj,line=-0.6,cex=0.75)
axis_fun(2,ymaj,ymin,ymaj,cex=0.75);box(lwd=1)
mtext(side=1,line=1,"TP AGM (\u03BCg L\u207B\u00B9)",cex=0.8)
mtext(side=2,line=2,"TP FWM (\u03BCg L\u207B\u00B9)",cex=0.8)
mtext(side=3,adj=0,"H6 - South Central Inshore (MW01, MI01, PR01, BL02)",cex=0.75)
# legend("topleft",legend=c("Annual Value","Overall Mean Value","MBLM Line","1:1 Line"),
#        pch=c(21,22,NA,NA),
#        lty=c(NA,NA,1,2),
#        lwd=c(0.01,0.01,1,1),
#        col=c("black","black","indianred1","black"),
#        pt.bg=c("grey","dodgerblue1",NA,NA),
#        pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)  
legend("topleft",legend=c("Annual Value","Overall Mean Value","1:1 Line"),
       pch=c(21,22,NA),
       lty=c(NA,NA,2),
       lwd=c(0.01,0.01,1),
       col=c("black","black","black"),
       pt.bg=c("grey","dodgerblue1",NA),
       pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5) 

xlim.val=c(0,4);by.x=1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,4);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TN.FWM~TN.GM.flow,H6.TN.annual,xlim=xlim.val,ylim=ylim.val,ann=F,axes=F,type="n",xaxs="i",yaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(H6.TN.annual,points(TN.GM.flow,TN.FWM,pch=21,bg="grey",lwd=0.1))
with(H6.TN.rescale,points(mean.GM,mean.FWM,pch=22,bg="dodgerblue1",lwd=0.1))
abline(0,1,lty=2)
# abline(mblm::mblm(TN.FWM~TN.GM.flow,H6.TN.annual),col="indianred1")
axis_fun(1,xmaj,xmin,xmaj,line=-0.6,cex=0.75)
axis_fun(2,ymaj,ymin,ymaj,cex=0.75);box(lwd=1)
mtext(side=1,line=1,"TN AGM (mg L\u207B\u00B9)",cex=0.8)
mtext(side=2,line=2,"TN FWM (mg L\u207B\u00B9)",cex=0.8)

dev.off()

## ENRH3 Rescale
TP.NNC=0.007*1000
TN.NNC=0.31

## TP
H3.TP.annual=subset(TP.annual,NNC_seg=="ENRH3"&sea.screen==1)
boxplot(TP.GM.flow~WQSite,subset(H3.TP.annual,sea.screen==1))

with(subset(H3.TP.annual,sea.screen==1),dunn.test(TP.GM.flow,WQSite))
with(subset(H3.TP.annual,sea.screen==1),dunn.test(TP.FWM,WQSite))

plot(TP.FWM~TP.GM.flow,subset(H3.TP.annual,sea.screen==1))
abline(0,1)
mod=mblm::mblm(TP.FWM~TP.GM.flow,subset(H3.TP.annual,sea.screen==1),repeated = T)
mod=lm(TP.FWM~TP.GM.flow+1,subset(H3.TP.annual,sea.screen==1))
summary(mod)
abline(mod,col="red")

H3.TP.annual$FWM_GM=with(H3.TP.annual,TP.FWM/TP.GM.flow)
H3.TP.rescale=ddply(H3.TP.annual,"WQSite",summarise,mean.GM=mean(TP.GM.flow,na.rm=T),mean.FWM=mean(TP.FWM,na.rm=T),RescaleFac=TP.NNC/mean.GM)
H3.TP.annual=merge(H3.TP.annual,H3.TP.rescale,"WQSite")
H3.TP.annual$FWM.RS=with(H3.TP.annual,TP.FWM*RescaleFac)

H3.TP.sum=ddply(H3.TP.annual,"WQSite",summarise,mean.GM=mean(TP.GM.flow,na.rm=T),mean.FWM=mean(TP.FWM,na.rm=T),sd.lnRSFWM=sd(log(FWM.RS),na.rm=T),N.lnRSFWM=N.obs(FWM.RS))
H3.TP.sum$df=with(H3.TP.sum,N.lnRSFWM-1)
H3.TP.sum$pool.var=with(H3.TP.sum,(sd.lnRSFWM^2)*df)

k.val=N.obs(H3.TP.rescale$RescaleFac)
n.siteyrs=with(H3.TP.annual,N.obs(FWM.RS))
dof=n.siteyrs-k.val
mean.RS=mean(H3.TP.annual$FWM.RS)
SE.RS=SE(H3.TP.annual$FWM.RS)
mean.ln.RS=mean(log(H3.TP.annual$FWM.RS))
sd.ln.RS=sd(log(H3.TP.annual$FWM.RS))
pooled.sd.ln.RS=with(H3.TP.sum,sqrt(sum(pool.var))/sum(df))
prob=0.10
Tp=abs(qt(prob,dof));
ann.FWM.lim=exp(mean.ln.RS+pooled.sd.ln.RS*Tp);ann.FWM.lim
LT.FWM.lim=mean.RS+pooled.sd.ln.RS*(Tp/sqrt(n.siteyrs));LT.FWM.lim

result.tab.H3=data.frame(Parameter=c("Downstream NNC","Number of Sites; k","Number of Site Years; N","Degree of Freedom; df","Mean Rescaled FWM; LTFWM","SE of Rescaled FWM; SE LTFWM","Annual Ln Mean FWM; m","Annual Ln SD FWM; std","Pooled Ln SD FWM; s","Assumed Probability; p","Students-t; Tp","Long Term FWM Limit"))
result.tab.H3=cbind(result.tab.H3,data.frame(Value.TP=c(TP.NNC,k.val,n.siteyrs,dof,mean.RS,SE.RS,mean.ln.RS,sd.ln.RS,pooled.sd.ln.RS,prob,Tp,LT.FWM.lim)))

## TN
H3.TN.annual=subset(TN.annual,NNC_seg=="ENRH3"&sea.screen==1)

boxplot(TN.GM.flow~WQSite,subset(H3.TN.annual,sea.screen==1))
with(subset(H3.TN.annual,sea.screen==1),dunn.test(TN.GM.flow,WQSite))

boxplot(TN.FWM~WQSite,subset(H3.TN.annual,sea.screen==1))
with(subset(H3.TN.annual,sea.screen==1),dunn.test(TN.FWM,WQSite))

plot(TN.FWM~TN.GM.flow,subset(H3.TN.annual,sea.screen==1))
abline(0,1)
mod=mblm::mblm(TN.FWM~TN.GM.flow,H3.TN.annual,repeated = F)
mod=lm(TN.FWM~TN.GM.flow+1,subset(H3.TN.annual,sea.screen==1))
summary(mod)
abline(mod,col="red")

H3.TN.annual$FWM_GM=with(H3.TN.annual,TN.FWM/TN.GM.flow)
H3.TN.rescale=ddply(H3.TN.annual,"WQSite",summarise,mean.GM=mean(TN.GM.flow,na.rm=T),mean.FWM=mean(TN.FWM,na.rm=T),RescaleFac=TN.NNC/mean.GM)
H3.TN.annual=merge(H3.TN.annual,H3.TN.rescale,"WQSite")
H3.TN.annual$FWM.RS=with(H3.TN.annual,TN.FWM*RescaleFac)

H3.TN.sum=ddply(H3.TN.annual,"WQSite",summarise,mean.GM=mean(TN.GM.flow,na.rm=T),mean.FWM=mean(TN.FWM,na.rm=T),sd.lnRSFWM=sd(log(FWM.RS),na.rm=T),N.lnRSFWM=N.obs(FWM.RS))
H3.TN.sum$df=with(H3.TN.sum,N.lnRSFWM-1)
H3.TN.sum$pool.var=with(H3.TN.sum,(sd.lnRSFWM^2)*df)

k.val=N.obs(H3.TN.rescale$RescaleFac)
n.siteyrs=with(H3.TN.annual,N.obs(FWM.RS))
dof=n.siteyrs-k.val
mean.RS=mean(H3.TN.annual$FWM.RS,na.rm=T)
SE.RS=SE(H3.TN.annual$FWM.RS)
mean.ln.RS=mean(log(H3.TN.annual$FWM.RS),na.rm=T)
sd.ln.RS=sd(log(H3.TN.annual$FWM.RS),na.rm=T)
pooled.sd.ln.RS=with(H3.TN.sum,sqrt(sum(pool.var,na.rm=T))/sum(df,na.rm=T))
prob=0.10
Tp=abs(qt(prob,dof));
ann.FWM.lim=exp(mean.ln.RS+pooled.sd.ln.RS*Tp);ann.FWM.lim
LT.FWM.lim=mean.RS+pooled.sd.ln.RS*(Tp/sqrt(n.siteyrs));LT.FWM.lim

result.tab.H3=cbind(result.tab.H3,data.frame(Value.TN=c(TN.NNC,k.val,n.siteyrs,dof,mean.RS,SE.RS,mean.ln.RS,sd.ln.RS,pooled.sd.ln.RS,prob,Tp,LT.FWM.lim)))
# write.csv(result.tab.H3,paste0(export.path,"H3NNC_rescale.csv"),row.names = F)

result.tab.H3%>%
  flextable()%>%
  colformat_num(i=1:4,j=2,digits=0)%>%
  colformat_num(i=c(5:7,11),j=2,digits=2)%>%
  colformat_num(i=8:9,j=2,digits=3)%>%
  colformat_num(i=10,j=2,digits=1)%>%
  colformat_num(i=12,j=2,digits=0)%>%
  colformat_num(i=1,j=3,digits=2)%>%
  colformat_num(i=2:4,j=3,digits=0)%>%
  colformat_num(i=c(5:7,11),j=3,digits=2)%>%
  colformat_num(i=8:9,j=3,digits=3)%>%
  colformat_num(i=10,j=3,digits=1)%>%
  colformat_num(i=12,j=3,digits=2)%>%
  width(width=c(2.25,1,1))%>%
  add_header(Value.TP="Total\nPhosphorus\n(\u03BCg L\u207B\u00B9)",
             Value.TN="Total\nNitrogen\n(mg L\u207B\u00B9)")%>%
  set_header_labels("Value.TP"="Value","Value.TN"="Value")%>%
  align(j=2:3,align="center",part="all")%>%
  font(fontname = "serif",part="all")%>%
  fontsize(size=10,part="body")%>%
  fontsize(size=12,part="header")%>%
  bold(part="header")

# png(filename=paste0(plot.path,"H3_FWMGM.png"),width=4.25,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1.25,2,1,0.25),oma=c(1.5,1.5,0.25,0.5));
layout(matrix(1:2,2,1))
xlim.val=c(0,12);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,22);by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TP.FWM~TP.GM.flow,H3.TP.annual,xlim=xlim.val,ylim=ylim.val,ann=F,axes=F,type="n",xaxs="i",yaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(H3.TP.annual,points(TP.GM.flow,TP.FWM,pch=21,bg="grey",lwd=0.01))
with(H3.TP.rescale,points(mean.GM,mean.FWM,pch=22,bg="dodgerblue1",lwd=0.01))
abline(0,1,lty=2)
# abline(mblm::mblm(TP.FWM~TP.GM.flow,H3.TP.annual),col="indianred1")
axis_fun(1,xmaj,xmin,xmaj,line=-0.6,cex=0.75)
axis_fun(2,ymaj,ymin,ymaj,cex=0.75);box(lwd=1)
mtext(side=1,line=1,"TP AGM (\u03BCg L\u207B\u00B9)",cex=0.8)
mtext(side=2,line=2,"TP FWM (\u03BCg L\u207B\u00B9)",cex=0.8)
mtext(side=3,adj=0,"H3 - North Central Inshore (CD01A, SP01)",cex=0.75)
# legend("topleft",legend=c("Annual Value","Overall Mean Value","MBLM Line","1:1 Line"),
#        pch=c(21,22,NA,NA),
#        lty=c(NA,NA,1,2),
#        lwd=c(0.01,0.01,1,1),
#        col=c("black","black","indianred1","black"),
#        pt.bg=c("grey","dodgerblue1",NA,NA),
#        pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)  
legend("topleft",legend=c("Annual Value","Overall Mean Value","1:1 Line"),
       pch=c(21,22,NA),
       lty=c(NA,NA,2),
       lwd=c(0.01,0.01,1),
       col=c("black","black","black"),
       pt.bg=c("grey","dodgerblue1",NA),
       pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
xlim.val=c(0,1);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,1);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TN.FWM~TN.GM.flow,H3.TN.annual,xlim=xlim.val,ylim=ylim.val,ann=F,axes=F,type="n",xaxs="i",yaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(H3.TN.annual,points(TN.GM.flow,TN.FWM,pch=21,bg="grey",lwd=0.1))
with(H3.TN.rescale,points(mean.GM,mean.FWM,pch=22,bg="dodgerblue1",lwd=0.1))
abline(0,1,lty=2)
# abline(mblm::mblm(TN.FWM~TN.GM.flow,H3.TN.annual),col="indianred1")
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.6,cex=0.75)
axis_fun(2,ymaj,ymin,format(ymaj),cex=0.75);box(lwd=1)
mtext(side=1,line=1,"TN AGM (mg L\u207B\u00B9)",cex=0.8)
mtext(side=2,line=2,"TN FWM (mg L\u207B\u00B9)",cex=0.8)

dev.off()


## ENRH9 Rescale
TP.NNC=0.010*1000
TN.NNC=0.29

## TP
H9.TP.annual=subset(TP.annual,NNC_seg=="ENRH9"&sea.screen==1)
boxplot(TP.GM.flow~WQSite,subset(H9.TP.annual,sea.screen==1))

with(subset(H9.TP.annual,sea.screen==1),dunn.test(TP.GM.flow,WQSite))
with(subset(H9.TP.annual,sea.screen==1),dunn.test(TP.FWM,WQSite))

plot(TP.FWM~TP.GM.flow,subset(H9.TP.annual,sea.screen==1))
abline(0,1)
mod=mblm::mblm(TP.FWM~TP.GM.flow,subset(H9.TP.annual,sea.screen==1),repeated = T)
mod=lm(TP.FWM~TP.GM.flow+1,subset(H9.TP.annual,sea.screen==1))
summary(mod)
abline(mod,col="red")

H9.TP.annual$FWM_GM=with(H9.TP.annual,TP.FWM/TP.GM.flow)
H9.TP.rescale=ddply(H9.TP.annual,"WQSite",summarise,mean.GM=mean(TP.GM.flow,na.rm=T),mean.FWM=mean(TP.FWM,na.rm=T),RescaleFac=TP.NNC/mean.GM)
H9.TP.annual=merge(H9.TP.annual,H9.TP.rescale,"WQSite")
H9.TP.annual$FWM.RS=with(H9.TP.annual,TP.FWM*RescaleFac)

H9.TP.sum=ddply(H9.TP.annual,"WQSite",summarise,mean.GM=mean(TP.GM.flow,na.rm=T),mean.FWM=mean(TP.FWM,na.rm=T),sd.lnRSFWM=sd(log(FWM.RS),na.rm=T),N.lnRSFWM=N.obs(FWM.RS))
H9.TP.sum$df=with(H9.TP.sum,N.lnRSFWM-1)
H9.TP.sum$pool.var=with(H9.TP.sum,(sd.lnRSFWM^2)*df)

k.val=N.obs(H9.TP.rescale$RescaleFac)
n.siteyrs=with(H9.TP.annual,N.obs(FWM.RS))
dof=n.siteyrs-k.val
mean.RS=mean(H9.TP.annual$FWM.RS)
SE.RS=SE(H9.TP.annual$FWM.RS)
mean.ln.RS=mean(log(H9.TP.annual$FWM.RS))
sd.ln.RS=sd(log(H9.TP.annual$FWM.RS))
pooled.sd.ln.RS=with(H9.TP.sum,sqrt(sum(pool.var))/sum(df))
prob=0.10
Tp=abs(qt(prob,dof));
ann.FWM.lim=exp(mean.ln.RS+pooled.sd.ln.RS*Tp);ann.FWM.lim
LT.FWM.lim=mean.RS+pooled.sd.ln.RS*(Tp/sqrt(n.siteyrs));LT.FWM.lim

result.tab.H9=data.frame(Parameter=c("Downstream NNC","Number of Sites; k","Number of Site Years; N","Degree of Freedom; df","Mean Rescaled FWM; LTFWM","SE of Rescaled FWM; SE LTFWM","Annual Ln Mean FWM; m","Annual Ln SD FWM; std","Pooled Ln SD FWM; s","Assumed Probability; p","Students-t; Tp","Long Term FWM Limit"))
result.tab.H9=cbind(result.tab.H9,data.frame(Value.TP=c(TP.NNC,k.val,n.siteyrs,dof,mean.RS,SE.RS,mean.ln.RS,sd.ln.RS,pooled.sd.ln.RS,prob,Tp,LT.FWM.lim)))

## TN
H9.TN.annual=subset(TN.annual,NNC_seg=="ENRH9"&sea.screen==1)

boxplot(TN.GM.flow~WQSite,subset(H9.TN.annual,sea.screen==1))
with(subset(H9.TN.annual,sea.screen==1),dunn.test(TN.GM.flow,WQSite))

boxplot(TN.FWM~WQSite,subset(H9.TN.annual,sea.screen==1))
with(subset(H9.TN.annual,sea.screen==1),dunn.test(TN.FWM,WQSite))

plot(TN.FWM~TN.GM.flow,subset(H9.TN.annual,sea.screen==1))
abline(0,1)
mod=mblm::mblm(TN.FWM~TN.GM.flow,H9.TN.annual,repeated = F)
mod=lm(TN.FWM~TN.GM.flow+1,subset(H9.TN.annual,sea.screen==1))
summary(mod)
abline(mod,col="red")

H9.TN.annual$FWM_GM=with(H9.TN.annual,TN.FWM/TN.GM.flow)
H9.TN.rescale=ddply(H9.TN.annual,"WQSite",summarise,mean.GM=mean(TN.GM.flow,na.rm=T),mean.FWM=mean(TN.FWM,na.rm=T),RescaleFac=TN.NNC/mean.GM)
H9.TN.annual=merge(H9.TN.annual,H9.TN.rescale,"WQSite")
H9.TN.annual$FWM.RS=with(H9.TN.annual,TN.FWM*RescaleFac)

H9.TN.sum=ddply(H9.TN.annual,"WQSite",summarise,mean.GM=mean(TN.GM.flow,na.rm=T),mean.FWM=mean(TN.FWM,na.rm=T),sd.lnRSFWM=sd(log(FWM.RS),na.rm=T),N.lnRSFWM=N.obs(FWM.RS))
H9.TN.sum$df=with(H9.TN.sum,N.lnRSFWM-1)
H9.TN.sum$pool.var=with(H9.TN.sum,(sd.lnRSFWM^2)*df)

k.val=N.obs(H9.TN.rescale$RescaleFac)
n.siteyrs=with(H9.TN.annual,N.obs(FWM.RS))
dof=n.siteyrs-k.val
mean.RS=mean(H9.TN.annual$FWM.RS,na.rm=T)
SE.RS=SE(H9.TN.annual$FWM.RS)
mean.ln.RS=mean(log(H9.TN.annual$FWM.RS),na.rm=T)
sd.ln.RS=sd(log(H9.TN.annual$FWM.RS),na.rm=T)
pooled.sd.ln.RS=with(H9.TN.sum,sqrt(sum(pool.var,na.rm=T))/sum(df,na.rm=T))
prob=0.10
Tp=abs(qt(prob,dof));
ann.FWM.lim=exp(mean.ln.RS+pooled.sd.ln.RS*Tp);ann.FWM.lim
LT.FWM.lim=mean.RS+pooled.sd.ln.RS*(Tp/sqrt(n.siteyrs));LT.FWM.lim

result.tab.H9=cbind(result.tab.H9,data.frame(Value.TN=c(TN.NNC,k.val,n.siteyrs,dof,mean.RS,SE.RS,mean.ln.RS,sd.ln.RS,pooled.sd.ln.RS,prob,Tp,LT.FWM.lim)))
# write.csv(result.tab.H9,paste0(export.path,"H9NNC_rescale.csv"),row.names = F)

result.tab.H9%>%
  flextable()%>%
  colformat_num(i=1:4,j=2,digits=0)%>%
  colformat_num(i=c(5:7,11),j=2,digits=2)%>%
  colformat_num(i=8:9,j=2,digits=3)%>%
  colformat_num(i=10,j=2,digits=1)%>%
  colformat_num(i=12,j=2,digits=0)%>%
  colformat_num(i=1,j=3,digits=2)%>%
  colformat_num(i=2:4,j=3,digits=0)%>%
  colformat_num(i=c(5:7,11),j=3,digits=2)%>%
  colformat_num(i=8:9,j=3,digits=3)%>%
  colformat_num(i=10,j=3,digits=1)%>%
  colformat_num(i=12,j=3,digits=2)%>%
  width(width=c(2.25,1,1))%>%
  add_header(Value.TP="Total\nPhosphorus\n(\u03BCg L\u207B\u00B9)",
             Value.TN="Total\nNitrogen\n(mg L\u207B\u00B9)")%>%
  set_header_labels("Value.TP"="Value","Value.TN"="Value")%>%
  align(j=2:3,align="center",part="all")%>%
  font(fontname = "serif",part="all")%>%
  fontsize(size=10,part="body")%>%
  fontsize(size=12,part="header")%>%
  bold(part="header")

# png(filename=paste0(plot.path,"H9_FWMGM.png"),width=4.25,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1.25,2,1,0.25),oma=c(1.5,1.5,0.25,0.5));
layout(matrix(1:2,2,1))
xlim.val=c(0,20);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,30);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TP.FWM~TP.GM.flow,H9.TP.annual,xlim=xlim.val,ylim=ylim.val,ann=F,axes=F,type="n",xaxs="i",yaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(H9.TP.annual,points(TP.GM.flow,TP.FWM,pch=21,bg="grey",lwd=0.01))
with(H9.TP.rescale,points(mean.GM,mean.FWM,pch=22,bg="dodgerblue1",lwd=0.01))
abline(0,1,lty=2)
# abline(mblm::mblm(TP.FWM~TP.GM.flow,H9.TP.annual),col="indianred1")
axis_fun(1,xmaj,xmin,xmaj,line=-0.6,cex=0.75)
axis_fun(2,ymaj,ymin,ymaj,cex=0.75);box(lwd=1)
mtext(side=1,line=1,"TP AGM (\u03BCg L\u207B\u00B9)",cex=0.8)
mtext(side=2,line=2,"TP FWM (\u03BCg L\u207B\u00B9)",cex=0.8)
mtext(side=3,adj=0,"H9 - Southern North Bay (MR01, LR01)",cex=0.75)
# legend("topleft",legend=c("Annual Value","Overall Mean Value","MBLM Line","1:1 Line"),
#        pch=c(21,22,NA,NA),
#        lty=c(NA,NA,1,2),
#        lwd=c(0.01,0.01,1,1),
#        col=c("black","black","indianred1","black"),
#        pt.bg=c("grey","dodgerblue1",NA,NA),
#        pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)  
legend("topleft",legend=c("Annual Value","Overall Mean Value","1:1 Line"),
       pch=c(21,22,NA),
       lty=c(NA,NA,2),
       lwd=c(0.01,0.01,1),
       col=c("black","black","black"),
       pt.bg=c("grey","dodgerblue1",NA),
       pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

xlim.val=c(0,1);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,1.5);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TN.FWM~TN.GM.flow,H9.TN.annual,xlim=xlim.val,ylim=ylim.val,ann=F,axes=F,type="n",xaxs="i",yaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(H9.TN.annual,points(TN.GM.flow,TN.FWM,pch=21,bg="grey",lwd=0.1))
with(H9.TN.rescale,points(mean.GM,mean.FWM,pch=22,bg="dodgerblue1",lwd=0.1))
abline(0,1,lty=2)
# abline(mblm::mblm(TN.FWM~TN.GM.flow,H9.TN.annual),col="indianred1")
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.6,cex=0.75)
axis_fun(2,ymaj,ymin,format(ymaj),cex=0.75);box(lwd=1)
mtext(side=1,line=1,"TN AGM (mg L\u207B\u00B9)",cex=0.8)
mtext(side=2,line=2,"TN FWM (mg L\u207B\u00B9)",cex=0.8)

dev.off()


## ENRH5 Rescale
TP.NNC=0.012*1000
TN.NNC=0.30

## TP
H5.TP.annual=subset(TP.annual,NNC_seg=="ENRH5"&sea.screen==1)
boxplot(TP.GM.flow~WQSite,subset(H5.TP.annual,sea.screen==1))

with(subset(H5.TP.annual,sea.screen==1),dunn.test(TP.GM.flow,WQSite))
with(subset(H5.TP.annual,sea.screen==1),dunn.test(TP.FWM,WQSite))

plot(TP.FWM~TP.GM.flow,subset(H5.TP.annual,sea.screen==1))
abline(0,1)
mod=mblm::mblm(TP.FWM~TP.GM.flow,subset(H5.TP.annual,sea.screen==1),repeated = T)
mod=lm(TP.FWM~TP.GM.flow+1,subset(H5.TP.annual,sea.screen==1))
summary(mod)
abline(mod,col="red")

H5.TP.annual$FWM_GM=with(H5.TP.annual,TP.FWM/TP.GM.flow)
H5.TP.rescale=ddply(H5.TP.annual,"WQSite",summarise,mean.GM=mean(TP.GM.flow,na.rm=T),mean.FWM=mean(TP.FWM,na.rm=T),RescaleFac=TP.NNC/mean.GM)
H5.TP.annual=merge(H5.TP.annual,H5.TP.rescale,"WQSite")
H5.TP.annual$FWM.RS=with(H5.TP.annual,TP.FWM*RescaleFac)

H5.TP.sum=ddply(H5.TP.annual,"WQSite",summarise,mean.GM=mean(TP.GM.flow,na.rm=T),mean.FWM=mean(TP.FWM,na.rm=T),sd.lnRSFWM=sd(log(FWM.RS),na.rm=T),N.lnRSFWM=N.obs(FWM.RS))
H5.TP.sum$df=with(H5.TP.sum,N.lnRSFWM-1)
H5.TP.sum$pool.var=with(H5.TP.sum,(sd.lnRSFWM^2)*df)

k.val=N.obs(H5.TP.rescale$RescaleFac)
n.siteyrs=with(H5.TP.annual,N.obs(FWM.RS))
dof=n.siteyrs-k.val
mean.RS=mean(H5.TP.annual$FWM.RS)
SE.RS=SE(H5.TP.annual$FWM.RS)
mean.ln.RS=mean(log(H5.TP.annual$FWM.RS))
sd.ln.RS=sd(log(H5.TP.annual$FWM.RS))
pooled.sd.ln.RS=with(H5.TP.sum,sqrt(sum(pool.var))/sum(df))
prob=0.10
Tp=abs(qt(prob,dof));
ann.FWM.lim=exp(mean.ln.RS+pooled.sd.ln.RS*Tp);ann.FWM.lim
LT.FWM.lim=mean.RS+pooled.sd.ln.RS*(Tp/sqrt(n.siteyrs));LT.FWM.lim

result.tab.H5=data.frame(Parameter=c("Downstream NNC","Number of Sites; k","Number of Site Years; N","Degree of Freedom; df","Mean Rescaled FWM; LTFWM","SE of Rescaled FWM; SE LTFWM","Annual Ln Mean FWM; m","Annual Ln SD FWM; std","Pooled Ln SD FWM; s","Assumed Probability; p","Students-t; Tp","Long Term FWM Limit"))
result.tab.H5=cbind(result.tab.H5,data.frame(Value.TP=c(TP.NNC,k.val,n.siteyrs,dof,mean.RS,SE.RS,mean.ln.RS,sd.ln.RS,pooled.sd.ln.RS,prob,Tp,LT.FWM.lim)))

## TN
H5.TN.annual=subset(TN.annual,NNC_seg=="ENRH5"&sea.screen==1)

boxplot(TN.GM.flow~WQSite,subset(H5.TN.annual,sea.screen==1))
with(subset(H5.TN.annual,sea.screen==1),dunn.test(TN.GM.flow,WQSite))

boxplot(TN.FWM~WQSite,subset(H5.TN.annual,sea.screen==1))
with(subset(H5.TN.annual,sea.screen==1),dunn.test(TN.FWM,WQSite))

plot(TN.FWM~TN.GM.flow,subset(H5.TN.annual,sea.screen==1))
abline(0,1)
mod=mblm::mblm(TN.FWM~TN.GM.flow,H5.TN.annual,repeated = F)
mod=lm(TN.FWM~TN.GM.flow+1,subset(H5.TN.annual,sea.screen==1))
summary(mod)
abline(mod,col="red")

H5.TN.annual$FWM_GM=with(H5.TN.annual,TN.FWM/TN.GM.flow)
H5.TN.rescale=ddply(H5.TN.annual,"WQSite",summarise,mean.GM=mean(TN.GM.flow,na.rm=T),mean.FWM=mean(TN.FWM,na.rm=T),RescaleFac=TN.NNC/mean.GM)
H5.TN.annual=merge(H5.TN.annual,H5.TN.rescale,"WQSite")
H5.TN.annual$FWM.RS=with(H5.TN.annual,TN.FWM*RescaleFac)

H5.TN.sum=ddply(H5.TN.annual,"WQSite",summarise,mean.GM=mean(TN.GM.flow,na.rm=T),mean.FWM=mean(TN.FWM,na.rm=T),sd.lnRSFWM=sd(log(FWM.RS),na.rm=T),N.lnRSFWM=N.obs(FWM.RS))
H5.TN.sum$df=with(H5.TN.sum,N.lnRSFWM-1)
H5.TN.sum$pool.var=with(H5.TN.sum,(sd.lnRSFWM^2)*df)

k.val=N.obs(H5.TN.rescale$RescaleFac)
n.siteyrs=with(H5.TN.annual,N.obs(FWM.RS))
dof=n.siteyrs-k.val
mean.RS=mean(H5.TN.annual$FWM.RS,na.rm=T)
SE.RS=SE(H5.TN.annual$FWM.RS)
mean.ln.RS=mean(log(H5.TN.annual$FWM.RS),na.rm=T)
sd.ln.RS=sd(log(H5.TN.annual$FWM.RS),na.rm=T)
pooled.sd.ln.RS=with(H5.TN.sum,sqrt(sum(pool.var,na.rm=T))/sum(df,na.rm=T))
prob=0.10
Tp=abs(qt(prob,dof));
ann.FWM.lim=exp(mean.ln.RS+pooled.sd.ln.RS*Tp);ann.FWM.lim
LT.FWM.lim=mean.RS+pooled.sd.ln.RS*(Tp/sqrt(n.siteyrs));LT.FWM.lim

result.tab.H5=cbind(result.tab.H5,data.frame(Value.TN=c(TN.NNC,k.val,n.siteyrs,dof,mean.RS,SE.RS,mean.ln.RS,sd.ln.RS,pooled.sd.ln.RS,prob,Tp,LT.FWM.lim)))
# write.csv(result.tab.H5,paste0(export.path,"H5NNC_rescale.csv"),row.names = F)

result.tab.H5%>%
  flextable()%>%
  colformat_num(i=1:4,j=2,digits=0)%>%
  colformat_num(i=c(5:7,11),j=2,digits=2)%>%
  colformat_num(i=8:9,j=2,digits=3)%>%
  colformat_num(i=10,j=2,digits=1)%>%
  colformat_num(i=12,j=2,digits=0)%>%
  colformat_num(i=1,j=3,digits=2)%>%
  colformat_num(i=2:4,j=3,digits=0)%>%
  colformat_num(i=c(5:7,11),j=3,digits=2)%>%
  colformat_num(i=8:9,j=3,digits=3)%>%
  colformat_num(i=10,j=3,digits=1)%>%
  colformat_num(i=12,j=3,digits=2)%>%
  width(width=c(2.25,1,1))%>%
  add_header(Value.TP="Total\nPhosphorus\n(\u03BCg L\u207B\u00B9)",
             Value.TN="Total\nNitrogen\n(mg L\u207B\u00B9)")%>%
  set_header_labels("Value.TP"="Value","Value.TN"="Value")%>%
  align(j=2:3,align="center",part="all")%>%
  font(fontname = "serif",part="all")%>%
  fontsize(size=10,part="body")%>%
  fontsize(size=12,part="header")%>%
  bold(part="header")

# png(filename=paste0(plot.path,"H5_FWMGM.png"),width=4.25,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1.25,2,1,0.25),oma=c(1.5,1.5,0.25,0.5));
layout(matrix(1:2,2,1))
xlim.val=c(0,25);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,30);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TP.FWM~TP.GM.flow,H5.TP.annual,xlim=xlim.val,ylim=ylim.val,ann=F,axes=F,type="n",xaxs="i",yaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(H5.TP.annual,points(TP.GM.flow,TP.FWM,pch=21,bg="grey",lwd=0.01))
with(H5.TP.rescale,points(mean.GM,mean.FWM,pch=22,bg="dodgerblue1",lwd=0.01))
abline(0,1,lty=2)
# abline(mblm::mblm(TP.FWM~TP.GM.flow,H5.TP.annual),col="indianred1")
axis_fun(1,xmaj,xmin,xmaj,line=-0.6,cex=0.75)
axis_fun(2,ymaj,ymin,ymaj,cex=0.75);box(lwd=1)
mtext(side=1,line=1,"TP AGM (\u03BCg L\u207B\u00B9)",cex=0.8)
mtext(side=2,line=2,"TP FWM (\u03BCg L\u207B\u00B9)",cex=0.8)
mtext(side=3,adj=0,"H5 - Northern North Bay (BS01,SK01)",cex=0.75)
# legend("topleft",legend=c("Annual Value","Overall Mean Value","MBLM Line","1:1 Line"),
#        pch=c(21,22,NA,NA),
#        lty=c(NA,NA,1,2),
#        lwd=c(0.01,0.01,1,1),
#        col=c("black","black","indianred1","black"),
#        pt.bg=c("grey","dodgerblue1",NA,NA),
#        pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)  
legend("topleft",legend=c("Annual Value","Overall Mean Value","1:1 Line"),
       pch=c(21,22,NA),
       lty=c(NA,NA,2),
       lwd=c(0.01,0.01,1),
       col=c("black","black","black"),
       pt.bg=c("grey","dodgerblue1",NA),
       pt.cex=1,ncol=1,cex=0.6,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

xlim.val=c(0,0.4);by.x=0.1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,0.4);by.y=0.1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(TN.FWM~TN.GM.flow,H5.TN.annual,xlim=xlim.val,ylim=ylim.val,ann=F,axes=F,type="n",xaxs="i",yaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(H5.TN.annual,points(TN.GM.flow,TN.FWM,pch=21,bg="grey",lwd=0.1))
with(H5.TN.rescale,points(mean.GM,mean.FWM,pch=22,bg="dodgerblue1",lwd=0.1))
abline(0,1,lty=2)
# abline(mblm::mblm(TN.FWM~TN.GM.flow,H5.TN.annual),col="indianred1")
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.6,cex=0.75)
axis_fun(2,ymaj,ymin,format(ymaj),cex=0.75);box(lwd=1)
mtext(side=1,line=1,"TN AGM (mg L\u207B\u00B9)",cex=0.8)
mtext(side=2,line=2,"TN FWM (mg L\u207B\u00B9)",cex=0.8)

dev.off()


