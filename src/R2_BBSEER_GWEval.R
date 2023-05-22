## 
## BBSEER
## WQ Eval
## 
##
## Code was compiled by Paul Julian
## contact info: paul.julian@floridadep.gov

## BAD ## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

## Preps DSSRip
options(dss_override_location="C:\\projects\\dssrip\\monolith")
options(dss_config_filename="C:\\projects\\dssrip\\dssrip2.config")
options(dss_default_config="monolith-win-x86_64")
options(dss_allowed_states="untested") 
options(dssrip_debug=T)

## Libraries
#devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);
library(openxlsx)
library(plyr)
library(reshape2)
library(dssrip)

#GIS Libraries
library(rgdal)
library(rgeos)
library(tmap)
library(raster)

library(Hmisc)
library(dunn.test)

library(flextable)
library(magrittr)

library(tmap)
tmap_mode("view")


wd="C:/Julian_LaCie/_GitHub/BBSEER_WQ"

fold.loc=paste0(wd,c("/Plots/WQLoadEval/","/Exports/","/Data/","/src/","/GIS"))
# Folder.Maker(paths);#One and done. Creates folders in working DIR
plot.path=fold.loc[1]
export.path=fold.loc[2]
data.path=fold.loc[3]
GIS.path=fold.loc[5]

gen.GIS="C:/Julian_LaCie/_GISData"
db.path=paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""); 

utm17=CRS("+init=epsg:26917")
wgs84=CRS("+init=epsg:4326")


# GIS ---------------------------------------------------------------------
shore=spTransform(readOGR(paste0(gen.GIS,"/SFER_GIS_Geodatabase.gdb"),"SFWMD_Shoreline"),utm17)
canals=spTransform(readOGR(paste0(gen.GIS,"/SFER_GIS_Geodatabase.gdb"),"SFWMD_Canals"),utm17)
# basins=spTransform(readOGR(paste0(gen.GIS,"/AHED_release/AHED_20171102.gdb"),"WATERSHED"),utm17)
roads=spTransform(readOGR(paste0(gen.GIS,"/FDOT"),"FDOT_Roads_SFWMDClip"),utm17)
evpa=spTransform(readOGR(paste0(gen.GIS,"/SFER_GIS_Geodatabase.gdb"),"EPA_Boundary"),utm17)
enp.shore=spTransform(readOGR(paste0(gen.GIS,"/SFER_GIS_Geodatabase.gdb"),"Shoreline_ENPClip"),utm17)
wcas=spTransform(readOGR(paste0(gen.GIS,"/SFER_GIS_Geodatabase.gdb"),"WCAs"),utm17)
ENP=spTransform(readOGR(paste0(gen.GIS,"/SFER_GIS_Geodatabase.gdb"),"ENP"),utm17)

sfwmd.mon=spTransform(readOGR(paste0(gen.GIS,"/SFWMD_Mointoring_20230302"),"DBHYDRO_SITE_STATION"),utm17)
wmd.struct=spTransform(readOGR(paste0(gen.GIS,"/AHED_release/20230405/AHED.gdb"),"STRUCTURE"),utm17)

bbseer=spTransform(readOGR(GIS.path,"BBSEER_PRJBND_09092020"),utm17)
bbseer.rnd2.struct=spTransform(readOGR(paste0(GIS.path,"/Round2/structures"),"StructRef_All"),utm17)

NPS=spTransform(readOGR(paste0(gen.GIS,"/DOI/NPSBoundary"),"nps_boundary"),utm17)
BNP=subset(NPS,UNIT_CODE=="BISC")

FDEP_AP=spTransform(readOGR(paste0(gen.GIS,"/FDEP"),"AquaticPreserves"),utm17)
nnc=spTransform(readOGR(paste0(gen.GIS,"/FDEP"),"Estuary_NNC"),utm17)

trans.lines=spTransform(readOGR(paste0(export.path,"/GIS"),"TransLines"),utm17)

# plot(crop(wmd.struct,bbseer))

head(sfwmd.mon)
plot(crop(subset(sfwmd.mon,ACTIVITY_S=='Groundwater'),bbseer))



# USGS Data ---------------------------------------------------------------
bbseer.bbox=as.numeric(bbox(spTransform(bbseer,wgs84)))
library(dataRetrieval)

test=parameterCdFile
test[grepl("phosphor",test$srsname ),]
# 00650/653

as.numeric(bbox(spTransform(bbseer,wgs84)))
whatNWISsites(
  bBox = c(-83.0, 36.5, -81.0, 38.5),
  parameterCd = c("00010", "00060"),
  hasDataTypeCd = "dv"
)
sites <- whatNWISsites(
  bBox = c(-80.59, 25.19, -80.12, 25.96)
)
sites=subset(sites,site_tp_cd=="GW")
bbseer.usgs.sites.shp=SpatialPointsDataFrame(sites[,c("dec_long_va","dec_lat_va")],
                                             data=sites,
                                             proj4string=wgs84)

tm_shape(bbseer.usgs.sites.shp)+tm_dots()


# test=readNWISqw(sites$site_no,"00650")
dates=date.fun(c("2012-05-01","2022-05-01"))
test=readWQPqw(paste0("USGS-", sites$site_no), "00650",
               startDate = dates[1],
               endDate = dates[2])
range(test$ActivityStartDate)

# -------------------------------------------------------------------------
## check this dataset https://hub.arcgis.com/datasets/MDC::public-ground-water-samples-result/explore?location=25.404114%2C-79.920724%2C9.84
## downloaded and placed in data folder

# WIN - Groundwater data --------------------------------------------------
## Downloaded from FDEP WIN - 2012 - 2022; groundwater; miami-dade county; general water quality,nutrients and field obs
dat.qual=data.frame(QUALIFIER=c(NA,"!","A","D","E","F","I","R","T","U","*","?","B","H","J","K","L","M","N","O","Q","V","Y","Z"),
                    FATALYN=c("N",rep("N",9),rep("Y",14)))

tmp=list.files(paste0(data.path,"WIN/GW/"))
win.dat1=read.table(paste0(data.path,"WIN/GW/",tmp[1]),skip=8,sep="|",header=T)
win.dat2=read.table(paste0(data.path,"WIN/GW/",tmp[2]),skip=8,sep="|",header=T)
win.dat3=read.table(paste0(data.path,"WIN/GW/",tmp[3]),skip=8,sep="|",header=T)
win.dat4=read.table(paste0(data.path,"WIN/GW/",tmp[4]),skip=8,sep="|",header=T)
win.dat5=read.table(paste0(data.path,"WIN/GW/",tmp[5]),skip=8,sep="|",header=T)
win.dat=rbind(win.dat1,win.dat2,win.dat3,win.dat4,win.dat5)

win.dat$Activity.Start.Date.Time=date.fun(as.character(win.dat$Activity.Start.Date.Time),form="%m/%d/%Y %H:%M:%S")
win.dat$Date.EST=date.fun(win.dat$Activity.Start.Date.Time)
range(win.dat$Date.EST)

win.dat$Station.ID=win.dat$Monitoring.Location.ID

# QA/QC qualifiers 
quals=as.character(unique(win.dat$Value.Qualifier))
spl=strsplit(quals,split="")
quals=data.frame(Value.Qualifier=quals,q1=sapply(spl,"[",1),q2=sapply(spl,"[",2),q3=sapply(spl,"[",3))
quals$Fatal=with(quals,ifelse(q1%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q2%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q3%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER,"Y","N"))

# unique(win.dat$Matrix)

head(win.dat)
win.dat.clean=win.dat
win.dat.clean$HalfMDL=with(win.dat.clean,ifelse(Sample.Collection.Type=="Field Testing-Discrete",DEP.Result.Value.Number,
                                                ifelse(DEP.Result.Value.Number<=DEP.MDL,DEP.MDL/2,DEP.Result.Value.Number)))
win.dat.clean=merge(win.dat.clean,quals[,c("Value.Qualifier","Fatal")],"Value.Qualifier",all.x=T)
win.dat.clean=subset(win.dat.clean,Fatal=="N")

unique(win.dat.clean$DEP.Analyte.Name)
win.parameters=data.frame(DEP.Analyte.Name=c("Nitrogen- Total Kjeldahl", "Turbidity", "Depth to Water (from measuring point)", 
                                             "Dissolved Oxygen Saturation", "Specific Conductance", "Dissolved Oxygen", 
                                             "Hardness- Calculated (CaCO3)", "Sulfate", "Alkalinity (CaCO3)", 
                                             "pH", "Purge Volume", "Nitrogen- Total", "Carbon- Organic", "Temperature, Water", 
                                             "Phosphorus- Total", "Nitrate-Nitrite (N)", "Ammonia (N)", "Residues- Filterable (TDS)", 
                                             "Depth to Water (from land surface elevation)", "Fluoride", "Chloride", 
                                             "Color- True", "Orthophosphate (P)"),
                          param=c("TKN","Turb","Depth.measure",
                                  "DOSat","SPC",'DO',
                                  "Hard","SO4","Alk",
                                  "pH","PurgeVol","TN","OC","Temp",
                                  "TP","NOx","NH4","TDS",
                                  "Depth.land","Fluoride","Cl","Color","SRP"))
win.dat.clean=merge(win.dat.clean,win.parameters,"DEP.Analyte.Name")


win.dat.clean.locs=ddply(win.dat.clean,c("Organization.ID","Monitoring.Location.ID","DEP.Latitude","DEP.Longitude"),summarise,
                         N.vals.allparams=N.obs(HalfMDL),
                         min.date=min(Date.EST),max.date=max(Date.EST))
win.dat.clean.locs.shp=spTransform(SpatialPointsDataFrame(win.dat.clean.locs[,c("DEP.Longitude","DEP.Latitude")],data=win.dat.clean.locs,proj4string=wgs84),utm17)

win.locs2=cbind(data.frame(Station.ID=win.dat.clean.locs.shp$Monitoring.Location.ID),coordinates(win.dat.clean.locs.shp))
colnames(win.locs2)=c("Station.ID","UTMX","UTMY")

plot(win.locs2[,2:3])


tm_shape(bbseer)+tm_borders("red")+
  tm_shape(win.dat.clean.locs.shp)+tm_dots()+
  tm_shape(subset(win.dat.clean.locs.shp,Monitoring.Location.ID=="6490"))+tm_dots("blue")


win.dat.clean.xtab=dcast(win.dat.clean,Station.ID+Date.EST~param,value.var="HalfMDL",mean)
win.dat.clean.xtab$TN=with(win.dat.clean.xtab,TN_Combine(NOx,TKN,TN))
win.dat.clean.xtab$DIN=with(win.dat.clean.xtab,NOx+NH4)

head(win.dat.clean.xtab)
vars=c("Station.ID","Date.EST","TP","SRP","TN","NOx","Depth.measure","Cl","SPC","Temp")
win.dat.clean.xtab2=win.dat.clean.xtab[,vars]
win.dat.clean.xtab2$source="FDEP_WIN"

# DBHYDRO -----------------------------------------------------------------
sfwmd.mon$START_DATE=date.fun(as.Date(sfwmd.mon$START_DATE,"%Y/%m/%d"))
plot(crop(subset(sfwmd.mon,ACTIVITY_S=='Groundwater'&START_DATE>date.fun("2013-05-01")),bbseer))
wmd.gw.bbseer=crop(subset(sfwmd.mon,ACTIVITY_S=='Groundwater'&START_DATE>date.fun("2013-05-01")),bbseer)
wmd.gw.bbseer@data

dates=date.fun(c("2012-05-01","2022-05-01"))
wmd.gw=DBHYDRO_WQ(dates[1],dates[2],wmd.gw.bbseer$STATION,test_number=NULL,matrix ="GW")

unique(wmd.gw$Test.Name)[grepl("DEPTH",unique(wmd.gw$Test.Name))]

head(wmd.gw)
ddply(wmd.gw,c("Test.Number","Test.Name","Units"),summarise,N.val=N.obs(HalfMDL))
ddply(wmd.gw,c("Station.ID","Date.EST"),summarise,
      Depth=mean(as.numeric(Depth),na.rm=T),
      TDepth=mean(as.numeric(T.Depth),na.rm=T),
      Upper.Depth=mean(as.numeric(Upper.Depth),na.rm=T),
      Lower.Depth=mean(as.numeric(Lower.Depth),na.rm=T))

params=data.frame(Test.Number=c(18,21,80,25,23,32,99,9,7),
                  param=c("NOx","TKN","TN","TP","SRP","Cl","Depth.measure","SPC","Temp"))
wmd.gw2=merge(wmd.gw,params,"Test.Number")

wmd.gw2.xtab=dcast(wmd.gw2,Station.ID+Date.EST~param,value.var="HalfMDL",mean)
wmd.gw2.xtab$TN=with(wmd.gw2.xtab,TN_Combine(NOx,TKN,TN))
wmd.gw2.xtab$source="SFWMD"

vars=c("Station.ID","Date.EST","TN","NOx","TP","Cl","SPC","Temp","Depth.measure","source")
gw.dat=rbind(win.dat.clean.xtab2[,vars],wmd.gw2.xtab[,vars])
gw.dat$month=as.numeric(format(gw.dat$Date.EST,"%m"))
gw.dat$CY=as.numeric(format(gw.dat$Date.EST,"%Y"))

plot(TP~month,gw.dat,log="xy")

##
station.invent=ddply(gw.dat,c("Station.ID"),summarise,min.date=min(Date.EST),max.date=max(Date.EST),N.TP=N.obs(TP),N.TN=N.obs(TN))

wmd.locs=cbind(data.frame(Station.ID=wmd.gw.bbseer$STATION),coordinates(wmd.gw.bbseer))
colnames(wmd.locs)=c("Station.ID","UTMX","UTMY")
plot(wmd.locs[,2:3])

win.locs2$source="FDEP_WIN"
wmd.locs$source="SFWMD"

sites.locs=rbind(win.locs2,wmd.locs)
sites.locs.shp=SpatialPointsDataFrame(sites.locs[,c("UTMX","UTMY")],data=sites.locs,proj4string=utm17)
sites.locs2=crop(sites.locs.shp,gBuffer(bbseer,width=100))

tm_shape(bbseer)+tm_borders("red")+
  tm_shape(sites.locs.shp)+tm_dots("source")


station.invent2=subset(station.invent,Station.ID%in%sites.locs2$Station.ID&N.TP>=1)

gw.dat2=subset(gw.dat,Station.ID%in%station.invent2$Station.ID&is.na(TP)==F)
gw.dat2=merge(gw.dat2,sites.locs2@data[,c("Station.ID","UTMX","UTMY")],"Station.ID")
head(gw.dat2)

N.TP.vals=ddply(gw.dat2,"Station.ID",summarise,N.TP=N.obs(TP))
sum(N.TP.vals$N.TP)

tm_shape(bbseer)+tm_borders("red")+
  tm_shape(subset(sites.locs.shp,Station.ID%in%unique(gw.dat2$Station.ID)))+tm_dots("source")

# png(filename=paste0(plot.path,"GW_Data_GAM_Model.png"),width=5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(0.5,0.5,0.5,0.5));
# layout(matrix(1:2,1,2,byrow=F),widths=c(1,0.5))

bbox.lims=bbox(gBuffer(bbseer,width=-100))
plot(shore,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,border=NA)
plot(wcas,add=T,col="darkseagreen2",border=NA)
plot(enp.shore,add=T,col="darkseagreen2",border="grey",lwd=0.1)
plot(ENP,add=T,col=NA,border="white",lty=2,lwd=2)
plot(subset(FDEP_AP,SHORT_NAME=="Biscayne Bay"),add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5)
plot(BNP,add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5,lty=2)
plot(shore,col=NA,border="grey",add=T,lwd=0.05)
plot(canals,add=T,col="lightblue")
tmp.locs=subset(sites.locs.shp,Station.ID%in%unique(gw.dat2$Station.ID))
tmp.locs$source=factor(tmp.locs$source,levels=c("FDEP_WIN","SFWMD"))
cols=c("dodgerblue","violetred1")
plot(tmp.locs,add=T,pch=21,col="grey",bg=cols[tmp.locs$source],lwd=0.1,cex=1.25)
plot(bbseer,add=T,border="red",lty=2)
box(lwd=1)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=0.75,seg.len=4,outer=F);

# plot(0:1,0:1,ann=F,axes=F,type="n")
legend("topleft",legend=c("FDEP WIN","SFWMD DBHYDRO"),
       pch=c(21),lty=0,lwd=c(0.1),
       pt.bg=cols,col="grey",
       pt.cex=1.25,ncol=1,cex=0.75,bty="o",box.col=adjustcolor("grey",0.5),box.lwd=0,y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj=0,title=" Groundwater Data (Sources)")
box(lwd=1)
dev.off()

plot(TP~Cl,gw.dat2,log="xy")
plot(TP~SPC,gw.dat2,log="xy")
plot(TP~Depth.measure,gw.dat2,log="xy")


## GAM
library(mgcv)
library(DHARMa)
library(gratia)

m1a=bam(log(TP)~
          s(month,bs="cc",k=12)+
          s(CY,k=10)+
          s(UTMX,UTMY, bs = 'ds', m = c(1, 0.5),k=30)+
          ti(month, CY,bs = c('cc', 'tp'),k=c(12,10)) +
          ti(UTMX,UTMY, month, d = c(2,1), bs = c('ds','cc'),
             m = list(c(1, 0.5), NA)) +
          ti(UTMX,UTMY, CY, d = c(2,1), bs = c('ds','tp'),
             m = list(c(1, 0.5), NA)),
        data=gw.dat2,nthreads=12,discrete=T,
        method="fREML"
)
testResiduals(simulateResiduals(m1a))
summary(m1a)

# http://r.qcbs.ca/workshop08/book-en/gam-model-checking.html
k.check(m1a)
layout(matrix(1:4,2,2))
gam.check(m1a)

# write.csv(gw.dat2,paste0(export.path,"BBSEER_20230522_GWWQData.csv"),row.names = F)
# save(m1a,file=paste0(export.path,"BBSEER_GW_GAM.RData"))

# draw(m1a)
plot(m1a)

gw.dat2$Pred.TP=exp(predict(m1a,gw.dat2))
plot(TP~Pred.TP,gw.dat2);abline(0,1)

range(gw.dat2$TP)

# png(filename=paste0(plot.path,"TPGAM_obspred.png"),width=5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,1,0.5,0.5),oma=c(3,2.5,0.75,0.25),lwd=0.5);

xlim.val=c(0.001,0.4);xmaj=log.scale.fun(xlim.val,"major");xmin=log.scale.fun(xlim.val,"minor")# by.x=0.1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=xlim.val;ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")# by.y=by.x;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)

plot(TP~Pred.TP,gw.dat2,xlim=xlim.val,ylim=ylim.val,ann=F,axes=F,type="n",log="xy");
abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.1)
points(TP~Pred.TP,gw.dat2,pch=21,bg=adjustcolor("dodgerblue1",0.5),col=adjustcolor("dodgerblue1",0.5),cex=0.75)
abline(0,1,col="indianred1",lty=2,lwd=1.5)
axis_fun(1,xmaj,xmin,xmaj*1000)
axis_fun(2,ymaj,ymin,ymaj*1000);box(lwd=1)
mtext(side=2,line=2.5,"Observed TP (\u03BCg L\u207B\u00B9)")
mtext(side=1,line=2,"Predicted TP (\u03BCg L\u207B\u00B9)")
mtext(side=3,adj=0,line=-2,paste0(" GAM R\u00B2 adj: ",round(summary(m1a)$r.sq,2),"\n POR: May 2012 - May 2022"),cex=0.8)
dev.off()


reg.ext=extent(bbseer)
date.rng=seq(dates[1],dates[2],"1 months")
pdat.sp<-data.frame(expand.grid(
  UTMX=seq(reg.ext[1],reg.ext[2],by=1000),
  UTMY=seq(reg.ext[3],reg.ext[4],by=1000),
  Date.EST=date.rng
))
pdat.sp$month=as.numeric(format(pdat.sp$Date.EST,"%m"))
pdat.sp$CY=as.numeric(format(pdat.sp$Date.EST,"%Y"))

fit.mod <- predict(m1a, pdat.sp,nthreads=12)
pred.mod <- cbind(pdat.sp, Fitted = exp(fit.mod))

mod.fit=predict(m1a,newdata=pdat.sp,type="terms",se.fit = T,nthreads=12)
tmp.fit=mod.fit$fit
colnames(tmp.fit)=paste("fit",c("month","CY","UTMXY","CYmonth","monthUTMXY","CYUTMXY"),sep=".")
tmp.SE=mod.fit$se.fit
colnames(tmp.SE)=paste("SE",c("month","CY","UTMXY","CYmonth","monthUTMXY","CYUTMXY"),sep=".")

pdat.sp=cbind(pdat.sp,tmp.fit,tmp.SE)
head(tmp.fit)
head(pdat.sp)

df.res <- df.residual(m1a)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)

pdat.sp <- transform(pdat.sp,
                     upper.month = fit.month + (crit.t * SE.month),
                     lower.month = fit.month - (crit.t * SE.month),
                     upper.CY = fit.CY + (crit.t * SE.CY),
                     lower.CY = fit.CY - (crit.t * SE.CY))

GAM.month.eff=ddply(pdat.sp,c('month'),summarise,
                  fit=mean(fit.month),
                  UCI=mean(upper.month,na.rm=T),
                  LCI=mean(lower.month,na.rm=T))

pred.org=predict.gam(m1a,type="terms")
partial.resids<-pred.org+residuals(m1a)

head(partial.resids)
range(partial.resids[1])

# png(filename=paste0(plot.path,"TP_GW_GAM.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,1.75,1,1),oma=c(1,2,1,0.25));
layout(matrix(c(1,1,2,2,3,4),2,3,byrow=F),widths=c(1,1,1),heights=c(0.75,0.5))

ylim.val=c(-1,1.5);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(1,12);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
plot(fit~month,GAM.month.eff,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
points(m1a$model$month,partial.resids[,1],pch=19,col=adjustcolor("dodgerblue1",0.5))
abline(h=0)
lines(fit~month,GAM.month.eff,lwd=2)
lines(UCI ~ month, data = GAM.month.eff, lty = "dashed")
lines(LCI ~ month, data = GAM.month.eff, lty = "dashed")
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0,"s(month)")
mtext(side=1,line=1.5,"Month")
mtext(side=2,line=2,"Effect")

ylim.val=c(-1.5,1.5);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(2012,2022);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
plot(fit.CY~CY,pdat.sp,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
points(m1a$model$CY,partial.resids[,2],pch=19,col=adjustcolor("dodgerblue1",0.5))
abline(h=0)
lines(fit.CY~CY,pdat.sp,lwd=2)
lines(upper.CY ~ CY, data = pdat.sp, lty = "dashed")
lines(lower.CY ~ CY, data = pdat.sp, lty = "dashed")
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0,"s(Year)")
mtext(side=1,line=1.5,"Year")

par(mar=c(0.1,0.1,1,0.1))
tmp=pdat.sp[,c("fit.UTMXY","UTMX","UTMY")]
coordinates(tmp)<-~UTMX + UTMY
gridded(tmp)<-TRUE
tmp=as(tmp,"RasterLayer")
proj4string(tmp)<-utm17
# tmp.m=raster::mask(tmp,crop(shore,bbseer))
tmp.m=raster::mask(tmp,bbseer)
tmp.m=raster::mask(tmp.m,crop(shore,bbseer))

bbox.lims=bbox(gBuffer(bbseer,width=-100))
plot(shore,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,border=NA)
plot(wcas,add=T,col="darkseagreen2",border=NA)
plot(enp.shore,add=T,col="darkseagreen2",border="grey",lwd=0.1)
plot(ENP,add=T,col=NA,border="white",lty=2,lwd=2)

breaks.val=seq(-1.75,3,0.25)
pal=viridis::plasma(length(breaks.val)-1,direction = -1)# hcl.colors(length(b)-1, "Spectral")
image(tmp.m,breaks=breaks.val,col=pal,axes=F,ann=F,add=T)

plot(shore,col=NA,border="grey",add=T,lwd=0.05)
plot(bbseer,add=T,border="red",lty=2)
plot(canals,add=T,col="lightblue")
plot(subset(sites.locs.shp,Station.ID%in%station.invent2$Station.ID),add=T,pch=21,bg=NA,col="black",lwd=0.1)
mtext(side=3,adj=0,"s(UTMX,UTMY)")
box(lwd=1)

legend_image=as.raster(matrix(rev(pal),ncol=1))
plot(0:1,0:1,ann=F,axes=F,type="n")
top.val=0.9;bot.val=0.1;mid.v.val=bot.val+(top.val-bot.val)/2
x.max=0.5;x.min=0.3;mid.val=x.min+(x.max-x.min)/2
rasterImage(legend_image,x.min,bot.val,x.max,top.val,xpd=NA)
text(x=x.max, y = c(bot.val,top.val), 
     labels =format(round(range(breaks.val),2)),
     cex=1,adj=0,pos=4,offset=0.5,xpd=NA)
text(x=mid.val,y=top.val,"Effect",pos=3,cex=1,xpd=NA)
dev.off()



# Animation ---------------------------------------------------------------

pdat.sp2=pdat.sp# subset(pdat.sp,month%in%c(1,5,11))

unique(pdat.sp$month)
unique(pdat.sp$CY)

fit.mod <- predict(m1a, pdat.sp2)
pred.mod <- cbind(pdat.sp2, Fitted = exp(fit.mod)*1000)
 
# time.vals=expand.grid(
#   month=1:12,
#   CY=seq(2012,2022,1))
# time.vals=time.vals[2:(nrow(time.vals)-2),]

time.vals=data.frame(date=unique(pdat.sp2$Date.EST))
time.vals$month=as.numeric(format(time.vals$date,"%m"))
time.vals$CY=as.numeric(format(time.vals$date,"%Y"))

time.vals=time.vals[seq(1,nrow(time.vals),3),]

for(i in 1:nrow(time.vals)){
  tmp=subset(pred.mod,CY==time.vals$CY[i]&month==time.vals$month[i])[,c("Fitted","UTMX","UTMY")]
  coordinates(tmp)<- ~UTMX + UTMY
  gridded(tmp)<-TRUE
  # rasterDF<-raster::raster(tmp,layer=1,values=T)
  tmp=as(tmp,"RasterLayer")
  proj4string(tmp)<-utm17
  tmp.m=raster::mask(tmp,bbseer)
  tmp.m=raster::mask(tmp.m,crop(shore,bbseer))
  # assign(paste0("GAM.GW.TP.",time.vals$CY[i],"_",time.vals$month[i]),tmp.m)

png(filename=paste0(plot.path,"GW_animation/","GAM_GW_TP_",time.vals$CY[i],"_",time.vals$month[i],".png"),width=6.5,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(0.5,0.5,0.5,0.5));
layout(matrix(1:2,1,2),widths=c(1,0.5))

bbox.lims=bbox(gBuffer(bbseer,width=-100))
plot(shore,ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,border=NA)
plot(wcas,add=T,col="darkseagreen2",border=NA)
plot(enp.shore,add=T,col="darkseagreen2",border="grey",lwd=0.1)
plot(ENP,add=T,col=NA,border="white",lty=2,lwd=2)

breaks.val=c(seq(0,150,10),2000)
pal=c(viridis::magma(length(breaks.val)-2,direction = -1),"red")
image(tmp.m,breaks=breaks.val,col=pal,axes=F,ann=F,add=T)

plot(shore,col=NA,border="grey",add=T,lwd=0.05)
plot(bbseer,add=T,border="red",lty=2,lwd=1.5)
plot(canals,add=T,col="lightblue")
plot(subset(sites.locs.shp,Station.ID%in%station.invent2$Station.ID),add=T,pch=21,bg=NA,col="black",lwd=0.1)
mtext(side=3,adj=0,line=-1.25,paste(" ",month.abb[time.vals$month[i]]," ",time.vals$CY[i],sep=""))
box(lwd=1)
mapmisc::scaleBar(utm17,"bottomright",bty="n",cex=0.75,seg.len=4,outer=F);

legend_image=as.raster(matrix(rev(pal[1:(length(pal)-1)]),ncol=1))
plot(0:1,0:1,ann=F,axes=F,type="n")
top.val=0.8;bot.val=0.1;mid.v.val=bot.val+(top.val-bot.val)/2
x.max=0.5;x.min=0.3;mid.val=x.min+(x.max-x.min)/2
rasterImage(legend_image,x.min,bot.val,x.max,top.val,xpd=NA)
text(x=x.max, y = c(bot.val,top.val), 
     labels =format(round(range(breaks.val[breaks.val<2000]),2)),
     cex=1,adj=0,pos=4,offset=0.5,xpd=NA)
legend_image=as.raster(matrix(pal[length(pal)],ncol=1))
top.val=0.9;bot.val=0.8;mid.v.val=bot.val+(top.val-bot.val)/2
x.max=0.5;x.min=0.3;mid.val=x.min+(x.max-x.min)/2
rasterImage(legend_image,x.min,bot.val,x.max,top.val,xpd=NA)
text(x=x.max, y = mid.v.val, 
     labels =paste(max(breaks.val[breaks.val<2000]),2000,sep="-"),
     cex=1,adj=0,pos=4,offset=0.5,xpd=NA)


text(x=mid.val,y=top.val,"GW TP (\u03BCg L\u207b\u00B9)",pos=3,cex=1,xpd=NA)
dev.off()
print(i)
}


ray.files=paste0(plot.path,"GW_animation/","GAM_GW_TP_",time.vals$CY,"_",time.vals$month,".png")
gifski::gifski(ray.files, 
               delay = 60 / 100, 
               gif_file  = paste0(plot.path,"GW_animation/BBSEER_GW_TP_GAMPred.gif"),
               loop = T)
