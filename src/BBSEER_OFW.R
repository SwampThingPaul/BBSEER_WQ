## 
## BBSEER
##
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
ogrListLayers(paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""))
ogrListLayers(paste(gen.GIS,"/AHED_release/AHED_20171102.gdb",sep=""))

## GIS Data
shore=spTransform(readOGR(paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""),"SFWMD_Shoreline"),utm17)
canal=spTransform(readOGR(paste(gen.GIS,"/SFER_GIS_Geodatabase.gdb",sep=""),"SFWMD_Canals"),utm17)

study.area=spTransform(readOGR("C:/Julian_LaCie/_GitHub/BBSEER_WQ/GIS","BBSEER_PRJBND_09092020"),utm17)
# ofw=spTransform(readOGR("C:/Julian_LaCie/_GitHub/BBSEER_WQ/GIS","OFW"),utm17)
# bbseer.ofw=ofw[study.area,]
# writeOGR(bbseer.ofw,GIS.path,"bbseer_ofw",driver="ESRI Shapefile")
ofw2=spTransform(readOGR("C:/Julian_LaCie/_GitHub/BBSEER_WQ/GIS","bbseer_ofw"),utm17)

# estuary.nnc=spTransform(readOGR("C:/Julian_LaCie/_GISData/FDEP","Estuary_NNC"),utm17)
# bbseer.estuary.nnc=estuary.nnc[gBuffer(study.area,width=5000),]
# writeOGR(bbseer.estuary.nnc,GIS.path,"bbseer_nnc",driver="ESRI Shapefile")
bbseer.estuary.nnc=spTransform(readOGR("C:/Julian_LaCie/_GitHub/BBSEER_WQ/GIS","bbseer_nnc"),utm17)

wmd.mon=spTransform(readOGR(paste0(gen.GIS,"/SFWMD_Monitoring_20200221"),"Environmental_Monitoring_Stations"),utm17)

wq.wmd.mon=subset(wmd.mon,ACTIVITY_S=="Surface Water Grab")
wq.wmd.mon$START_DATE=as.Date(substring(wq.wmd.mon$START_DATE,1,10))
wq.wmd.mon$END_DATE=as.Date(substring(wq.wmd.mon$END_DATE,1,10))

bbseer.mon=wq.wmd.mon[gBuffer(study.area,width=5000),]
# bbseer.mon=subset(bbseer.mon,(END_DATE-START_DATE)>(365*2))
plot(bbseer.mon)
plot(subset(bbseer.mon,START_DATE<as.Date("1980-05-01")))

plot(subset(ofw2,ALT_NAME=="BISCAYNE BAY (CAPE FLORIDA)"),add=T)
subset(ofw2,NAME=="BISCAYNE BAY AQUATIC PRESERVE")@data

library(tmap)
tmap_mode("view")
tm_shape(subset(bbseer.mon,START_DATE<as.Date("1980-05-01")))+tm_dots()+
  tm_shape(subset(ofw2,ALT_NAME=="BISCAYNE BAY (CAPE FLORIDA)"))+tm_polygons()

tm_shape(study.area)+tm_polygons(col="grey",alpha=0.5)+
  tm_shape(ofw2)+tm_polygons(col="indianred1",alpha=0.5)

unique(ofw2@data$NAME)

# plot(study.area)
# plot(subset(ofw2,NAME=="John Pennekamp Coral Reef State Park"),add=T)


# -------------------------------------------------------------------------
dates=date.fun(c("1978-05-01","1982-05-01"))
parameters=data.frame(Test.Number=c(18,21,80,25,61,112,179,20,23,9,8,7),param=c("NOx","TKN","TN","TP","Chla","Chla","Chla","NH4","SRP","SPC","DO","Temp"))


# BBAP (North) ------------------------------------------------------------
# tiff(filename=paste0(plot.path,"map_BBAP_north.tiff"),width=4,height=5,units="in",res=200,compression=c("lzw"),bg="white")
# png(filename=paste0(plot.path,"map_BBAP_north.png"),width=4,height=5,units="in",res=200,bg="white")
par(family="serif",oma=c(0.5,0.5,0.5,0.5),mar=c(0.1,0.1,0.1,0.1))
# layout(matrix(c(1:2),1,2,byrow=T),widths = c(1,0.4))

bbox.lims=bbox(subset(ofw2,ALT_NAME=="BISCAYNE BAY (CAPE FLORIDA)"))
plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.05,asp=1)
plot(canal,add=T,col="lightblue")
plot(subset(ofw2,ALT_NAME=="BISCAYNE BAY (CAPE FLORIDA)"),add=T,col=adjustcolor("grey",0.5),border="red")
plot(subset(bbseer.mon,START_DATE<as.Date("1980-05-01")),add=T,pch=21,bg="dodgerblue1",lwd=0.1)
plot(study.area,add=T,lty=2)
box(lwd=1)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=1,seg.len=4)

# plot(0:1,0:1,ann=F,axes=F,type="n")
legend("bottomright",legend=c("BBSEER Study Area","OFW Boundry", "SFWMD Monitoring"),
       pch=c(NA,22,21),lty=c(2,NA,NA),lwd=c(1,0.1),
       pt.bg=c(NA,adjustcolor("grey",0.5),"dodgerblue1"),col=c("black","red","black"),
       pt.cex=1.25,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=1)
dev.off()

# WQ Data
BBAP.n.dat=DBHYDRO_WQ(dates[1],dates[2],c("SP01","MR01","LR01","BS01"),parameters$Test.Number)
BBAP.n.dat=merge(BBAP.n.dat,parameters,"Test.Number")

BBAP.n.dat.xtab=cast(BBAP.n.dat,Station.ID+Date.EST~param,value="HalfMDL",mean)
BBAP.n.dat.xtab$WY=WY(BBAP.n.dat.xtab$Date.EST)
BBAP.n.dat.xtab$TP.ugL=BBAP.n.dat.xtab$TP*1000

plot(TP.ugL~Date.EST,subset(BBAP.n.dat.xtab,Station.ID=="MR01"))

BBAP.n.TP.GM=ddply(BBAP.n.dat.xtab,c('Station.ID','WY'),summarise,TP.GM=exp(mean(log(TP.ugL),na.rm=T)),N.val=N.obs(TP.ugL))
