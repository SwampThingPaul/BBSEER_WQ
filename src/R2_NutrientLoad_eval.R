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
# -------------------------------------------------------------------------
struct.wq.list=read.xlsx(paste0(data.path,"BBSEER_Structure_WQ_Region.xlsx"))

struct.wq.list$Station.ID=struct.wq.list$WQ_up
struct.wq.list
struct.wq.list$Station.ID=with(struct.wq.list,ifelse(STRUCTURE=="S21A","PR03",Station.ID))
struct.wq.list$Station.ID=with(struct.wq.list,ifelse(STRUCTURE=="S703","PR03",Station.ID))
struct.wq.list$Station.ID=with(struct.wq.list,ifelse(STRUCTURE=="S20F","MW04",Station.ID))
struct.wq.list$Station.ID=with(struct.wq.list,ifelse(STRUCTURE%in%c("S23A","S23B","S706A"),"PR03",Station.ID))
struct.wq.list$Station.ID=with(struct.wq.list,ifelse(STRUCTURE%in%c("S706B","S706C"),"MI02",Station.ID))
struct.wq.list$Station.ID=with(struct.wq.list,ifelse(STRUCTURE%in%c("S708","S23C"),"MI02",Station.ID))
struct.wq.list$Station.ID=with(struct.wq.list,ifelse(STRUCTURE%in%c("S23D","S23E","S23F","S712A","S712B"),"MW04",Station.ID))
struct.wq.list

struct.wq.list$RSM.site=struct.wq.list$STRUCTURE
struct.wq.list$RSM.site=with(struct.wq.list,ifelse(STRUCTURE%in%c("S23A","S23B"),"S23A-S23B",RSM.site))
struct.wq.list$RSM.site=with(struct.wq.list,ifelse(STRUCTURE%in%c("S23D","S23E","S23F"),"S23D-S23E-S23F",RSM.site))
struct.wq.list$RSM.site=with(struct.wq.list,ifelse(STRUCTURE%in%c("S712A","S712B"),"S712A-S712B",RSM.site))
struct.wq.list$RSM.site=with(struct.wq.list,ifelse(STRUCTURE%in%c("S706B","S706C"),"S706B-S706C",RSM.site))
struct.wq.list$RSM.site=with(struct.wq.list,ifelse(STRUCTURE%in%c("S708","S23C"),"S708-S23C",RSM.site))

S700.rep=data.frame()
for(i in 1:4){
  tmp=subset(struct.wq.list,STRUCTURE=="S700")
  tmp$RSM.site=paste0("S700_P",LETTERS[i])
  S700.rep=rbind(S700.rep,tmp)
}
S700.rep

struct.wq.list=rbind(
  subset(struct.wq.list,STRUCTURE!="S700"),
  S700.rep
)

# write.xlsx(struct.wq.list,paste0(data.path,"BBSEER_Structure_WQ_Region_v2.xlsx"))
## S703 is downstream of S21A - remove S21A from analysis?

# General North to South order
RSM.site.order=c("S29", "G58", "S28", "S27", "S26", "S25B", "S25", "G93", "S22", "S700",
  "S123", "S21", "S703", "S23A-S23B", "S706A", "S706B-S706C", 
  "S20G", "S708-S23C", "S23D-S23E-S23F", "S20F", "S712A-S712B", 
  "S197")

# png(filename=paste0(plot.path,"WQ_BBSEER_structmap.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(1,1,1,1));
layout(matrix(1:2,1,2),widths=c(1,0.75))

bbox.lims=bbox(bbseer)
plot(shore,col="cornsilk",border="grey",bg="lightblue",
     ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(wcas,add=T,col="darkseagreen2",border=NA)
plot(enp.shore,add=T,col="darkseagreen2",border="grey",lwd=0.1)
plot(ENP,add=T,col=NA,border="white",lty=2,lwd=2)
plot(subset(FDEP_AP,SHORT_NAME=="Biscayne Bay"),add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5)
plot(BNP,add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5,lty=2)
plot(canals,add=T,col="lightblue")
plot(subset(bbseer.rnd2.struct,Structure %in%RSM.site.order),add=T,pch=21,bg="red",lwd=0.1)
plot(subset(wmd.struct,NAME%in%RSM.site.order),add=T,pch=21,bg="dodgerblue1",lwd=0.1)
plot(bbseer,add=T,border="red",lty=2)
plot(subset(trans.lines,transect%in%c("ST01_ol","ST02_ol","ST03_ol")),
     col=rev(wesanderson::wes_palette("Zissou1",3,"continuous")),lwd=3,add=T)

box(lwd=1)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=0.75,seg.len=4,outer=F);

plot(0:1,0:1,ann=F,axes=F,type="n")
legend("center",legend=c("Existing Structures","BBSEER Round 2 Structures",paste0("ST0",1:3)),
       lty=c(0,0,1,1,1),lwd=c(0.1,0.1,2,2,2),col=c(rep("black",2),rev(wesanderson::wes_palette("Zissou1",3,"continuous"))),
       pch=c(rep(21,2),rep(NA,3)),pt.bg=c("dodgerblue1","red"),pt.cex=1.5,
       ncol=1,cex=0.75,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj = 0,title="Coastal Water\nControl Structures\n(to Biscayne Bay)")

dev.off()

# png(filename=paste0(plot.path,"WQ_BBSEER_structmap_NNC.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(1,1,1,1));
layout(matrix(1:2,1,2),widths=c(1,0.75))

bbox.lims=bbox(bbseer)
plot(shore,col="cornsilk",border="grey",bg="lightblue",
     ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(wcas,add=T,col="darkseagreen2",border=NA)
plot(enp.shore,add=T,col="darkseagreen2",border="grey",lwd=0.1)
plot(ENP,add=T,col=NA,border="white",lty=2,lwd=2)
# plot(subset(FDEP_AP,SHORT_NAME=="Biscayne Bay"),add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5)
# plot(BNP,add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5,lty=2)
plot(canals,add=T,col="lightblue")
plot(bbseer,add=T,border="red",lty=2)
plot(nnc,add=T,lty=2)
plot(subset(nnc,ESTUARY_SE%in%paste0("ENRH",c(5,9,3,6,2))),add=T,col="grey80",border="white",lty=2)
raster::text(subset(nnc,ESTUARY_SE%in%paste0("ENRH",c(5,9,3,6,2))),subset(nnc,ESTUARY_SE%in%paste0("ENRH",c(5,9,3,6,2)))$ESTUARY_SE,halo=T,col="red",font=2,cex=0.75,pos=4)
plot(subset(bbseer.rnd2.struct,Structure %in%RSM.site.order),add=T,pch=21,bg="red",lwd=0.1)
plot(subset(wmd.struct,NAME%in%RSM.site.order),add=T,pch=21,bg="dodgerblue1",lwd=0.1)

box(lwd=1)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=0.75,seg.len=4,outer=F);

plot(0:1,0:1,ann=F,axes=F,type="n")
legend("center",legend=c("Existing Structures","BBSEER Round 2 Structures","NNC Segments"),
       lty=c(0),lwd=c(0.1),col="black",
       pch=c(rep(21,2),22),pt.bg=c("dodgerblue1","red","grey80"),pt.cex=1.5,
       ncol=1,cex=0.75,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj = 0,title="Coastal Water\nControl Structures\n(to Biscayne Bay)")

dev.off()


# -------------------------------------------------------------------------
## WQ data collated and combined in BBSEER_WQ_Eval.R file
wq.dat=read.csv(paste0(export.path,"/20230504_BBSEERWQ_data.csv"))
wq.dat$Date.EST=date.fun(wq.dat$Date.EST)

wq.dat$month=as.numeric(format(wq.dat$Date.EST,'%m'))
wq.dat$CY=as.numeric(format(wq.dat$Date.EST,'%Y'))
wq.dat$WY=WY(wq.dat$Date.EST)
tmp=ddply(wq.dat,"Station.ID",summarise,min.date=min(Date.EST),max.date=max(Date.EST))
subset(tmp,Station.ID%in%unique(struct.wq.list$Station.ID))

## Last 10 years of data
wq.dat2=subset(wq.dat,WY%in%seq(2013,2022,1))

vars=c("Station.ID","Date.EST","WY","month","CY","TN","TP")
wq.dat2.melt=melt(wq.dat2[,vars],id.vars=vars[!(vars%in%c("TN","TP"))])

month.mean.val=ddply(subset(wq.dat2.melt,Station.ID%in%c(unique(struct.wq.list$Station.ID))),
      c("Station.ID","variable","month"),summarise,
      mean.val=mean(value,na.rm=T),
      N.val=N.obs(value,na.rm=T),
      SD.val=sd(value,na.rm=T),
      CV.val=cv.per(value))
# write.csv(month.mean.val,paste0(export.path,"STRUCT_MeanMonthlyWQ.csv"),row.names = F)



dcast(subset(month.mean.val,variable=="TP"),Station.ID~month,value.var = "mean.val",mean)
dcast(subset(month.mean.val,variable=="TP"),Station.ID~month,value.var = "N.val",mean)

dcast(subset(month.mean.val,variable=="TN"),Station.ID~month,value.var = "mean.val",mean)
dcast(subset(month.mean.val,variable=="TN"),Station.ID~month,value.var = "N.val",mean)


plot(mean.val~SD.val,subset(month.mean.val,variable=="TP"));abline(0,1)
boxplot(subset(month.mean.val,variable=="TP")$CV.val)
plot(mean.val~SD.val,subset(month.mean.val,variable=="TN"));abline(0,1)
boxplot(subset(month.mean.val,variable=="TN")$CV.val)


# subset(month.mean.val,Station.ID=="AC03")
# zoo::na.approx(subset(month.mean.val,Station.ID=="AC03")$mean.val)

## Linearly interpolated missing monthly mean data.
month.mean.val$mean.val=with(month.mean.val,ifelse(Station.ID=="AC03"&variable=="TN",zoo::na.approx(mean.val),mean.val))
month.mean.val.xtab=dcast(month.mean.val,Station.ID+month~variable,value.var = "mean.val",mean)


# png(filename=paste0(plot.path,"WQ_monthlymean_TP.png"),width=6.5,height=5.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,1,0.5,0.5),oma=c(2,3,1,0.5),lwd=0.5);
layout(matrix(c(1:16),4,4,byrow=T))

ylim.val=c(0,0.1);by.y=0.020;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
ylim.val=c(0.001,0.2);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
xlim.val=c(1,12);by.x=3;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
wq.sites=unique(month.mean.val$Station.ID)
for(i in 1:16){
  tmp=subset(month.mean.val,Station.ID==wq.sites[i]&variable=="TP")
  plot(mean.val~month,tmp,xlim=xlim.val,ylim=ylim.val,ann=F,axes=F,type="n",log="y")
  abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
  with(tmp,pt_line(month,mean.val,2,"indianred1",1,21,"indianred1",cex=1.25))
  # with(tmp,pt_line_error(month,mean.val,SD.val,2,"indianred1",1,21,"indianred1",cex=1.25,length=0.025))
  
  if(i%in%c(1,5,9,13)){axis_fun(2,ymaj,ymin,format(ymaj*1000))}else{axis_fun(2,ymaj,ymin,NA)}
  if(i%in%c(13:16)){axis_fun(1,xmaj,xmin,month.abb[xmaj],line=-0.5)}else{axis_fun(1,xmaj,xmin,NA)}
  mtext(side=3,adj=0,line=-1.2,paste0(" ",wq.sites[i]),cex=0.7)
  box(lwd=1)
  
}
mtext(side=1,outer=T,line=1,"Month")
mtext(side=2,outer=T,line=1.5,"Arithmetic Mean TP (\u03BCg L\u207B\u00B9) ")

dev.off()

# png(filename=paste0(plot.path,"WQ_monthlymean_TN.png"),width=6.5,height=5.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,1,0.5,0.5),oma=c(2,3,1,0.5),lwd=0.5);
layout(matrix(c(1:16),4,4,byrow=T))

ylim.val=c(0,6);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
# ylim.val=c(0.001,0.2);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
xlim.val=c(1,12);by.x=3;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
wq.sites=unique(month.mean.val$Station.ID)
for(i in 1:16){
  tmp=subset(month.mean.val,Station.ID==wq.sites[i]&variable=="TN")
  plot(mean.val~month,tmp,xlim=xlim.val,ylim=ylim.val,ann=F,axes=F,type="n")
  abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
  with(tmp,pt_line(month,mean.val,2,"dodgerblue1",1,21,"dodgerblue1",cex=1.25))
  # with(tmp,pt_line_error(month,mean.val,SD.val,2,"indianred1",1,21,"indianred1",cex=1.25,length=0.025))
  
  if(i%in%c(1,5,9,13)){axis_fun(2,ymaj,ymin,format(ymaj))}else{axis_fun(2,ymaj,ymin,NA)}
  if(i%in%c(13:16)){axis_fun(1,xmaj,xmin,month.abb[xmaj],line=-0.5)}else{axis_fun(1,xmaj,xmin,NA)}
  mtext(side=3,adj=0,line=-1.2,paste0(" ",wq.sites[i]),cex=0.7)
  box(lwd=1)
  
}
mtext(side=1,outer=T,line=1,"Month")
mtext(side=2,outer=T,line=1.5,"Arithmetic Mean TN (mg L\u207B\u00B9) ")

dev.off()

# RSMGL flows  ------------------------------------------------------------
dss.cat.fun=function(x,open.dss.fun=T){
  dss_out=if(open.dss.fun==T){opendss(x)}else{x}
  rslt=data.frame(path=getCatalogedPathnames(dss_out))
  str.val=strsplit(rslt$path,"/")
  rslt=data.frame(SITE=sapply(str.val,"[",3),TYPE=sapply(str.val,"[",4),
                  DateVal=sapply(str.val,"[",5))
  rslt=ddply(rslt,c("SITE","TYPE"),summarise,N.val=N.obs(SITE))
  return(rslt)
}


alts=c("FWO","FWOi","ECB22",paste0("ALT",c(21,22,23,24)))
n.alts=length(alts)
cols.alts=c("grey50","grey50","grey10",wesanderson::wes_palette("Zissou1",4,"continuous"))

test.FWO=dss.cat.fun(paste0(data.path,"Round2_RSMGL/",alts[1],"/RSMGL_output.dss"))
subset(test.FWO,SITE%in%struct.wq.list$RSM.site)
length(unique(struct.wq.list$RSM.site))
nrow(subset(test.FWO,SITE%in%struct.wq.list$RSM.site))
unique(struct.wq.list$RSM.site)[unique(struct.wq.list$RSM.site)%in%test.FWO$SITE==F]

test.Alt22=dss.cat.fun(paste0(data.path,"Round2_RSMGL/",alts[6],"/RSMGL_output.dss"))
test.Alt22.rslt=subset(test.Alt22,SITE%in%struct.wq.list$RSM.site)
length(unique(struct.wq.list$RSM.site))
nrow(subset(test.Alt22,SITE%in%struct.wq.list$RSM.site))
unique(struct.wq.list$RSM.site)[unique(struct.wq.list$RSM.site)%in%test.Alt22$SITE==F];# S197 is removed

RSM.sites=unique(struct.wq.list$RSM.site)
BBSEER.flow=data.frame()

for(j in 1:n.alts){
  dss_out=opendss(paste0(data.path,"Round2_RSMGL/",alts[j],"/RSMGL_output.dss"))  
  test=dss.cat.fun(dss_out,open.dss.fun=F)
  for(i in 1:length(RSM.sites)){
    
    if(nrow(subset(test,SITE==RSM.sites[i]))==1){
      paths=paste0("/RSMGL/",RSM.sites[i],"/FLOW//1DAY/SIMULATED/")  
      tmp=data.frame(getFullTSC(dss_out,paths))
      tmp$Date=date.fun(rownames(tmp))
      rownames(tmp)<-NULL
      tmp$SITE=RSM.sites[i]
      tmp$Alt=alts[j]
      BBSEER.flow=rbind(tmp,BBSEER.flow)
      print(i)
    }else{
      next
    }
  }
}

range(BBSEER.flow$FLOW)
BBSEER.flow$month=as.numeric(format(BBSEER.flow$Date,'%m'))
BBSEER.flow$CY=as.numeric(format(BBSEER.flow$Date,'%Y'))
BBSEER.flow$Q.acft=cfs.to.acftd(BBSEER.flow$FLOW)
BBSEER.flow$Q.Ld=BBSEER.flow$FLOW*2.447e6

struct.wq.list2=ddply(struct.wq.list,c("RSM.site","Station.ID","ENR","TRANSECT1"),summarise,
                      N.val=N.obs(RSM.site))

BBSEER.flow.mon.tot=ddply(BBSEER.flow,c("Alt","SITE","month","CY"),summarise,
                          mon.Q.Lmon=sum(Q.Ld,na.rm=T),
                          mon.Q.AcFtmon=sum(Q.acft,na.rm=T))


BBSEER.flow.mon.tot=merge(BBSEER.flow.mon.tot,struct.wq.list2,
                          by.x="SITE",by.y="RSM.site")
head(BBSEER.flow.mon.tot)
head(month.mean.val.xtab)

BBSEER.flow.mon.tot=merge(BBSEER.flow.mon.tot,month.mean.val.xtab,c("Station.ID","month"))
BBSEER.flow.mon.tot$TPload.kg=with(BBSEER.flow.mon.tot,TP*mon.Q.Lmon)*1e-6
BBSEER.flow.mon.tot$TNload.kg=with(BBSEER.flow.mon.tot,TN*mon.Q.Lmon)*1e-6

BBSEER.flow.mon.tot$SITE=with(BBSEER.flow.mon.tot,ifelse(SITE%in%paste0("S700_P",LETTERS[1:4]),"S700",SITE))

unique(BBSEER.flow.mon.tot$SITE)
### By site -----------------------------------------------------------------
ann.Q.load=ddply(BBSEER.flow.mon.tot,c("SITE","Alt","CY","ENR","TRANSECT1"),summarise,
      TFlow.Lyr=sum(mon.Q.Lmon,na.rm=T),
      TFlow.kacftyr=sum(mon.Q.AcFtmon,na.rm=T)/1000,
      TPload.kgyr=sum(TPload.kg,na.rm=T),
      TNload.kgyr=sum(TNload.kg,na.rm=T))
ann.Q.load=subset(ann.Q.load,SITE!="S21A")
ann.Q.load$TP.fwm=with(ann.Q.load,(TPload.kgyr*1e9)/TFlow.Lyr)
ann.Q.load$TN.fwm=with(ann.Q.load,(TNload.kgyr*1e6)/TFlow.Lyr)
ann.Q.load$Alt=factor(ann.Q.load$Alt,levels=alts)
ann.Q.load$SITE=factor(ann.Q.load$SITE,levels=RSM.site.order)
# write.csv(ann.Q.load,paste0(export.path,"BBSEER_ann_struct_flows.csv"),row.names = F)

unique(ann.Q.load$SITE)

mean.ann.Q.load=ddply(ann.Q.load,c("SITE","Alt","ENR","TRANSECT1"),summarise,
                      mean.Q.kacft=mean(TFlow.kacftyr,na.rm=T),
                      mean.TPLoad=mean(TPload.kgyr,na.rm=T),
                      mean.TNLoad=mean(TNload.kgyr,na.rm=T))
mean.ann.Q.load=merge(mean.ann.Q.load,
                      data.frame(expand.grid(SITE=unique(ann.Q.load$SITE),Alt=alts)),
                      c("SITE","Alt"),all.y=T)
unique(mean.ann.Q.load$SITE)

RSM.site.order
unique(ann.Q.load$SITE)
subset(mean.ann.Q.load,SITE=="S23A-S23B")
# png(filename=paste0(plot.path,"WQ_Rnd2_struct_Q.png"),width=8,height=5.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,1,0.5,0.5),oma=c(4,3,1,0.5),lwd=0.5);
layout(matrix(c(1:22,23,23),4,6,byrow=T))

ylim.val=c(0,300);by.y=100;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:22){
tmp=subset(mean.ann.Q.load,SITE==RSM.site.order[i])
x=barplot(tmp$mean.Q.kacft,col=adjustcolor(cols.alts,0.5),
        border=cols.alts,
        space=0.05,
        ylim=ylim.val,axes=F,ann=F)
text(x,tmp$mean.Q.kacft,format(round(tmp$mean.Q.kacft,1)),pos=3,cex=0.5,offset=0.25)
if(i%in%c(1,7,13,19)){axis_fun(2,ymaj,ymin,format(ymaj))}else{axis_fun(2,ymaj,ymin,NA)}
if(i%in%c(17:22)){axis_fun(1,x,x,alts,las=2)}else{axis_fun(1,x,x,NA)}
mtext(side=3,adj=0,line=-1.2,paste0(" ",RSM.site.order[i]),cex=0.7)
box(lwd=1)
}
mtext(side=1,outer=T,line=2.5,"Alternative")
mtext(side=2,outer=T,line=1.5,"Mean Discharge (x1000 Ac-Ft Yr\u207B\u00B9)")

plot(0:1,0:1,ann=F,axes=F,type="n")
legend(0.5,0.4,legend=c(alts),
       lty=c(0),lwd=c(0.1),col=c(cols.alts),pch=22,pt.bg=adjustcolor(cols.alts,0.5),
       ncol=3,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()


range(mean.ann.Q.load$mean.TPLoad,na.rm=T)
# png(filename=paste0(plot.path,"WQ_Rnd2_struct_TPLoad.png"),width=8,height=5.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,1,0.5,0.5),oma=c(4,3.5,1,0.5),lwd=0.5);
layout(matrix(c(1:22,23,23),4,6,byrow=T))

ylim.val=c(0,6000);by.y=1000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:22){
  tmp=subset(mean.ann.Q.load,SITE==RSM.site.order[i])
  x=barplot(tmp$mean.TPLoad,col=adjustcolor(cols.alts,0.5),
            border=cols.alts,
            space=0.05,
            ylim=ylim.val,axes=F,ann=F)
  text(x,tmp$mean.TPLoad,format(round(tmp$mean.TPLoad,0)),pos=3,cex=0.5,offset=0.25)
  if(i%in%c(1,7,13,19)){axis_fun(2,ymaj,ymin,format(ymaj))}else{axis_fun(2,ymaj,ymin,NA)}
  if(i%in%c(17:22)){axis_fun(1,x,x,alts,las=2)}else{axis_fun(1,x,x,NA)}
  mtext(side=3,adj=0,line=-1.2,paste0(" ",RSM.site.order[i]),cex=0.7)
  box(lwd=1)
}
mtext(side=1,outer=T,line=2.5,"Alternative")
mtext(side=2,outer=T,line=2,"Mean TP Load (kg Yr\u207B\u00B9)")

plot(0:1,0:1,ann=F,axes=F,type="n")
legend(0.5,0.4,legend=c(alts),
       lty=c(0),lwd=c(0.1),col=c(cols.alts),pch=22,pt.bg=adjustcolor(cols.alts,0.5),
       ncol=3,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

range(mean.ann.Q.load$mean.TNLoad,na.rm=T)
# png(filename=paste0(plot.path,"WQ_Rnd2_struct_TNLoad.png"),width=8,height=5.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,1,0.5,0.5),oma=c(4,3.5,1,0.5),lwd=0.5);
layout(matrix(c(1:22,23,23),4,6,byrow=T))

ylim.val=c(0,60e4);by.y=20e4;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:22){
  tmp=subset(mean.ann.Q.load,SITE==RSM.site.order[i])
  x=barplot(tmp$mean.TNLoad,col=adjustcolor(cols.alts,0.5),
            border=cols.alts,
            space=0.05,
            ylim=ylim.val,axes=F,ann=F)
  text(x,tmp$mean.TNLoad,format(round(tmp$mean.TNLoad/1e4,0)),pos=3,cex=0.5,offset=0.25)
  if(i%in%c(1,7,13,19)){axis_fun(2,ymaj,ymin,format(ymaj/1e4))}else{axis_fun(2,ymaj,ymin,NA)}
  if(i%in%c(17:22)){axis_fun(1,x,x,alts,las=2)}else{axis_fun(1,x,x,NA)}
  mtext(side=3,adj=0,line=-1.2,paste0(" ",RSM.site.order[i]),cex=0.7)
  box(lwd=1)
}
mtext(side=1,outer=T,line=2.5,"Alternative")
mtext(side=2,outer=T,line=2,"Mean TN Load (x10\u2074 kg Yr\u207B\u00B9)")

plot(0:1,0:1,ann=F,axes=F,type="n")
legend(0.5,0.4,legend=c(alts),
       lty=c(0),lwd=c(0.1),col=c(cols.alts),pch=22,pt.bg=adjustcolor(cols.alts,0.5),
       ncol=3,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()


## Percent difference
mean.TP.loads=ddply(ann.Q.load,c("SITE","Alt"),summarise,mean.TPLoad=mean(TPload.kgyr))
mean.TP.loads.PerDif.FWOi=merge(
  mean.TP.loads,
  dcast(mean.TP.loads,SITE~Alt,value.var="mean.TPLoad",mean)[,c("SITE","FWOi")],
  "SITE"
)
mean.TP.loads.PerDif.FWOi$PerDiff=with(mean.TP.loads.PerDif.FWOi,(mean.TPLoad-FWOi)/mean.TPLoad)*100

tmp=dcast(mean.TP.loads.PerDif.FWOi,SITE~Alt,value.var="PerDiff",mean)
tmp$y.val=nrow(tmp):1
tmp[is.na(tmp)]=0
tmp$ALT21[is.infinite(tmp$ALT21)]=0

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}

scale(tmp$ALT21)
scale_values(tmp$ALT21)

alt.vals=paste0("ALT2",1:4)

tmp2=tmp[,c("SITE",alt.vals)]
tmp2[,2:5]=sapply(tmp2[,2:5],scale_values)
tmp2$y.val=nrow(tmp2):1

# png(filename=paste0(plot.path,"WQ_Rnd2_TPPerDiff_site.png"),width=8,height=5.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.5,0.5),oma=c(3,5,1,0.5),lwd=0.5);
layout(matrix(1:4,1,4))
ylim.val=c(1,22);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/by.y)
xlim.val=c(-100,100);by.x=50;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
xlim.val2=c(-110,110)

for(i in 1:4){
plot(tmp$y.val~tmp[,alt.vals[i]],ylim=ylim.val,xlim=xlim.val2,axes=F,ann=F,type="n")
abline(v=xmaj,lty=2,col="grey",lwd=0.1)
abline(v=0)
segments(rep(0,nrow(tmp)),tmp$y.val,tmp[,alt.vals[i]],tmp$y.val,lwd=1.5)
points(tmp$y.val~tmp[,alt.vals[i]],pch=21,bg="firebrick1",cex=1.25)
text(ifelse(tmp[tmp[,alt.vals[i]]<0,alt.vals[i]]<(-100),-100,tmp[tmp[,alt.vals[i]]<0,alt.vals[i]]),
     tmp[tmp[,alt.vals[i]]<0,"y.val"],
     round(tmp[tmp[,alt.vals[i]]<0,alt.vals[i]],1),
     pos=ifelse(tmp[tmp[,alt.vals[i]]<0,alt.vals[i]]<(-100),4,2),offset=0.5)
text(tmp[tmp[,alt.vals[i]]>=0,alt.vals[i]],
     tmp[tmp[,alt.vals[i]]>=0,"y.val"],
     round(tmp[tmp[,alt.vals[i]]>=0,alt.vals[i]],1),
     pos=ifelse(tmp[tmp[,alt.vals[i]]>=0,alt.vals[i]]>80,2,4),offset=0.5)
if(i==1){axis_fun(2,ymaj,ymaj,rev(RSM.site.order))}else{axis_fun(2,ymaj,ymaj,NA)}
axis_fun(1,xmaj,xmin,xmaj);box(lwd=1)
mtext(side=3,adj=0,alt.vals[i])
}
mtext(side=1,outer=T,line=1,"TP Load Avg Percent Difference to FWOi")
dev.off()

par(family="serif",mar=c(1,2,0.5,0.5),oma=c(3,5,1,0.5),lwd=0.5);
layout(matrix(1:4,1,4))
ylim.val=c(1,22);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/by.y)
xlim.val=c(0,1);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
for(i in 1:4){
  plot(tmp2$y.val~tmp2[,alt.vals[i]],ylim=ylim.val,xlim=xlim.val,axes=F,ann=F,type="n")
  abline(v=xmaj,lty=2,col="grey",lwd=0.1)
  
  segments(rep(0,nrow(tmp2)),tmp2$y.val,tmp2[,alt.vals[i]],tmp2$y.val,lwd=1.5)
  points(tmp2$y.val~tmp2[,alt.vals[i]],pch=21,bg="firebrick1",cex=1.25)
  if(i==1){axis_fun(2,ymaj,ymaj,rev(RSM.site.order))}else{axis_fun(2,ymaj,ymaj,NA)}
  axis_fun(1,xmaj,xmin,format(xmaj));box(lwd=1)
  mtext(side=3,adj=0,alt.vals[i])
}
mtext(side=1,outer=T,line=1,"Rescaled TP Load Avg Percent Difference to FWOi")
dev.off()

## TN
mean.TN.loads=ddply(ann.Q.load,c("SITE","Alt"),summarise,mean.TNLoad=mean(TNload.kgyr))
mean.TN.loads.PerDif.FWOi=merge(
  mean.TN.loads,
  dcast(mean.TN.loads,SITE~Alt,value.var="mean.TNLoad",mean)[,c("SITE","FWOi")],
  "SITE"
)
mean.TN.loads.PerDif.FWOi$PerDiff=with(mean.TN.loads.PerDif.FWOi,(mean.TNLoad-FWOi)/mean.TNLoad)*100

tmp=dcast(mean.TN.loads.PerDif.FWOi,SITE~Alt,value.var="PerDiff",mean)
tmp$y.val=nrow(tmp):1
tmp[is.na(tmp)]=0
tmp$ALT21[is.infinite(tmp$ALT21)]=0

# png(filename=paste0(plot.path,"WQ_Rnd2_TNPerDiff_site.png"),width=8,height=5.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.5,0.5),oma=c(3,5,1,0.5),lwd=0.5);
layout(matrix(1:4,1,4))
ylim.val=c(1,22);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/by.y)
xlim.val=c(-100,100);by.x=50;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
xlim.val2=c(-110,110)

for(i in 1:4){
  plot(tmp$y.val~tmp[,alt.vals[i]],ylim=ylim.val,xlim=xlim.val2,axes=F,ann=F,type="n")
  abline(v=xmaj,lty=2,col="grey",lwd=0.1)
  abline(v=0)
  segments(rep(0,nrow(tmp)),tmp$y.val,tmp[,alt.vals[i]],tmp$y.val,lwd=1.5)
  points(tmp$y.val~tmp[,alt.vals[i]],pch=21,bg="firebrick1",cex=1.25)
  text(ifelse(tmp[tmp[,alt.vals[i]]<0,alt.vals[i]]<(-100),-100,tmp[tmp[,alt.vals[i]]<0,alt.vals[i]]),
       tmp[tmp[,alt.vals[i]]<0,"y.val"],
       round(tmp[tmp[,alt.vals[i]]<0,alt.vals[i]],1),
       pos=ifelse(tmp[tmp[,alt.vals[i]]<0,alt.vals[i]]<(-100),4,2),offset=0.5)
  text(tmp[tmp[,alt.vals[i]]>=0,alt.vals[i]],
       tmp[tmp[,alt.vals[i]]>=0,"y.val"],
       round(tmp[tmp[,alt.vals[i]]>=0,alt.vals[i]],1),
       pos=ifelse(tmp[tmp[,alt.vals[i]]>=0,alt.vals[i]]>80,2,4),offset=0.5)
  if(i==1){axis_fun(2,ymaj,ymaj,rev(RSM.site.order))}else{axis_fun(2,ymaj,ymaj,NA)}
  axis_fun(1,xmaj,xmin,xmaj);box(lwd=1)
  mtext(side=3,adj=0,alt.vals[i])
}
mtext(side=1,outer=T,line=1,"TN Load Avg Percent Difference to FWOi")
dev.off()



# ENR  --------------------------------------------------------------------
ann.ENR.Q.load=ddply(BBSEER.flow.mon.tot,c("ENR","Alt","CY"),summarise,
                 TFlow.Lyr=sum(mon.Q.Lmon,na.rm=T),
                 TFlow.kacftyr=sum(mon.Q.AcFtmon,na.rm=T)/1000,
                 TPload.kgyr=sum(TPload.kg,na.rm=T),
                 TNload.kgyr=sum(TNload.kg,na.rm=T))
ann.ENR.Q.load$ENR=factor(ann.ENR.Q.load$ENR,levels=paste0("ENRH",c(5,9,3,6,2)))
unique(ann.ENR.Q.load$ENR)
plot(TPload.kgyr~CY,subset(ann.ENR.Q.load,ENR=="ENRH2"&Alt=="FWO"))

ENR.vals=paste0("ENRH",c(5,9,3,6,2))
alts

xlim.val=c(1965,2016);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0,12000);by.y=5000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)

# png(filename=paste0(plot.path,"WQ_Rnd2_TP_ENR_TS.png"),width=8,height=5.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.5,1,0.5,0.5),oma=c(3,4,1,0.5),lwd=0.5);
layout(matrix(1:35,5,7))

for(j in 1:length(alts)){
  for(i in 1:length(ENR.vals)){
    plot(TPload.kgyr~CY,ann.ENR.Q.load,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
    abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
    lines(TPload.kgyr~CY,subset(ann.ENR.Q.load,ENR==ENR.vals[i]&Alt==alts[j]),col="red",lwd=1.2)
    if(i==length(ENR.vals)){axis_fun(1,xmaj,xmin,xmaj,line=-0.5)}else{axis_fun(1,xmaj,xmin,NA)}
    if(j==1){axis_fun(2,ymaj,ymin,ymaj)}else{axis_fun(2,ymaj,ymin,NA)}
    if(j==1){mtext(side=3,adj=0,line=-1.25,paste0(" ",ENR.vals[i]),cex=0.75)}
    if(i==1){mtext(side=3,alts[j])}
    box(lwd=1)
  }
}
mtext(side=2,outer=T,line=2.5,"TP Load (kg yr\u207B\u00B9)")
mtext(side=1,outer=T,line=1.5,"Calendar Year")
dev.off()


mean.ann.enr.Q.load=ddply(ann.ENR.Q.load,c("ENR","Alt"),summarise,
                            mean.Q=mean(TFlow.kacftyr,na.rm=T),
                            mean.TPLoad=mean(TPload.kgyr,na.rm=T),
                            mean.TNLoad=mean(TNload.kgyr,na.rm=T))
mean.ann.enr.Q.load$Alt=factor(mean.ann.enr.Q.load$Alt,levels=alts)
mean.ann.enr.Q.load$ENR=factor(mean.ann.enr.Q.load$ENR,levels=paste0("ENRH",c(5,9,3,6,2)))
mean.ann.enr.Q.load=mean.ann.enr.Q.load[order(mean.ann.enr.Q.load$Alt),]


ENR.vals=paste0("ENRH",c(5,9,3,6))
# png(filename=paste0(plot.path,"WQ_Rnd2_ENR_QLoads.png"),width=6.5,height=5.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.5,1),oma=c(4,1,1,0.5),lwd=0.5);
layout(matrix(1:12,4,3,byrow=F))

ylim.val=c(0,800);by.y=200;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:length(ENR.vals)){
  tmp=subset(mean.ann.enr.Q.load,ENR==ENR.vals[i])
  x=barplot(tmp$mean.Q,col=adjustcolor(cols.alts,0.5),
            border=cols.alts,
            space=0.05,
            ylim=ylim.val,axes=F,ann=F)
  text(x,tmp$mean.Q,format(round(tmp$mean.Q,1)),pos=3,cex=0.75,offset=0.25)
  axis_fun(2,ymaj,ymin,format(ymaj))
  if(i==4){axis_fun(1,x,x,alts,las=2)}else{axis_fun(1,x,x,NA)}
  mtext(side=3,adj=0,line=-1.2,paste0(" ",ENR.vals[i]),cex=0.7)
  box(lwd=1)
  if(i==1){mtext(side=3,adj=0,"Structure Only Discharge & Load")}
}
mtext(side=2,outer=T,line=-0.5,"Mean Discharge (x10\u00B3 Ac-Ft Yr\u207B\u00B9)")

ylim.val=c(0,15000);by.y=5000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:length(ENR.vals)){
  tmp=subset(mean.ann.enr.Q.load,ENR==ENR.vals[i])
  x=barplot(tmp$mean.TPLoad,col=adjustcolor(cols.alts,0.5),
            border=cols.alts,
            space=0.05,
            ylim=ylim.val,axes=F,ann=F)
  text(x,tmp$mean.TPLoad,format(round(tmp$mean.TPLoad/1e3,1)),pos=3,cex=0.75,offset=0.25)
  axis_fun(2,ymaj,ymin,format(ymaj/1e3))
  if(i==4){axis_fun(1,x,x,alts,las=2)}else{axis_fun(1,x,x,NA)}
  mtext(side=3,adj=0,line=-1.2,paste0(" ",ENR.vals[i]),cex=0.7)
  box(lwd=1)
}
mtext(side=2,outer=T,line=-17,"Mean TP Load (x10\u00B3 kg Yr\u207B\u00B9)")

ylim.val=c(0,150e4);by.y=50e4;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:length(ENR.vals)){
  tmp=subset(mean.ann.enr.Q.load,ENR==ENR.vals[i])
  x=barplot(tmp$mean.TNLoad,col=adjustcolor(cols.alts,0.5),
            border=cols.alts,
            space=0.05,
            ylim=ylim.val,axes=F,ann=F)
  text(x,tmp$mean.TNLoad,format(round(tmp$mean.TNLoad/1e4,1)),pos=3,cex=0.75,offset=0.25)
  axis_fun(2,ymaj,ymin,format(ymaj/1e4))
  if(i==4){axis_fun(1,x,x,alts,las=2)}else{axis_fun(1,x,x,NA)}
  mtext(side=3,adj=0,line=-1.2,paste0(" ",ENR.vals[i]),cex=0.7)
  box(lwd=1)
  
}
mtext(side=2,outer=T,line=-33,"Mean TN Load (x10\u2074 kg Yr\u207B\u00B9)")
mtext(side=1,outer=T,line=2.5,"Alternative")
dev.off()

### 
mean.TP.loads=ddply(subset(ann.ENR.Q.load,ENR%in%ENR.vals),c("ENR","Alt"),summarise,
                    mean.Q=mean(TFlow.kacftyr),
                    mean.TPLoad=mean(TPload.kgyr),
                    mean.TNLoad=mean(TNload.kgyr))
mean.TP.loads.PerDif.FWOi=merge(
  mean.TP.loads,
  dcast(mean.TP.loads,ENR~Alt,value.var="mean.TPLoad",mean)[,c("ENR","FWOi")],
  "ENR"
)
mean.TP.loads.PerDif.FWOi$PerDiff=with(mean.TP.loads.PerDif.FWOi,(mean.TPLoad-FWOi)/mean.TPLoad)*100
mean.TP.loads.PerDif.FWOi.sum=dcast(mean.TP.loads.PerDif.FWOi,ENR~Alt,value.var="PerDiff",mean)
mean.TP.loads.PerDif.FWOi.sum$y.val=nrow(mean.TP.loads.PerDif.FWOi.sum):1

mean.TN.loads=ddply(subset(ann.ENR.Q.load,ENR%in%ENR.vals),c("ENR","Alt"),summarise,
                    mean.Q=mean(TFlow.kacftyr),
                    mean.TNLoad=mean(TNload.kgyr),
                    mean.TNLoad=mean(TNload.kgyr))
mean.TN.loads.PerDif.FWOi=merge(
  mean.TN.loads,
  dcast(mean.TN.loads,ENR~Alt,value.var="mean.TNLoad",mean)[,c("ENR","FWOi")],
  "ENR"
)
mean.TN.loads.PerDif.FWOi$PerDiff=with(mean.TN.loads.PerDif.FWOi,(mean.TNLoad-FWOi)/mean.TNLoad)*100
mean.TN.loads.PerDif.FWOi.sum=dcast(mean.TN.loads.PerDif.FWOi,ENR~Alt,value.var="PerDiff",mean)
mean.TN.loads.PerDif.FWOi.sum$y.val=nrow(mean.TN.loads.PerDif.FWOi.sum):1


# png(filename=paste0(plot.path,"WQ_Rnd2_LoadPerDiff_ENR.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(3,1,0.5,0.5),oma=c(0,4.5,1,0.5),lwd=0.5);
layout(matrix(1:8,2,4,byrow=T))
ylim.val=c(1,4);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/by.y)
xlim.val=c(-100,100);by.x=50;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)

for(i in 1:4){
  plot(mean.TP.loads.PerDif.FWOi.sum$y.val~mean.TP.loads.PerDif.FWOi.sum[,alt.vals[i]],ylim=ylim.val,xlim=xlim.val,axes=F,ann=F,type="n")
  abline(v=xmaj,lty=2,col="grey",lwd=0.1)
  abline(v=0)
  segments(rep(0,nrow(mean.TP.loads.PerDif.FWOi.sum)),
           mean.TP.loads.PerDif.FWOi.sum$y.val,
           mean.TP.loads.PerDif.FWOi.sum[,alt.vals[i]],
           mean.TP.loads.PerDif.FWOi.sum$y.val,lwd=1.5)
  points(mean.TP.loads.PerDif.FWOi.sum$y.val~mean.TP.loads.PerDif.FWOi.sum[,alt.vals[i]],pch=21,bg="firebrick1",cex=1.25)
  if(i==1){axis_fun(2,ymaj,ymaj,rev(ENR.vals))}else{axis_fun(2,ymaj,ymaj,NA)}
  axis_fun(1,xmaj,xmin,xmaj,line=-0.5);box(lwd=1)
  mtext(side=3,adj=0,alt.vals[i])
}
mtext(side=1,outer=T,line=-18,"TP Load Avg Percent Difference to FWOi")

for(i in 1:4){
  plot(mean.TN.loads.PerDif.FWOi.sum$y.val~mean.TN.loads.PerDif.FWOi.sum[,alt.vals[i]],ylim=ylim.val,xlim=xlim.val,axes=F,ann=F,type="n")
  abline(v=xmaj,lty=2,col="grey",lwd=0.1)
  abline(v=0)
  segments(rep(0,nrow(mean.TN.loads.PerDif.FWOi.sum)),
           mean.TN.loads.PerDif.FWOi.sum$y.val,
           mean.TN.loads.PerDif.FWOi.sum[,alt.vals[i]],
           mean.TN.loads.PerDif.FWOi.sum$y.val,lwd=1.5)
  points(mean.TN.loads.PerDif.FWOi.sum$y.val~mean.TN.loads.PerDif.FWOi.sum[,alt.vals[i]],pch=21,bg="firebrick1",cex=1.25)
  if(i==1){axis_fun(2,ymaj,ymaj,rev(ENR.vals))}else{axis_fun(2,ymaj,ymaj,NA)}
  axis_fun(1,xmaj,xmin,xmaj,line=-0.5);box(lwd=1)
}
mtext(side=1,outer=T,line=-1,"TN Load Avg Percent Difference to FWOi")
mtext(side=2,outer=T,line=3,"Estuary NNC Segment")
dev.off()




# transect  --------------------------------------------------------------------
ann.trans.Q.load=ddply(BBSEER.flow.mon.tot,c("TRANSECT1","Alt","CY"),summarise,
                     TFlow.Lyr=sum(mon.Q.Lmon,na.rm=T),
                     TFlow.kacftyr=sum(mon.Q.AcFtmon,na.rm=T)/1000,
                     TPload.kgyr=sum(TPload.kg,na.rm=T),
                     TNload.kgyr=sum(TNload.kg,na.rm=T))

ann.trans.Q.load$TRANSECT1=factor(ann.trans.Q.load$TRANSECT1,levels=paste0("ST0",c(1,2,3,5)))
unique(ann.trans.Q.load$TRANSECT1)
plot(TPload.kgyr~CY,subset(ann.trans.Q.load,TRANSECT1=="ST01"&Alt=="FWO"))

trans.vals=paste0("ST0",c(1,2,3,5))
alts

xlim.val=c(1965,2016);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0,20000);by.y=5000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)

# png(filename=paste0(plot.path,"WQ_Rnd2_TP_trans1_TS.png"),width=8,height=5.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.5,1,0.5,0.5),oma=c(3,4,1,0.5),lwd=0.5);
layout(matrix(1:28,4,7))

for(j in 1:length(alts)){
  for(i in 1:length(trans.vals)){
    plot(TPload.kgyr~CY,ann.trans.Q.load,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
    abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
    lines(TPload.kgyr~CY,subset(ann.trans.Q.load,TRANSECT1==trans.vals[i]&Alt==alts[j]),col="red",lwd=1.2)
    if(i==length(trans.vals)){axis_fun(1,xmaj,xmin,xmaj,line=-0.5)}else{axis_fun(1,xmaj,xmin,NA)}
    if(j==1){axis_fun(2,ymaj,ymin,ymaj)}else{axis_fun(2,ymaj,ymin,NA)}
    if(j==1){mtext(side=3,adj=0,line=-1.25,paste0(" ",trans.vals[i]),cex=0.75)}
    if(i==1){mtext(side=3,alts[j])}
    box(lwd=1)
  }
}
mtext(side=2,outer=T,line=2.5,"TP Load (kg yr\u207B\u00B9)")
mtext(side=1,outer=T,line=1.5,"Calendar Year")
dev.off()


mean.ann.trans.Q.load=ddply(subset(ann.trans.Q.load,TRANSECT1%in%c("ST01","ST02","ST03")),c("TRANSECT1","Alt"),summarise,
                    mean.Q=mean(TFlow.kacftyr,na.rm=T),
                    mean.TPLoad=mean(TPload.kgyr,na.rm=T),
                    mean.TNLoad=mean(TNload.kgyr,na.rm=T))
mean.ann.trans.Q.load$Alt=factor(mean.ann.trans.Q.load$Alt,levels=alts)
mean.ann.trans.Q.load=mean.ann.trans.Q.load[order(mean.ann.trans.Q.load$TRANSECT1,mean.ann.trans.Q.load$Alt),]

trans.vals=paste0("ST0",c(1,2,3))
# png(filename=paste0(plot.path,"WQ_Rnd2_trans_QLoads.png"),width=6.5,height=5.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.5,1),oma=c(4,1,1,0.5),lwd=0.5);
layout(matrix(1:9,3,3,byrow=F))

ylim.val=c(0,800);by.y=200;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:length(trans.vals)){
  tmp=subset(mean.ann.trans.Q.load,TRANSECT1==trans.vals[i])
  x=barplot(tmp$mean.Q,col=adjustcolor(cols.alts,0.5),
            border=cols.alts,
            space=0.05,
            ylim=ylim.val,axes=F,ann=F)
  text(x,tmp$mean.Q,format(round(tmp$mean.Q,1)),pos=3,cex=0.75,offset=0.25)
  axis_fun(2,ymaj,ymin,format(ymaj))
  if(i==3){axis_fun(1,x,x,alts,las=2)}else{axis_fun(1,x,x,NA)}
  mtext(side=3,adj=0,line=-1.2,paste0(" ",trans.vals[i]),cex=0.7)
  box(lwd=1)
  if(i==1){mtext(side=3,adj=0,"Structure Only Discharge & Load")}
  if(i==2){mtext(side=2,line=2.5,"Mean Discharge (x10\u00B3 Ac-Ft Yr\u207B\u00B9)")}
}

ylim.val=c(0,15000);by.y=5000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:length(trans.vals)){
  tmp=subset(mean.ann.trans.Q.load,TRANSECT1==trans.vals[i])
  x=barplot(tmp$mean.TPLoad,col=adjustcolor(cols.alts,0.5),
            border=cols.alts,
            space=0.05,
            ylim=ylim.val,axes=F,ann=F)
  text(x,tmp$mean.TPLoad,format(round(tmp$mean.TPLoad/1e3,1)),pos=3,cex=0.75,offset=0.25)
  axis_fun(2,ymaj,ymin,format(ymaj/1e3))
  if(i==3){axis_fun(1,x,x,alts,las=2)}else{axis_fun(1,x,x,NA)}
  mtext(side=3,adj=0,line=-1.2,paste0(" ",trans.vals[i]),cex=0.7)
  box(lwd=1)
  if(i==2){mtext(side=2,line=2.5,"Mean TP Load (x10\u00B3 kg Yr\u207B\u00B9)")}
}

ylim.val=c(0,150e4);by.y=50e4;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:length(trans.vals)){
  tmp=subset(mean.ann.trans.Q.load,TRANSECT1==trans.vals[i])
  x=barplot(tmp$mean.TNLoad,col=adjustcolor(cols.alts,0.5),
            border=cols.alts,
            space=0.05,
            ylim=ylim.val,axes=F,ann=F)
  text(x,tmp$mean.TNLoad,format(round(tmp$mean.TNLoad/1e4,1)),pos=3,cex=0.75,offset=0.25)
  axis_fun(2,ymaj,ymin,format(ymaj/1e4))
  if(i==3){axis_fun(1,x,x,alts,las=2)}else{axis_fun(1,x,x,NA)}
  mtext(side=3,adj=0,line=-1.2,paste0(" ",trans.vals[i]),cex=0.7)
  box(lwd=1)
  if(i==2){mtext(side=2,line=2.5,"Mean TN Load (x10\u2074 kg Yr\u207B\u00B9)")}
}
mtext(side=1,outer=T,line=2.5,"Alternative")
dev.off()

### 
mean.TP.loads=ddply(subset(ann.trans.Q.load,TRANSECT1%in%c("ST01","ST02","ST03")),c("TRANSECT1","Alt"),summarise,
                    mean.Q=mean(TFlow.kacftyr),
                    mean.TPLoad=mean(TPload.kgyr),
                    mean.TNLoad=mean(TNload.kgyr))
mean.TP.loads.PerDif.FWOi=merge(
  mean.TP.loads,
  dcast(mean.TP.loads,TRANSECT1~Alt,value.var="mean.TPLoad",mean)[,c("TRANSECT1","FWOi")],
  "TRANSECT1"
)
mean.TP.loads.PerDif.FWOi$PerDiff=with(mean.TP.loads.PerDif.FWOi,(mean.TPLoad-FWOi)/mean.TPLoad)*100

mean.TP.loads.PerDif.FWOi.sum=dcast(mean.TP.loads.PerDif.FWOi,TRANSECT1~Alt,value.var="PerDiff",mean)
mean.TP.loads.PerDif.FWOi.sum$y.val=nrow(tmp):1

mean.TN.loads.PerDif.FWOi=merge(
  mean.TP.loads,
  dcast(mean.TP.loads,TRANSECT1~Alt,value.var="mean.TNLoad",mean)[,c("TRANSECT1","FWOi")],
  "TRANSECT1"
)
mean.TN.loads.PerDif.FWOi$PerDiff=with(mean.TN.loads.PerDif.FWOi,(mean.TNLoad-FWOi)/mean.TNLoad)*100

mean.TN.loads.PerDif.FWOi.sum=dcast(mean.TN.loads.PerDif.FWOi,TRANSECT1~Alt,value.var="PerDiff",mean)
mean.TN.loads.PerDif.FWOi.sum$y.val=nrow(tmp):1


# png(filename=paste0(plot.path,"WQ_Rnd2_LoadPerDiff_trans.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(3,1,0.5,0.5),oma=c(0,4,1,0.5),lwd=0.5);
layout(matrix(1:8,2,4,byrow=T))
ylim.val=c(1,3);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/by.y)
xlim.val=c(-75,50);by.x=25;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)


for(i in 1:4){
  plot(mean.TP.loads.PerDif.FWOi.sum$y.val~mean.TP.loads.PerDif.FWOi.sum[,alt.vals[i]],ylim=ylim.val,xlim=xlim.val,axes=F,ann=F,type="n")
  abline(v=xmaj,lty=2,col="grey",lwd=0.1)
  abline(v=0)
  segments(rep(0,nrow(mean.TP.loads.PerDif.FWOi.sum)),
           mean.TP.loads.PerDif.FWOi.sum$y.val,
           mean.TP.loads.PerDif.FWOi.sum[,alt.vals[i]],
           mean.TP.loads.PerDif.FWOi.sum$y.val,lwd=1.5)
  points(mean.TP.loads.PerDif.FWOi.sum$y.val~mean.TP.loads.PerDif.FWOi.sum[,alt.vals[i]],pch=21,bg="firebrick1",cex=1.25)
  if(i==1){axis_fun(2,ymaj,ymaj,rev(trans.vals))}else{axis_fun(2,ymaj,ymaj,NA)}
  axis_fun(1,xmaj,xmin,xmaj,line=-0.5);box(lwd=1)
  mtext(side=3,adj=0,alt.vals[i])
}
mtext(side=1,outer=T,line=-18,"TP Load Avg Percent Difference to FWOi")

for(i in 1:4){
  plot(mean.TN.loads.PerDif.FWOi.sum$y.val~mean.TN.loads.PerDif.FWOi.sum[,alt.vals[i]],ylim=ylim.val,xlim=xlim.val,axes=F,ann=F,type="n")
  abline(v=xmaj,lty=2,col="grey",lwd=0.1)
  abline(v=0)
  segments(rep(0,nrow(mean.TN.loads.PerDif.FWOi.sum)),
           mean.TN.loads.PerDif.FWOi.sum$y.val,
           mean.TN.loads.PerDif.FWOi.sum[,alt.vals[i]],
           mean.TN.loads.PerDif.FWOi.sum$y.val,lwd=1.5)
  points(mean.TN.loads.PerDif.FWOi.sum$y.val~mean.TN.loads.PerDif.FWOi.sum[,alt.vals[i]],pch=21,bg="firebrick1",cex=1.25)
  if(i==1){axis_fun(2,ymaj,ymaj,rev(trans.vals))}else{axis_fun(2,ymaj,ymaj,NA)}
  axis_fun(1,xmaj,xmin,xmaj,line=-0.5);box(lwd=1)
}
mtext(side=1,outer=T,line=-1,"TN Load Avg Percent Difference to FWOi")
mtext(side=2,outer=T,line=2,"Transect")
dev.off()
