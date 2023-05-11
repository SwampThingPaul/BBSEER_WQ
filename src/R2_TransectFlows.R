## BBSEER - SRS Inflow evalution
##
## Code was compiled by Paul Julian
## contact info: pjulian@evergladesfoundation.org

## BAD 
## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
## Clears Everything...start fresh.
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

# GIS Libraries
library(rgdal)
library(rgeos)
library(raster)
library(gstat)
library(tmap)

## Paths
wd="C:/Julian_LaCie/_GitHub/BBSEER_WQ"

paths=paste0(wd,c("/Plots/","/Exports/","/Data/","/src/","/GIS"))
# Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]
gis.path=paths[5]
GIS.path.gen="C:/Julian_LaCie/_GISData"
# Helper variables
nad83.pro=CRS("+init=epsg:4269")
utm17=CRS("+init=epsg:26917")
wgs84=CRS("+init=epsg:4326")
harnNAD83=CRS("+init=epsg:2881")

tmap_mode("view")


# -------------------------------------------------------------------------
alts=c("FWO","FWOi","ECB22",paste0("ALT",c(21,22,23,24)))
n.alts=length(alts)

cols.alts=c("grey50","grey50","grey10",wesanderson::wes_palette("Zissou1",4,"continuous"))


# transect locs -----------------------------------------------------------
library(rvest)
rslt=rvest::read_html(paste0(data.path,"flow_transects_ECB22.xml"))
rslt=rvest::html_nodes(rslt,"flowgage")
# html_text(rslt[1])
# html_attr(rslt,"section")
# html_attr(rslt,"label")

trans.loc=data.frame()
for(i in 1:length(rslt)){
  nodes.vals=as.numeric(strsplit(gsub("[{}]","",html_text(rslt[i])), split = "\\s+")[[1]])
  nodes.vals=nodes.vals[is.na(nodes.vals)==F]
  sec.val=html_attr(rslt[i],"section")
  lab.val=html_attr(rslt[i],"label")
  
  dat.split=data.frame(transect=lab.val,
           section=sec.val,
           grid.nodes=nodes.vals)
  trans.loc=rbind(trans.loc,dat.split)
}

mesh=spTransform(readOGR(paste0(gis.path,"/Round2"),"RSMGL_mesh"),utm17)
# plot(mesh)

meshnodes=spTransform(readOGR(paste0(gis.path,"/Round2"),"meshnodes"),utm17)
# plot(meshnodes,pch=21)
meshnodes=cbind(meshnodes,coordinates(meshnodes))
colnames(meshnodes@data)=c("MESH2D", "nodeid", "x_lat", "y_long", "z", "UTMX", "UTMY")

meshnodes2=meshnodes@data
meshnodes2$nodeid=as.numeric(meshnodes2$nodeid)
meshnodes2=merge(meshnodes2,subset(trans.loc,section=="ol"),by.x="nodeid",by.y="grid.nodes",all.y=T)

subset(meshnodes2,transect=="T27_ol")
meshnodes2=meshnodes2[order(meshnodes2$nodeid,meshnodes2$transect),]
meshnodes2=meshnodes2[order(meshnodes2$transect,meshnodes2$UTMY,meshnodes2$UTMX),]

meshnodes2$dUTMY=with(meshnodes2,ave(UTMY,transect,FUN=function(x)c(0,diff(x))))
meshnodes2$dUTMX=with(meshnodes2,ave(UTMX,transect,FUN=function(x)c(0,diff(x))))
meshnodes2$dist.m=with(meshnodes2,sqrt((dUTMY^2)+(dUTMX^2)));#calculates distance between points
## subset(meshnodes2,transect=="T18N_ol")
meshnodes2=SpatialPointsDataFrame(coords=meshnodes2[,c("UTMX","UTMY")],data=meshnodes2,proj4string = utm17)
# writeOGR(meshnodes2,paste0(export.path,"GIS"),"TransPts",driver="ESRI Shapefile")
path=sp::split(meshnodes2,meshnodes2$transect)

trans.id=ddply(meshnodes2@data,c("transect"),summarise,N.val=N.obs(transect),max.dist=max(dist.m,na.rm=T),tot.dist=sum(dist.m,na.rm=T))

##Convert spatial points to spatial lines for each hurricane
sp_lines=SpatialLinesDataFrame(SpatialLines(list(Lines(list(Line(path[[1]])),unique(path[[1]]@data$transect))),utm17),data.frame(row.names=trans.id$transect,transect=trans.id$transect))
pb=txtProgressBar(1,max=length(path),style=3)
for(i in 2:length(path)){
  tmp=SpatialLinesDataFrame(SpatialLines(list(Lines(list(Line(path[[i]])),unique(path[[i]]@data$transect))),utm17),data.frame(row.names=trans.id$transect,transect=trans.id$transect))
  #tmp=SpatialLines(list(Lines(list(Line(path[[i]])),unique(path[[i]]@data$Key))),CRS("+init=epsg:4269"))
  sp_lines=rbind(sp_lines,tmp)
  setTxtProgressBar(pb,i)
}
chk=data.frame(gIsValid(sp_lines,byid=T));#checks
chk$Key=rownames(chk)
colnames(chk)=c("geo.valid","Key")
subset(chk,geo.valid=="FALSE")
# writeOGR(sp_lines,paste0(export.path,"GIS"),"TransLines",driver="ESRI Shapefile",overwrite_layer = T)

plot(subset(sp_lines,transect=="T27_ol"))
plot(subset(sp_lines,transect=="ST01_ol"))
plot(subset(sp_lines,transect=="ST03_ol"))

# ogrListLayers(paste0(GIS.path.gen,"/AHED_release/20230405/AHED.gdb"))
# ogrListLayers(paste0(GIS.path.gen,"/AHED_release/AHED_20171102.gdb"))
canals=spTransform(readOGR(paste0(GIS.path.gen,"/SFER_GIS_Geodatabase.gdb"),"SFWMD_Canals"),utm17)
ENP=spTransform(readOGR(paste0(GIS.path.gen,"/SFER_GIS_Geodatabase.gdb"),"ENP"),utm17)
shore=spTransform(readOGR(paste0(GIS.path.gen,"/FWC"),"FWC_Shoreline"),utm17)

trans.mean.center=ddply(meshnodes2@data,c("transect"),summarise,mean.UTMX=mean(UTMX),mean.UTMY=mean(UTMY))
trans.mean.center=SpatialPointsDataFrame(coords=trans.mean.center[,c("mean.UTMX","mean.UTMY")],data=trans.mean.center,proj4string = utm17)
trans.mean.center$trans.lab=gsub("_ol","",trans.mean.center$transect)

# png(filename=paste0(plot.path,"Round2/T23_QuickMap.png"),width=5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",oma=c(1,1,1,1),mar=c(0.1,0.1,0.1,0.1),xpd=F)

bbox.lims=bbox(gBuffer(subset(sp_lines,grepl("T23",transect)==T),width=1000))

plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(mesh,add=T,border="grey")
plot(ENP,col=NA,add=T)
plot(canals,add=T,col="dodgerblue2",lwd=1.5)
plot(subset(sp_lines,grepl("T23",transect)==T),col=c("red","blue","lightgreen"),add=T,lwd=2)
raster::text(subset(trans.mean.center,grepl("T23",transect)==T),"trans.lab",pos=3)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=0.75,seg.len=4);box(lwd=1)
dev.off()


subset(sp_lines,grepl("T19",transect)==T)
bbox.lims=bbox(gBuffer(subset(sp_lines,grepl("T19",transect)==T),width=1000))
plot(shore,col="cornsilk",border="grey",bg="lightblue",ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(mesh,add=T,border="grey")
plot(ENP,col=NA,add=T)
plot(canals,add=T,col="dodgerblue2",lwd=1.5)
plot(subset(sp_lines,grepl("T19",transect)==T),col=c("red","blue","lightgreen"),add=T,lwd=2)

plot(subset(sp_lines,transect=="ST01_ol"))

# -------------------------------------------------------------------------
dss.cat.fun=function(x,open.dss.fun=T){
  dss_out=if(open.dss.fun==T){opendss(x)}else{x}
  rslt=data.frame(path=getCatalogedPathnames(dss_out))
  str.val=strsplit(rslt$path,"/")
  rslt=data.frame(SITE=sapply(str.val,"[",3),TYPE=sapply(str.val,"[",4),
                  DateVal=sapply(str.val,"[",5))
  rslt=ddply(rslt,c("SITE","TYPE"),summarise,N.val=N.obs(SITE))
  return(rslt)
}



dss_out=opendss(paste0(data.path,"Round2_RSMGL/",alts[4],"/transect_flows.dss"))
dss_cat=dss.cat.fun(dss_out,F)

unique(dss_cat$SITE)

subset(dss_cat,LOCATION=="T23A_TRANSECT")

RSM.sites=c(paste0("T23",LETTERS[1:3],"_TRANSECT"))

trans_flow.ol=data.frame()
for(j in 1:n.alts){
  dss_out=opendss(paste0(data.path,"Round2_RSMGL/",alts[j],"/transect_flows.dss"))  
  test=pathsToDataFrame(getCatalogedPathnames(dss_out), simplify=T)
  for(i in 1:length(RSM.sites)){
    
    if(nrow(subset(test,LOCATION==RSM.sites[i]))>0){
      paths=paste0("/RSMGL/",RSM.sites[i],"/OLFLOW//1DAY/SIMULATED/")  
      tmp=data.frame(getFullTSC(dss_out,paths))
      tmp$Date=date.fun(rownames(tmp))
      rownames(tmp)<-NULL
      tmp$transect=RSM.sites[i]
      tmp$Alt=alts[j]
      trans_flow.ol=rbind(tmp,trans_flow.ol)
      print(i)
    }else{
      next
    }
  }
}
trans_flow.ol$Alt=factor(trans_flow.ol$Alt,levels=alts)
trans_flow.ol$month=as.numeric(format(trans_flow.ol$Date,"%m"))
trans_flow.ol$hydro.season=with(trans_flow.ol,ifelse(month%in%seq(6,10,1),"A_Wet","B_Dry"));# FL.Hydroseason(trans_flow.ol$Date)
trans_flow.ol$CY=as.numeric(format(trans_flow.ol$Date,"%Y"))
trans_flow.ol$OLFLOW.kacft=(trans_flow.ol$OLFLOW*2.29569e-5)/1000 ;# no time conversion

# trans_flow.ol$Q.pos=with(trans_flow.ol,ifelse(OLFLOW.kacft>0,OLFLOW.kacft,NA))
# trans_flow.ol$Q.neg=with(trans_flow.ol,ifelse(OLFLOW.kacft<0,OLFLOW.kacft,NA))

trans.flow.ann=ddply(trans_flow.ol,c("Alt","CY","transect"),summarise,TFlow.kacft=sum(OLFLOW.kacft,na.rm=T))
trans.flow.ann=dcast(trans.flow.ann,Alt~transect,value.var="TFlow.kacft",mean)

trans.flow.seasonal=ddply(trans_flow.ol,c("Alt","CY","hydro.season","transect"),summarise,TFlow.kacft=sum(OLFLOW.kacft,na.rm=T))
trans.flow.seasonal=dcast(trans.flow.seasonal,Alt+transect~hydro.season,value.var="TFlow.kacft",mean)

trans.flow.month=ddply(trans_flow.ol,c("Alt","CY","month","transect"),summarise,TFlow.kAcft=sum(OLFLOW.kacft,na.rm=T))
trans.flow.month=ddply(trans.flow.month,c("Alt","month","transect"),summarise,mean.val=mean(TFlow.kAcft,na.rm=T))

# trans.flow.ann=ddply(trans_flow.ol,c("Alt","CY","transect"),summarise,TFlow.pos=sum(Q.pos,na.rm=T),TFlow.neg=sum(Q.neg,na.rm=T))
# trans.flow.ann.pos=dcast(trans.flow.ann,Alt~transect,value.var="TFlow.pos",mean)
# trans.flow.ann.neg=dcast(trans.flow.ann,Alt~transect,value.var="TFlow.neg",mean)
# 
# trans.flow.seasonal=ddply(trans_flow.ol,c("Alt","CY","hydro.season","transect"),summarise,TFlow.pos=sum(Q.pos,na.rm=T),TFlow.neg=sum(Q.neg,na.rm=T))
# trans.flow.seasonal.pos=dcast(trans.flow.seasonal,Alt+transect~hydro.season,value.var="TFlow.pos",mean)
# trans.flow.seasonal.neg=dcast(trans.flow.seasonal,Alt+transect~hydro.season,value.var="TFlow.neg",mean)
RSM.sites.labs=paste0("T23",LETTERS[1:3])
# png(filename=paste0(plot.path,"Round2/T23_AnnMean.png"),width=6,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.5,0.5),oma=c(2,3,1,0.5),lwd=0.5);
layout(matrix(1:3,3,1))

for(i in 1:3){
ylim.val=c(-100,200);by.y=100;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot(trans.flow.ann[,RSM.sites[i]],
          col=adjustcolor(cols.alts,0.5),
          border=cols.alts,
          space=0.05,
          ylim=ylim.val,axes=F,ann=F)
if(i==1){mtext(side=3,adj=0,"Overland Flow")}
abline(h=0,lwd=1)
text(x,trans.flow.ann[,RSM.sites[i]],round(trans.flow.ann[,RSM.sites[i]],0),pos=ifelse(round(trans.flow.ann[,RSM.sites[i]],0)<0,1,3),offset=0.1,cex=0.75)
if(i==3){axis_fun(1,x,x,alts,line=-0.5,cex=0.8)}else{axis_fun(1,x,x,NA,line=-0.5,cex=0.8)}
axis_fun(2,ymaj,ymin,format(ymaj));
box(lwd=1)
mtext(side=3,adj=0,line=-1.25,paste0(" ",RSM.sites.labs[i]))
}

mtext(side=2,outer=T,line=1,"Average Annual Net Discharge Volume (x10\u2074 AcFt)")
mtext(side=1,outer=T,line=1,"Alternatives")
dev.off()


cols=c("lightblue","khaki3")
# png(filename=paste0(plot.path,"Round2/T23_Seasonal_AnnMean.png"),width=6,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.5,0.5),oma=c(2,3,1,0.5),lwd=0.5);
layout(matrix(1:3,3,1))

ylim.val=c(-100,200);by.y=100;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:length(RSM.sites)){
tmp=t(subset(trans.flow.seasonal,transect==RSM.sites[i])[,3:4])
x=barplot(tmp,col=cols,border="grey",
          space=0.05,
          ylim=ylim.val,axes=F,ann=F,names.arg = rep(NA,n.alts))
if(i==1){mtext(side=3,adj=0,"Overland Flow")}
abline(h=0,lwd=1)
text(x,colSums(tmp),format(round(colSums(tmp),1)),offset=0.1,pos=ifelse(colSums(tmp)<0,1,3))
text(x,tmp[1,]+tmp[2,]/2,format(round(tmp[2,],1)),font=3,cex=0.75)
text(x,tmp[1,]/2,format(round(tmp[1,],1)),font=3,cex=0.75)
if(i==3){axis_fun(1,x,x,alts,line=-0.5,cex=0.8)}else{axis_fun(1,x,x,NA,line=-0.5,cex=0.8)}
axis_fun(2,ymaj,ymin,format(ymaj));
box(lwd=1)
mtext(side=3,adj=0,line=-1.25,paste0(" ",RSM.sites.labs[i]))
if(i==1){
  legend("topright",legend=c("Wet (June - Oct)","Dry (Nov - May)"),
         lty=c(0),lwd=c(0.1),col=c(cols),pch=22,pt.bg=cols,
         ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
}
  
}
mtext(side=2,outer=T,line=1,"Average Annual Net Discharge Volume (x10\u2074 AcFt)")
mtext(side=1,outer=T,line=1,"Alternatives")
dev.off()


# png(filename=paste0(plot.path,"Round2/T23_MonthlyMean.png"),width=5,height=5.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.5,0.5),oma=c(2,3,1,0.5),lwd=0.5);
layout(matrix(1:3,3,1))

xlim.val=c(1,12);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(-30,30);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(j in 1:length(RSM.sites)){
plot(0:1,0:1,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 4:n.alts){
  # ylim.val=c(-0.15,9.5);by.y=1;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
  with(subset(trans.flow.month,Alt==alts[1]&transect==RSM.sites[j]),lines(month,mean.val,col=cols.alts[1],lty=1,lwd=1.5))
  with(subset(trans.flow.month,Alt==alts[2]&transect==RSM.sites[j]),lines(month,mean.val,col=cols.alts[2],lty=2,lwd=1.5))
  with(subset(trans.flow.month,Alt==alts[3]&transect==RSM.sites[j]),lines(month,mean.val,col=cols.alts[3],lty=1,lwd=1.5))
  with(subset(trans.flow.month,Alt==alts[i]&transect==RSM.sites[j]),lines(month,mean.val,col=cols.alts[i],lty=1,lwd=1.5))
}
if(j==1){mtext(side=3,adj=0,"Overland Flow")}
abline(h=0,lwd=1)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
if(j==1){
legend("bottomleft",legend=c(alts),
       lty=c(1,2,1,1,1,1,1),lwd=c(1.5),col=c(cols.alts),
       ncol=3,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
}
if(j==3){axis_fun(1,xmaj,xmin,month.abb[xmaj])}else{axis_fun(1,xmaj,xmin,NA)}
mtext(side=3,adj=0,line=-1.25,paste0(" ",RSM.sites.labs[j]))
}

mtext(side=2,outer=T,line=1,"Average Monthly Net Discharge Volume (x10\u2074 AcFt)")
mtext(side=1,line=2,"Month")
dev.off()

# png(filename=paste0(plot.path,"Round2/T24_daily_QDC.png"),width=7,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,0.75,0.5,1),oma=c(2.5,3.75,1,0.25));
layout(matrix(1:12,3,4,byrow=T))

xlim.val=c(0,1);by.x=0.5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(-2,6);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(j in 1:length(RSM.sites)){
  for(i in 4:n.alts){
    plot(0:1,0:1,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
    abline(h=ymaj,v=xmaj,lty=3,col="grey")
    with(ecdf_fun(subset(trans_flow.ol,transect==RSM.sites[j]&Alt==alts[1])$OLFLOW.kacft),lines(1-proportion,value,col=cols.alts[1],lty=1,lwd=1.5))
    with(ecdf_fun(subset(trans_flow.ol,transect==RSM.sites[j]&Alt==alts[2])$OLFLOW.kacft),lines(1-proportion,value,col=cols.alts[2],lty=2,lwd=1.5))
    with(ecdf_fun(subset(trans_flow.ol,transect==RSM.sites[j]&Alt==alts[3])$OLFLOW.kacft),lines(1-proportion,value,col=cols.alts[3],lty=1,lwd=1.5))
    with(ecdf_fun(subset(trans_flow.ol,transect==RSM.sites[j]&Alt==alts[i])$OLFLOW.kacft),lines(1-proportion,value,col=adjustcolor(cols.alts[i],0.5),lwd=2))
  abline(h=0,lwd=1)
    if(j==1){
    legend("topright",legend=c(alts[c(1:3,i)]),
         lty=c(1,2,1),lwd=c(1.5,1.5,1.5),col=c(cols.alts[1:3],adjustcolor(as.character(cols.alts[i]),0.5)),
         ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
  }
  if(j==3){axis_fun(1,xmaj,xmin,format(xmaj))}else{axis_fun(1,xmaj,xmin,NA)}
  if(i==4){axis_fun(2,ymaj,ymin,ymaj)}else{axis_fun(2,ymaj,ymin,NA)}
  # if(i==4){mtext(side=3, adj=0,"SRS Total Annual Discharge")}
  if(i==4){mtext(side=3,adj=0,line=-1.25,paste0(" ",RSM.sites.labs[j]))}
  box(lwd=1)
  if(i==4&j==1){mtext(side=3,adj=0,"Overland Flow")}
  }
}
mtext(side=2,line=1.5,outer=T,"Net Daily Discharge Volume (x10\u2074 AcFt)")
mtext(side=1,line=1,outer=T,"Proportion of Time \u2265 Discharge")
dev.off()


ecdf.rslts=data.frame()
prop.val=seq(0.1,0.9,0.1)
for(k in 1:length(RSM.sites)){
for(i in 1:n.alts){
test=ecdf_fun(subset(trans_flow.ol,transect==RSM.sites[k]&Alt==alts[i])$OLFLOW.kacft)
test$proportion=1-test$proportion

for(j in 1:length(prop.val)){
  tmp=min(subset(test,proportion<prop.val[j])$value)
  rslt=data.frame(Alt=alts[i],transect=RSM.sites[k],prop=prop.val[j],value=tmp)
  ecdf.rslts=rbind(ecdf.rslts,rslt)
}
}
}
tmp=dcast(ecdf.rslts,transect+prop~Alt,value.var = "value",mean)

tmp=subset(ecdf.rslts,Alt=="FWOi")[,c("prop","transect","value")]
colnames(tmp)=c("prop","transect","FWOi")

ecdf.rslts2=merge(subset(ecdf.rslts,Alt!="FWOi"),
                  tmp,c("prop","transect"))
ecdf.rslts2$diff=with(ecdf.rslts2,FWOi-value)

k=i
for(k in 1:length(RSM.sites)){
png(filename=paste0(plot.path,"Round2/",RSM.sites.labs[k],"_daily_QDCDiff_FWOi.png"),width=7,height=4.5,units="in",res=200,type="windows",bg="white")
  par(family="serif",mar=c(1,0.75,1,1),oma=c(2.5,4,1,0.25));
  layout(matrix(1:8,2,4,byrow=T))

xlim.val=c(0,1);by.x=0.5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
for(i in 4:n.alts){
ylim.val=c(-2,6);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(0:1,0:1,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(ecdf_fun(subset(trans_flow.ol,transect==RSM.sites[k]&Alt==alts[2])$OLFLOW.kacft),lines(1-proportion,value,col=cols.alts[2],lty=2,lwd=1.5))
with(ecdf_fun(subset(trans_flow.ol,transect==RSM.sites[k]&Alt==alts[i])$OLFLOW.kacft),lines(1-proportion,value,col=cols.alts[i],lty=1,lwd=1.5))
abline(h=0,lwd=1)
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
if(i==4){axis_fun(2,ymaj,ymin,ymaj)}else{axis_fun(2,ymaj,ymin,NA)}
box(lwd=1)
# mtext(side=3,adj=0,line=-1.25,paste0(" ",alts[i]))
if(i==4){mtext(side=2,line=2,outer=F,"Net Daily Discharge Volume\n(x10\u2074 AcFt)")}
  legend("topright",legend=c(alts[c(2,i)]),
         lty=c(2,1),lwd=c(1.5,1.5),col=c(cols.alts[2],as.character(cols.alts[i])),
         ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
if(i==4){mtext(side=3,adj=0,paste("  ",RSM.sites.labs[k]),line=-1.25)}
  if(i==4){mtext(side=3,adj=0,"Overland Flow")}
}

for(i in 4:n.alts){
ylim.val=c(-1,0);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(0:1,0:1,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(subset(ecdf.rslts2,transect==RSM.sites[k]&Alt==alts[i]),pt_line(prop,diff,2,"firebrick",1,21,"firebrick1",pt.lwd=0.01,cex=1.25))
with(subset(ecdf.rslts2,transect==RSM.sites[k]&Alt==alts[i]),text(prop,diff,round(diff,2),pos=ifelse(round(diff,1)>-0.2,1,3),cex=0.5))
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
if(i==4){axis_fun(2,ymaj,ymin,format(ymaj))}else{axis_fun(2,ymaj,ymin,NA)}
box(lwd=1)
if(i==4){mtext(side=2,line=2,"\u0394 Net Discharge\n(x10\u2074 AcFt)")}
}
mtext(side=1,line=1,outer=T,"Proportion of Time \u2265 Discharge")
dev.off()
}



# GW flow - transects -----------------------------------------------------
RSM.sites=c(paste0("T23",LETTERS[1:3],"_TRANSECT"))

trans_flow.gw=data.frame()
for(j in 1:n.alts){
  dss_out=opendss(paste0(data.path,"Round2_RSMGL/",alts[j],"/transect_flows.dss"))  
  test=pathsToDataFrame(getCatalogedPathnames(dss_out), simplify=T)
  for(i in 1:length(RSM.sites)){
    
    if(nrow(subset(test,LOCATION==RSM.sites[i]))>0){
      paths=paste0("/RSMGL/",RSM.sites[i],"/GWFLOW//1DAY/SIMULATED/")  
      tmp=data.frame(getFullTSC(dss_out,paths))
      tmp$Date=date.fun(rownames(tmp))
      rownames(tmp)<-NULL
      tmp$transect=RSM.sites[i]
      tmp$Alt=alts[j]
      trans_flow.gw=rbind(tmp,trans_flow.gw)
      print(i)
    }else{
      next
    }
  }
}
trans_flow.gw$Alt=factor(trans_flow.gw$Alt,levels=alts)
trans_flow.gw$month=as.numeric(format(trans_flow.gw$Date,"%m"))
trans_flow.gw$hydro.season=with(trans_flow.gw,ifelse(month%in%seq(6,10,1),"A_Wet","B_Dry"));# FL.Hydroseason(trans_flow.gw$Date)
trans_flow.gw$CY=as.numeric(format(trans_flow.gw$Date,"%Y"))
trans_flow.gw$GWFLOW.kacft=(trans_flow.gw$GWFLOW*2.29569e-5)/1000 ;# no time conversion

# trans_flow.gw$Q.pos=with(trans_flow.gw,ifelse(OLFLOW.kacft>0,OLFLOW.kacft,NA))
# trans_flow.gw$Q.neg=with(trans_flow.gw,ifelse(OLFLOW.kacft<0,OLFLOW.kacft,NA))

trans.flow.gw.ann=ddply(trans_flow.gw,c("Alt","CY","transect"),summarise,TFlow.kacft=sum(GWFLOW.kacft,na.rm=T))
trans.flow.gw.ann=dcast(trans.flow.gw.ann,Alt~transect,value.var="TFlow.kacft",mean)

trans.flow.gw.seasonal=ddply(trans_flow.gw,c("Alt","CY","hydro.season","transect"),summarise,TFlow.kacft=sum(GWFLOW.kacft,na.rm=T))
trans.flow.gw.seasonal=dcast(trans.flow.gw.seasonal,Alt+transect~hydro.season,value.var="TFlow.kacft",mean)

trans.flow.gw.month=ddply(trans_flow.gw,c("Alt","CY","month","transect"),summarise,TFlow.kAcft=sum(GWFLOW.kacft,na.rm=T))
trans.flow.gw.month=ddply(trans.flow.gw.month,c("Alt","month","transect"),summarise,mean.val=mean(TFlow.kAcft,na.rm=T))



cols=c("lightblue","khaki3")
# png(filename=paste0(plot.path,"Round2/T23_Seasonal_AnnMean_GW.png"),width=6,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.5,0.5),oma=c(2,3,1,0.5),lwd=0.5);
layout(matrix(1:3,3,1))

# ylim.val=c(-5,15);by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
ylim.val.min=c(-0.5,0,-4)
ylim.val.max=c(0,5,12)
by.y.val=c(0.2,1,4)
for(i in 1:length(RSM.sites)){
  ylim.val=c(ylim.val.min[i],ylim.val.max[i]);by.y=by.y.val[i];ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  tmp=t(subset(trans.flow.gw.seasonal,transect==RSM.sites[i])[,3:4])
  x=barplot(tmp,col=cols,border="grey",
            space=0.05,
            ylim=ylim.val,axes=F,ann=F,names.arg = rep(NA,n.alts))
  if(i==1){mtext(side=3,adj=0,"Groundwater Flow")}
  abline(h=0,lwd=1)
  text(x,colSums(tmp),format(round(colSums(tmp),2)),offset=0.1,pos=ifelse(colSums(tmp)<0,1,3))
  text(x,tmp[1,]+tmp[2,]/2,format(round(tmp[2,],2)),font=3,cex=0.75)
  text(x,tmp[1,]/2,format(round(tmp[1,],2)),font=3,cex=0.75)
  if(i==3){axis_fun(1,x,x,alts,line=-0.5,cex=0.8)}else{axis_fun(1,x,x,NA,line=-0.5,cex=0.8)}
  axis_fun(2,ymaj,ymin,format(ymaj));
  box(lwd=1)
  mtext(side=3,adj=0,line=-1.25,paste0(" ",RSM.sites.labs[i]))
  if(i==1){
    legend("bottomright",legend=c("Wet (June - Oct)","Dry (Nov - May)"),
           lty=c(0),lwd=c(0.1),col=c(cols),pch=22,pt.bg=cols,
           ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
  }
  
}
mtext(side=2,outer=T,line=1,"Average Annual Net Discharge Volume (x10\u2074 AcFt)")
mtext(side=1,outer=T,line=1,"Alternatives")
dev.off()
