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


## Paths
wd="C:/Julian_LaCie/_GitHub/BBSEER_WQ"

paths=paste0(wd,c("/Plots/","/Export/","/Data/","/src/"))
# Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]

## Functions
consec.startend=function(var){
  runs=rle(var)
  myruns = which(runs$values == TRUE)
  runs.lengths.cumsum = cumsum(runs$lengths)
  ends = runs.lengths.cumsum[myruns]
  newindex = ifelse(myruns>1, myruns-1, 0)
  starts = runs.lengths.cumsum[newindex] + 1
  if (0 %in% newindex) starts = c(1,starts)
  rslt=list(starts=starts,ends=ends)
  return(rslt)
}
poly.step=function(x,y.L, y.U, bg, col = bg, lty = 3, col.adj = 0.25, 
                   lwd = 1,...){
  n=length(x)
  xs <- rep(1:n, each = 2)[-2*n]
  ys <- c(1, rep(2:n, each = 2))
  xx <- x[xs]
  y.low <- y.L[ys]
  y.hi <- y.U[ys]
  
  xx=c(xx,rev(xx))
  yy=c(y.low,rev(y.hi))
  polygon(xx, yy, col = adjustcolor(bg, col.adj), border = col, 
          lty = lty, lwd = lwd)
}


# -------------------------------------------------------------------------
alts=c("FWO","FWOi","ECB22",paste0("ALT",c(21,22,23,24)))
n.alts=length(alts)

cols.alts=c("grey50","grey50","grey10",wesanderson::wes_palette("Zissou1",4,"continuous"))

dss_out=opendss(paste0(data.path,"Round2_RSMGL/",alts[3],"/RSMGL_output.dss"))
test=pathsToDataFrame(getCatalogedPathnames(dss_out), simplify=T)
subset(test,LOCATION=="S333")
subset(test,LOCATION=="S200")

subset(test,LOCATION=="S631")

subset(test,LOCATION=="S333_US")



dss_out=opendss(paste0(data.path,"Round2_RSMGL/",alts[2],"/RSMGL_output.dss"))
test=pathsToDataFrame(getCatalogedPathnames(dss_out), simplify=T)

test[grep("200",test$LOCATION),]

# SRS ---------------------------------------------------------------------
RSM.sites=c(paste0("S12",c("A","B","C","D")),"S333","S333N","S334","S355A","S355B","S356","S335",
            "S151","S631","S632","S633","S31","S337","G211","C111_TO_C111SC","S177","S18C","S178",paste0("S200",LETTERS[1:4]));#L29_DIVIDE from blue shanty east
alt.srs=data.frame()

for(j in 1:n.alts){
  dss_out=opendss(paste0(data.path,"Round2_RSMGL/",alts[j],"/RSMGL_output.dss"))  
  test=pathsToDataFrame(getCatalogedPathnames(dss_out), simplify=T)
  for(i in 1:length(RSM.sites)){
    
    if(nrow(subset(test,LOCATION==RSM.sites[i]))==1){
    paths=paste0("/RSMGL/",RSM.sites[i],"/FLOW//1DAY/SIMULATED/")  
    tmp=data.frame(getFullTSC(dss_out,paths))
    tmp$Date=date.fun(rownames(tmp))
    rownames(tmp)<-NULL
    tmp$SITE=RSM.sites[i]
    tmp$Alt=alts[j]
    alt.srs=rbind(tmp,alt.srs)
    print(i)
    }else{
      next
      }
  }
}

unique(alt.srs$SITE)
unique(alt.srs$Alt)

alt.srs$Q=cfs.to.acftd(alt.srs$FLOW)
alt.srs$Alt=factor(alt.srs$Alt,levels=alts)

alt.srs$CY=as.numeric(format(alt.srs$Date,'%Y'))

alt.srs.CY=ddply(alt.srs,c("CY","Alt","SITE"),summarise,
                 TFlow.acft=sum(Q))
# alt.srs.CY=ddply(alt.srs.CY,c("Alt","SITE"),summarise,
#                  mean.TFlow.acft=mean(TFlow.acft))

alt.srs.CY=dcast(alt.srs.CY,Alt~SITE,value.var = "TFlow.acft",mean)

ylim.val=c(0,200e3);by.y=50e3;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)

RSM.sites2=c(paste0("S12",c("A","B","C","D")),"S333","S333N","S355A","S355B","S356")
# png(filename=paste0(plot.path,"Round2/SRS_AnnMean.png"),width=6.5,height=5.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.5,0.5),oma=c(2,3,0.5,0.5),lwd=0.5);
layout(matrix(1:10,5,2))

ylim.val.max=c(50e3,60e3,150e3,200e3,600e3,80e3,18e3,18e3,200e3)
by.y.val=ylim.val.max/4

for(i in 1:4){
  ylim.val=c(0,ylim.val.max[i]);by.y=by.y.val[i];ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  x=barplot(alt.srs.CY[,RSM.sites2[i]],
            col=adjustcolor(cols.alts,0.5),
            border=cols.alts,
            space=0.05,
            ylim=ylim.val,axes=F,ann=F)
  if(i==4){axis_fun(1,x,x,alts,line=-0.5,cex=0.8)}else{axis_fun(1,x,x,NA,line=-0.5)}
  axis_fun(2,ymaj,ymin,format(ymaj/1e3));
  box(lwd=1)
  mtext(side=3,adj=0,line=-1.25,paste0(" ",RSM.sites2[i]))
  text(x,alt.srs.CY[,RSM.sites2[i]],format(round(alt.srs.CY[,RSM.sites2[i]]/1000,1)),
       pos=1,offset=0.1,cex=0.75)
}
plot(0:1,0:1,ann=F,axes=F,type="n")

for(i in 5:9){
  ylim.val=c(0,ylim.val.max[i]);by.y=by.y.val[i];ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  x=barplot(alt.srs.CY[,RSM.sites2[i]],
            col=adjustcolor(cols.alts,0.5),
            border=cols.alts,
            space=0.05,
            ylim=ylim.val,axes=F,ann=F)
  if(i==9){axis_fun(1,x,x,alts,line=-0.5,cex=0.8)}else{axis_fun(1,x,x,NA,line=-0.5)}
  axis_fun(2,ymaj,ymin,format(ymaj/1e3));
  box(lwd=1)
  mtext(side=3,adj=0,line=-1.25,paste0(" ",RSM.sites2[i]))
  text(x,alt.srs.CY[,RSM.sites2[i]],format(round(alt.srs.CY[,RSM.sites2[i]]/1000,1)),
       pos=ifelse(round(alt.srs.CY[,RSM.sites2[i]]/1000,1)<10,3,1),
       offset=0.1,cex=0.75)
}
mtext(side=2,outer=T,line=1,"Average Annual Discharge Volume (x1000 AcFt Y\u207B\u00B9)")
mtext(side=1,outer=T,"Alternatives")
dev.off()

flow.xtab=reshape2::dcast(alt.srs,Alt+Date~SITE,value.var="Q",mean)
flow.xtab2=flow.xtab
# flow.xtab$FedWY=WY(flow.xtab$Date,"Fed")
# flow.xtab=subset(flow.xtab,FedWY%in%seq(1966,2016,1))
dcast(alt.srs,Alt~SITE,value.var="Q",min)
dcast(alt.srs,Alt~SITE,value.var="Q",max)

dcast(alt.srs,Alt~SITE,value.var="FLOW",max)

## Settlement Agreement Calculation - S12s+[S333s + S355s + min(S356,S335)-S334]
## Adjsuted values are needed to assign flow to each respective strcuture for WQ and FWM calculation purposes
flow.xtab$minS356_S335=apply(flow.xtab[,c('S356','S335')],1,min);
flow.xtab$SA.ESRS=rowSums(flow.xtab[,c("S333","S333N","S355A","S355B","minS356_S335")],na.rm=T)
flow.xtab$S333R=with(flow.xtab,S333*(1-ifelse(SA.ESRS>0,ifelse(S334/SA.ESRS>1,1,S334/SA.ESRS),0)))
flow.xtab$S333NR=with(flow.xtab,S333N*(1-ifelse(SA.ESRS>0,ifelse(S334/SA.ESRS>1,1,S334/SA.ESRS),0)))
flow.xtab$S355A.adj=with(flow.xtab,S355A*(1-ifelse(SA.ESRS>0,ifelse(S334/SA.ESRS>1,1,S334/SA.ESRS),0)))
flow.xtab$S355B.adj=with(flow.xtab,S355B*(1-ifelse(SA.ESRS>0,ifelse(S334/SA.ESRS>1,1,S334/SA.ESRS),0)))
flow.xtab$S356.adj=with(flow.xtab,minS356_S335*(1-ifelse(SA.ESRS>0,ifelse(S334/SA.ESRS>1,1,S334/SA.ESRS),0)))
flow.xtab$TFlow.SA=rowSums(flow.xtab[,c("S12A","S12B","S12C","S12D","S333","S355A","S355B","minS356_S335")],na.rm=T)

flow.xtab2$total.ESRS=rowSums(flow.xtab2[,c("S333","S333N","S355A","S355B","S356")],na.rm=T); #total flow into SRS
flow.xtab2$total.ESRS.R=with(flow.xtab2,ifelse(total.ESRS-S334<0,0,total.ESRS-S334))
flow.xtab2$total.WSRS=rowSums(flow.xtab2[,paste0("S12",LETTERS[1:4])],na.rm=T);
flow.xtab2$CY=as.numeric(format(flow.xtab2$Date,"%Y"))




flow.xtab2.CY=ddply(flow.xtab2,c("CY","Alt"),summarise,
                    sum.ESRS=sum(total.ESRS.R,na.rm=T),
                    sum.WSRS=sum(total.WSRS,na.rm=T))
flow.xtab2.CY$Tflow=rowSums(flow.xtab2.CY[,c("sum.ESRS","sum.WSRS")],na.rm=T)

flow.xtab2.CY2=ddply(flow.xtab2.CY,c("Alt"),summarise,
                     mean.ESRS=mean(sum.ESRS),
                     mean.WSRS=mean(sum.WSRS),
                     mean.Tflow=mean(Tflow))
flow.xtab2.CY2$Per.ESRS=with(flow.xtab2.CY2,(mean.ESRS/mean.Tflow)*100)
flow.xtab2.CY2$Per.WSRS=with(flow.xtab2.CY2,(mean.WSRS/mean.Tflow)*100)



cols=wesanderson::wes_palette("Zissou1",2,"continuous")#viridis::plasma(3)[1:2]
# png(filename=paste0(plot.path,"Round2/SRS_flowdist.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
ylim.val=c(0,100);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)

par(family="serif",mar=c(2,2,0.5,1),oma=c(1,2,0.5,0.5),lwd=0.5);
layout(matrix(1:2,2,1),heights=c(1,0.4))

x=barplot(t(flow.xtab2.CY2[,c("Per.ESRS","Per.WSRS")]),
          col=adjustcolor(cols,0.5),
          border=cols,
          space=0.05,
          ylim=ylim.val,axes=F,ann=F)
with(flow.xtab2.CY2,text(x,(Per.ESRS/2),paste0(round(Per.ESRS,1),"%","\n(",round(mean.ESRS/1000),")"),cex=0.75))
with(flow.xtab2.CY2,text(x,((Per.ESRS-Per.WSRS)/2)+Per.ESRS,paste0(round(Per.WSRS,1),"%","\n(",round(mean.WSRS/1000),")"),cex=0.75))
axis_fun(1,x,x,alts,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Percent Distribution")
mtext(side=1,line=2,"Alternative")


par(mar=c(0,2,0,1))
plot(0:1,0:1,ann=F,axes=F,type="n")
legend("center",legend=c("East Shark River Slough (S333s + S355s + S356 - S334)","West Shark River Slough (S12s)"),
       pch=22,lty=0,lwd=0.1,pt.bg=adjustcolor(cols,0.5),pt.cex=2,col=cols,
       ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=1,xpd=NA,xjust=0,yjust=1,
       title.adj = 0,title=" Average Annual Flow Distribution")
mtext(side=1,line=-1,adj=0,"Percentage (Annual Average Discharge Volume x1000 AcFt Y\u207B\u00B9)",font=3,cex=0.8)
dev.off()


plot(ecdf(subset(flow.xtab2.CY,Alt=="FWO")$Tflow))
plot(ecdf(subset(flow.xtab2.CY,Alt=="FWOi")$Tflow),add=T)
plot(ecdf(subset(flow.xtab2.CY,Alt=="ECB22")$Tflow),add=T)

# png(filename=paste0(plot.path,"Round2/SRS_FDC.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,0.75,0.25,1),oma=c(2.5,3.75,1,0.25));
layout(matrix(1:4,1,4,byrow=T))

xlim.val=c(0,1);by.x=0.5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,1400e3);by.y=400e3;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)

for(i in 4:n.alts){
  plot(0:1,0:1,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(ecdf_fun(subset(flow.xtab2.CY,Alt==alts[1])$Tflow),lines(1-proportion,value,col=cols.alts[1],lty=1,lwd=1.5))
  with(ecdf_fun(subset(flow.xtab2.CY,Alt==alts[2])$Tflow),lines(1-proportion,value,col=cols.alts[2],lty=2,lwd=1.5))
  with(ecdf_fun(subset(flow.xtab2.CY,Alt==alts[3])$Tflow),lines(1-proportion,value,col=cols.alts[3],lty=1,lwd=1.5))
  with(ecdf_fun(subset(flow.xtab2.CY,Alt==alts[i])$Tflow),lines(1-proportion,value,col=adjustcolor(cols.alts[i],0.5),lwd=2))
  
  legend("topright",legend=c(alts[c(1:3,i)]),
         lty=c(1,2,1),lwd=c(1.5,1.5,1.5),col=c(cols.alts[1:3],adjustcolor(as.character(cols.alts[i]),0.5)),
         ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
  axis_fun(1,xmaj,xmin,format(xmaj))
  if(i==4){axis_fun(2,ymaj,ymin,ymaj/1e3)}else{axis_fun(2,ymaj,ymin,NA)}
  if(i==4){mtext(side=3, adj=0,"SRS Total Annual Discharge")}
  box(lwd=1)
}
mtext(side=2,line=1.5,outer=T,"Discharge (x1000 AcFt Y\u207B\u00B9)")
mtext(side=1,line=1,outer=T,"Proportion of Time \u2265 Discharge")
dev.off()


# S356 --------------------------------------------------------------------

S356.q=subset(alt.srs,SITE%in%c("S356"))
S356.q$hydro.season=FL.Hydroseason(S356.q$Date)
S356.q$month=as.numeric(format(S356.q$Date,"%m"))
# wca3b.flow.xtab=dcast(wca3b.q,Alt+CY~SITE,value.var = "Q",sum)

S356.flow.ann=ddply(S356.q,c("Alt","CY","SITE"),summarise,TFlow.kAcftY=sum(Q/1000,na.rm=T))
S356.flow.ann=dcast(S356.flow.ann,Alt~SITE,value.var="TFlow.kAcftY",mean)

S356.flow.seasonal=ddply(S356.q,c("Alt","CY","hydro.season","SITE"),summarise,TFlow.kAcftY=sum(Q/1000,na.rm=T))
S356.flow.seasonal=dcast(S356.flow.seasonal,Alt+SITE~hydro.season,value.var="TFlow.kAcftY",mean)

S356.flow.month=ddply(S356.q,c("Alt","CY","month","SITE"),summarise,TFlow.kAcftY=sum(Q/1000,na.rm=T))
S356.flow.month=ddply(S356.flow.month,c("Alt","month","SITE"),summarise,mean.val=mean(TFlow.kAcftY,na.rm=T))

S356.flow.month.avg=ddply(S356.q,c("Alt","month","SITE"),summarise,mean.cfs=mean(FLOW,na.rm=T))

ecdf.fun2=function(Q){
  Q.sort=sort(Q,decreasing=T)
  df=data.frame(proportion=(1/length(Q.sort))*1:length(Q.sort),value=Q.sort)
  return(df)
}

Q.val=subset(S356.q,Alt==alts[i])$Q
Q.sort=sort(Q.val,decreasing = T)
prop=(100/length(Q.val))*1:length(Q.val)/100

cols2=c("lightblue","khaki3")
# png(filename=paste0(plot.path,"Round2/S356_ann_month_qdf.png"),width=8,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.5,1),oma=c(3,2,0.75,0.25),lwd=0.5);
layout(matrix(1:3,1,3))

ylim.val=c(0,200);by.y=50;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
tmp=t(S356.flow.seasonal[,3:4])
x=barplot(t(S356.flow.seasonal[,3:4]),col=cols2,border="grey",
          space=0.05,
          ylim=ylim.val,axes=F,ann=F)
text(x,colSums(tmp),format(round(colSums(tmp),1)),offset=0.1,pos=3,cex=0.75)
text(x,tmp[1,]+tmp[2,]/2,format(round(tmp[2,],1)),font=3,cex=0.75)
text(x,tmp[1,]/2,format(round(tmp[1,],1)),font=3,cex=0.75)
axis_fun(1,x,x,alts,line=-0.5,cex=0.8,las=2)
axis_fun(2,ymaj,ymin,format(ymaj));
box(lwd=1)
mtext(side=3,adj=0," S-356")
mtext(side=2,outer=F,line=2.25,"Avg. Annual Discharge Volume (kAcFt Y\u207B\u00B9)",cex=0.8)
mtext(side=1,outer=F,line=2.6,"Alternatives")
legend("topright",legend=c("Wet (May - Oct)","Dry (Nov - April)"),
       lty=c(0),lwd=c(0.1),col="grey",pch=22,pt.bg=cols2,pt.cex=1.25,
       ncol=1,cex=0.75,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)


xlim.val=c(1,12);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0,400);by.y=100;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
plot(0:1,0:1,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
# ylim.val=c(-0.15,9.5);by.y=1;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
with(subset(S356.flow.month.avg,Alt==alts[1]),lines(month,mean.cfs,col=cols.alts[1],lty=1,lwd=1.5))
with(subset(S356.flow.month.avg,Alt==alts[2]),lines(month,mean.cfs,col=cols.alts[2],lty=2,lwd=1.5))
with(subset(S356.flow.month.avg,Alt==alts[3]),lines(month,mean.cfs,col=cols.alts[3],lty=1,lwd=1.5))
for(i in 4:n.alts){
  with(subset(S356.flow.month.avg,Alt==alts[i]),lines(month,mean.cfs,col=cols.alts[i],lty=1,lwd=1.5))
}
axis_fun(1,xmaj,xmin,month.abb[xmaj])
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,outer=F,line=2.25,"Avg. Daily Discharge (ft\u00B3 sec\u207B\u00B9)",cex=0.8)
mtext(side=1,line=2.6,"Month")
legend("topleft",legend=c(alts),
       lty=c(1,2,1,1,1,1,1),lwd=c(1.5),col=c(cols.alts),
       ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

xlim.val=c(0,1);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,650);by.y=200;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
  # ylim.val=c(-0.15,9.5);by.y=1;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
plot(0:1,0:1,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
lty.val=c(1,2,1)
for(i in 1:3){
  tmp=ecdf.fun2(subset(S356.q,Alt==alts[i])$FLOW)
  #tmp$proportion=1-tmp$proportion
  # tmp$value=exp(tmp$value)
  with(tmp,lines(proportion,value,col=cols.alts[i],lty=lty.val[i],lwd=1.5))
}
for(i in 4:n.alts){
  with(ecdf_fun(subset(S356.q,Alt==alts[i])$FLOW),lines(1-proportion,value,col=adjustcolor(cols.alts[i],0.5),lwd=2))
}
legend("topright",legend=c(alts),
       lty=c(1,2,1,1,1,1,1),lwd=c(1.5),col=c(cols.alts),
       ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
axis_fun(1,xmaj,xmin,format(xmaj))
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,outer=F,line=2.5,"Daily Discharge (ft\u00B3 sec\u207B\u00B9)",cex=0.8)
mtext(side=1,line=2.6,"Proportion of Time \u2265 Discharge")
dev.off()

tmp=ecdf_fun(log(subset(S356.q,Alt==alts[1])$Q/1000))
tmp$proportion=1-tmp$proportion
plot(proportion~exp(value),tmp)

# WCA3B flows -------------------------------------------------------------

wca3b.q=subset(alt.srs,SITE%in%c("S151","S631","S632","S633","S31","S337","S355A","S355B"))
# wca3b.flow.xtab=dcast(wca3b.q,Alt+CY~SITE,value.var = "Q",sum)

wca3b.flow.ann=ddply(wca3b.q,c("Alt","CY","SITE"),summarise,TFlow.kAcftY=sum(Q/1000,na.rm=T))
wca3b.flow.ann.avg=dcast(wca3b.flow.ann,Alt~SITE,value.var="TFlow.kAcftY",mean)


wca3b.flow.ann.avg[is.na(wca3b.flow.ann.avg)==T]=0
RSM.sites2=c("S151","S631","S632","S633","S31","S337","S355A","S355B")
# png(filename=paste0(plot.path,"Round2/WCA3B_AnnMean.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.5,0.5),oma=c(2,3,1,0.5),lwd=0.5);
layout(matrix(1:8,4,2))

for(i in 1:4){
  ylim.val=c(0,200);by.y=50;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  x=barplot(wca3b.flow.ann.avg[,RSM.sites2[i]],
            col=adjustcolor(cols.alts,0.5),
            border=cols.alts,
            space=0.05,
            ylim=ylim.val,axes=F,ann=F)
  text(x,wca3b.flow.ann.avg[,RSM.sites2[i]],format(round(wca3b.flow.ann.avg[,RSM.sites2[i]],1)),
       pos=ifelse(round(wca3b.flow.ann.avg[,RSM.sites2[i]],1)<10,3,1),
       offset=0.1,cex=0.75)
  if(i==4){axis_fun(1,x,x,alts,line=-0.5,cex=0.8)}else{axis_fun(1,x,x,NA,line=-0.5)}
  axis_fun(2,ymaj,ymin,format(ymaj));
  box(lwd=1)
  mtext(side=3,adj=0,line=-1.25,paste0(" ",RSM.sites2[i]))
  if(i==1){mtext("WCA-3B Inflow")}
}

for(i in 5:8){
  ylim.val=c(0,100);by.y=25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  x=barplot(wca3b.flow.ann.avg[,RSM.sites2[i]],
            col=adjustcolor(cols.alts,0.5),
            border=cols.alts,
            space=0.05,
            ylim=ylim.val,axes=F,ann=F)
  text(x,wca3b.flow.ann.avg[,RSM.sites2[i]],format(round(wca3b.flow.ann.avg[,RSM.sites2[i]],1)),
       pos=ifelse(round(wca3b.flow.ann.avg[,RSM.sites2[i]],1)<10,3,1),
       offset=0.1,cex=0.75)
  if(i==8){axis_fun(1,x,x,alts,line=-0.5,cex=0.8)}else{axis_fun(1,x,x,NA,line=-0.5)}
  axis_fun(2,ymaj,ymin,format(ymaj));
  box(lwd=1)
  mtext(side=3,adj=0,line=-1.25,paste0(" ",RSM.sites2[i]))
  if(i==5){mtext("WCA-3B Outflow")}
}
mtext(side=2,outer=T,line=1,"Average Annual Discharge Volume (x1000 AcFt Y\u207B\u00B9)")
mtext(side=1,outer=T,line=1,"Alternatives")
dev.off()

plot(ecdf(subset(wca3b.flow.xtab,Alt=="FWO")$S631))
plot(ecdf(subset(wca3b.flow.xtab,Alt=="FWOi")$S631),add=T)
plot(ecdf(subset(wca3b.flow.xtab,Alt=="ALT21")$S631),add=T)
plot(ecdf(subset(wca3b.flow.xtab,Alt=="ALT23")$S631),add=T)


# South Dade --------------------------------------------------------------

sd.q=subset(alt.srs,SITE%in%c("G211"))
sd.q$hydro.season=FL.Hydroseason(sd.q$Date)
sd.q$month=as.numeric(format(sd.q$Date,"%m"))
# wca3b.flow.xtab=dcast(wca3b.q,Alt+CY~SITE,value.var = "Q",sum)

sd.flow.ann=ddply(sd.q,c("Alt","CY","SITE"),summarise,TFlow.kAcftY=sum(Q/1000,na.rm=T))
sd.flow.ann=dcast(sd.flow.ann,Alt~SITE,value.var="TFlow.kAcftY",mean)

sd.flow.seasonal=ddply(sd.q,c("Alt","CY","hydro.season","SITE"),summarise,TFlow.kAcftY=sum(Q/1000,na.rm=T))
sd.flow.seasonal=dcast(sd.flow.seasonal,Alt+SITE~hydro.season,value.var="TFlow.kAcftY",mean)

sd.flow.month=ddply(sd.q,c("Alt","CY","month","SITE"),summarise,TFlow.kAcftY=sum(Q/1000,na.rm=T))
sd.flow.month=ddply(sd.flow.month,c("Alt","month","SITE"),summarise,mean.val=mean(TFlow.kAcftY,na.rm=T))

# png(filename=paste0(plot.path,"Round2/SouthDade_G211_MonthlyMean.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.5,0.5),oma=c(2,3,1,0.5),lwd=0.5);

plot(mean.val~month,subset(sd.flow.month,Alt=="FWO"),type="l")
lines(mean.val~month,subset(sd.flow.month,Alt=="FWOi"),type="l")
xlim.val=c(1,12);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0,20);by.y=2;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
plot(0:1,0:1,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 4:n.alts){
  # ylim.val=c(-0.15,9.5);by.y=1;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
  with(subset(sd.flow.month,Alt==alts[1]),lines(month,mean.val,col=cols.alts[1],lty=1,lwd=1.5))
  with(subset(sd.flow.month,Alt==alts[2]),lines(month,mean.val,col=cols.alts[2],lty=2,lwd=1.5))
  with(subset(sd.flow.month,Alt==alts[3]),lines(month,mean.val,col=cols.alts[3],lty=1,lwd=1.5))
  with(subset(sd.flow.month,Alt==alts[i]),lines(month,mean.val,col=cols.alts[i],lty=1,lwd=1.5))
}
legend("topleft",legend=c(alts),
       lty=c(1,2,1,1,1,1,1),lwd=c(1.5),col=c(cols.alts),
       ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
axis_fun(1,xmaj,xmin,month.abb[xmaj])
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,outer=F,line=2,"Average Monthly Discharge Volume\n(x1000 AcFt Y\u207B\u00B9)")
mtext(side=1,line=2,"Month")
mtext(side=3,adj=0,"G211")
dev.off()


cols=c("lightblue","khaki3")
# png(filename=paste0(plot.path,"Round2/SouthDade_G211Seasonal_AnnMean.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.5,0.5),oma=c(2,3,1,0.5),lwd=0.5);
ylim.val=c(0,200);by.y=50;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
tmp=t(sd.flow.seasonal[,3:4])
x=barplot(t(sd.flow.seasonal[,3:4]),col=cols,border=cols.alts,
        space=0.05,
        ylim=ylim.val,axes=F,ann=F)
text(x,colSums(tmp),format(round(colSums(tmp),1)),offset=0.1,pos=3)
text(x,tmp[1,]+tmp[2,]/2,format(round(tmp[2,],1)),font=3,cex=0.75)
text(x,tmp[1,]/2,format(round(tmp[1,],1)),font=3,cex=0.75)
axis_fun(1,x,x,alts,line=-0.5,cex=0.8)
axis_fun(2,ymaj,ymin,format(ymaj));
box(lwd=1)
mtext(side=3,adj=0,line=-1.25,paste0(" ",RSM.sites2[i]))
mtext(side=2,outer=F,line=2,"Average Annual Discharge Volume\n(x1000 AcFt Y\u207B\u00B9)")
mtext(side=1,outer=T,line=1,"Alternatives")
legend("topright",legend=c("Wet (May - Oct)","Dry (Nov - April)"),
       lty=c(0),lwd=c(0.1),col=c(cols),pch=22,pt.bg=cols,
       ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

RSM.sites2=c("G211")
# png(filename=paste0(plot.path,"Round2/SouthDadeQ1_AnnMean.png"),width=6.5,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.5,0.5),oma=c(2,3,1,0.5),lwd=0.5);
# layout(matrix(1:8,4,2))

# for(i in 1:4){
  ylim.val=c(0,200);by.y=50;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  x=barplot(sd.flow.ann[,RSM.sites2[1]],
            col=adjustcolor(cols.alts,0.5),
            border=cols.alts,
            space=0.05,
            ylim=ylim.val,axes=F,ann=F)
  text(x,sd.flow.ann[,RSM.sites2[i]],format(round(sd.flow.ann[,RSM.sites2[i]],1)),
       pos=ifelse(round(sd.flow.ann[,RSM.sites2[i]],1)<10,3,1),
       offset=0.1,cex=0.75)
  axis_fun(1,x,x,alts,line=-0.5,cex=0.8)
  axis_fun(2,ymaj,ymin,format(ymaj));
  box(lwd=1)
  mtext(side=3,adj=0,line=-1.25,paste0(" ",RSM.sites2[i]))
  mtext(side=2,outer=F,line=2,"Average Annual Discharge Volume\n(x1000 AcFt Y\u207B\u00B9)")
  mtext(side=1,outer=T,line=1,"Alternatives")
  dev.off()

  
frogpond.q=subset(alt.srs,SITE%in%c(paste0("S200",LETTERS[1:4])))
ddply(frogpond.q,"SITE",summarise,min.val=min(FLOW),max.val=max(FLOW))
  
par(family="serif",mar=c(1,2,0.5,0.5),oma=c(6,3,1,0.5))
boxplot(TFlow~SITE+Alt,ddply(frogpond.q,c("CY","Alt","SITE"),summarise,TFlow=sum(Q/1000)),las=3,xlab=NA)
mtext(side=2,line=3,"Discharge (x1000 Ac-Ft")

S200s=dcast(frogpond.q,CY+Alt~SITE,value.var = "Q",sum)
  
  
# Stages -------------------------------------------------------------------
RSM.stage.sites=c("S333_US","WCA3A_3GAVG","WCA3B_3B-71","S335U")
stg.dat=data.frame()

for(j in 1:n.alts){
  dss_out=opendss(paste0(data.path,"Round2_RSMGL/",alts[j],"/RSMGL_output.dss"))  
  test=pathsToDataFrame(getCatalogedPathnames(dss_out), simplify=T)
  for(i in 1:length(RSM.stage.sites)){
    
    if(nrow(subset(test,LOCATION==RSM.stage.sites[i]))==1){
      paths=paste0("/RSMGL/",RSM.stage.sites[i],"/STAGE//1DAY/SIMULATED/")  
      tmp=data.frame(getFullTSC(dss_out,paths))
      tmp$Date=date.fun(rownames(tmp))
      rownames(tmp)<-NULL
      tmp$SITE=RSM.stage.sites[i]
      tmp$Alt=alts[j]
      stg.dat=rbind(tmp,stg.dat)
      print(i)
    }else{
      next
    }
  }
}

range(subset(stg.dat,SITE=="WCA3A_3GAVG")$STAGE)
range(subset(stg.dat,SITE=="WCA3B_3B-71")$STAGE)

range(subset(stg.dat,SITE=="S335U")$STAGE)

SDC_seg_3A=ddply(subset(stg.dat,SITE=="WCA3A_3GAVG"),"Alt",summarise,
              SDC_point_10=min(subset(ecdf_fun(STAGE),(1-proportion)<0.10)$value),
              SDC_point_20=min(subset(ecdf_fun(STAGE),1-proportion<0.20)$value),
              SDC_point_30=min(subset(ecdf_fun(STAGE),1-proportion<0.30)$value),
              SDC_point_40=min(subset(ecdf_fun(STAGE),1-proportion<0.40)$value),
              SDC_point_50=min(subset(ecdf_fun(STAGE),1-proportion<0.50)$value),
              SDC_point_60=min(subset(ecdf_fun(STAGE),1-proportion<0.60)$value),
              SDC_point_70=min(subset(ecdf_fun(STAGE),1-proportion<0.70)$value),
              SDC_point_80=min(subset(ecdf_fun(STAGE),1-proportion<0.80)$value),
              SDC_point_90=min(subset(ecdf_fun(STAGE),1-proportion<0.90)$value))
SDC_seg_3B=ddply(subset(stg.dat,SITE=="WCA3B_3B-71"),"Alt",summarise,
              SDC_point_10=min(subset(ecdf_fun(STAGE),(1-proportion)<0.10)$value),
              SDC_point_20=min(subset(ecdf_fun(STAGE),1-proportion<0.20)$value),
              SDC_point_30=min(subset(ecdf_fun(STAGE),1-proportion<0.30)$value),
              SDC_point_40=min(subset(ecdf_fun(STAGE),1-proportion<0.40)$value),
              SDC_point_50=min(subset(ecdf_fun(STAGE),1-proportion<0.50)$value),
              SDC_point_60=min(subset(ecdf_fun(STAGE),1-proportion<0.60)$value),
              SDC_point_70=min(subset(ecdf_fun(STAGE),1-proportion<0.70)$value),
              SDC_point_80=min(subset(ecdf_fun(STAGE),1-proportion<0.80)$value),
              SDC_point_90=min(subset(ecdf_fun(STAGE),1-proportion<0.90)$value))

# png(filename=paste0(plot.path,"Round2/WCA3_STAGE_SDC.png"),width=7,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,0.75,0.5,1),oma=c(2.5,3.75,1,0.25));
layout(matrix(1:8,2,4,byrow=T))

xlim.val=c(0,1);by.x=0.5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
for(i in 4:n.alts){
  ylim.val=c(7,13);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  plot(0:1,0:1,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(ecdf_fun(subset(stg.dat,SITE=="WCA3A_3GAVG"&Alt==alts[1])$STAGE),lines(1-proportion,value,col=cols.alts[1],lty=1,lwd=1.5))
  with(ecdf_fun(subset(stg.dat,SITE=="WCA3A_3GAVG"&Alt==alts[2])$STAGE),lines(1-proportion,value,col=cols.alts[2],lty=2,lwd=1.5))
  with(ecdf_fun(subset(stg.dat,SITE=="WCA3A_3GAVG"&Alt==alts[3])$STAGE),lines(1-proportion,value,col=cols.alts[3],lty=1,lwd=1.5))
  with(ecdf_fun(subset(stg.dat,SITE=="WCA3A_3GAVG"&Alt==alts[i])$STAGE),lines(1-proportion,value,col=adjustcolor(cols.alts[i],0.5),lwd=2))
  
  legend("topright",legend=c(alts[c(1:3,i)]),
         lty=c(1,2,1),lwd=c(1.5,1.5,1.5),col=c(cols.alts[1:3],adjustcolor(as.character(cols.alts[i]),0.5)),
         ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
  axis_fun(1,xmaj,xmin,NA)
  if(i==4){axis_fun(2,ymaj,ymin,ymaj)}else{axis_fun(2,ymaj,ymin,NA)}
  # if(i==4){mtext(side=3, adj=0,"SRS Total Annual Discharge")}
  if(i==4){mtext(side=3,adj=0,"WCA-3A (3 Gauge Avg)",cex=0.9)}
  box(lwd=1)
}
for(i in 4:n.alts){
  ylim.val=c(5,9.5);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  plot(0:1,0:1,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(ecdf_fun(subset(stg.dat,SITE=="WCA3B_3B-71"&Alt==alts[1])$STAGE),lines(1-proportion,value,col=cols.alts[1],lty=1,lwd=1.5))
  with(ecdf_fun(subset(stg.dat,SITE=="WCA3B_3B-71"&Alt==alts[2])$STAGE),lines(1-proportion,value,col=cols.alts[2],lty=2,lwd=1.5))
  with(ecdf_fun(subset(stg.dat,SITE=="WCA3B_3B-71"&Alt==alts[3])$STAGE),lines(1-proportion,value,col=cols.alts[3],lty=1,lwd=1.5))
  with(ecdf_fun(subset(stg.dat,SITE=="WCA3B_3B-71"&Alt==alts[i])$STAGE),lines(1-proportion,value,col=adjustcolor(cols.alts[i],0.5),lwd=2))
  
  legend("topright",legend=c(alts[c(1:3,i)]),
         lty=c(1,2,1),lwd=c(1.5,1.5,1.5),col=c(cols.alts[1:3],adjustcolor(as.character(cols.alts[i]),0.5)),
         ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
  axis_fun(1,xmaj,xmin,format(xmaj))
  if(i==4){axis_fun(2,ymaj,ymin,ymaj)}else{axis_fun(2,ymaj,ymin,NA)}
  # if(i==4){mtext(side=3, adj=0,"SRS Total Annual Discharge")}
  if(i==4){mtext(side=3,adj=0,"WCA-3B 3-71",cex=0.9)}
  box(lwd=1)
}
mtext(side=2,line=1.5,outer=T,"Stage (FT, NGVD29)")
mtext(side=1,line=1,outer=T,"Proportion of Time \u2265 Stage")
dev.off()


# png(filename=paste0(plot.path,"Round2/S335HW_SDC.png"),width=7,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,0.75,0.5,1),oma=c(2.5,3.75,1,0.25));
layout(matrix(1:4,1,4,byrow=T))

xlim.val=c(0,1);by.x=0.5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
for(i in 4:n.alts){
  ylim.val=c(3.5,9.5);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  plot(0:1,0:1,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(ecdf_fun(subset(stg.dat,SITE=="S335U"&Alt==alts[1])$STAGE),lines(1-proportion,value,col=cols.alts[1],lty=1,lwd=1.5))
  with(ecdf_fun(subset(stg.dat,SITE=="S335U"&Alt==alts[2])$STAGE),lines(1-proportion,value,col=cols.alts[2],lty=2,lwd=1.5))
  with(ecdf_fun(subset(stg.dat,SITE=="S335U"&Alt==alts[3])$STAGE),lines(1-proportion,value,col=cols.alts[3],lty=1,lwd=1.5))
  with(ecdf_fun(subset(stg.dat,SITE=="S335U"&Alt==alts[i])$STAGE),lines(1-proportion,value,col=adjustcolor(cols.alts[i],0.5),lwd=2))
  
  legend("topright",legend=c(alts[c(1:3,i)]),
         lty=c(1,2,1),lwd=c(1.5,1.5,1.5),col=c(cols.alts[1:3],adjustcolor(as.character(cols.alts[i]),0.5)),
         ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
  axis_fun(1,xmaj,xmin,format(xmaj))
  if(i==4){axis_fun(2,ymaj,ymin,ymaj)}else{axis_fun(2,ymaj,ymin,NA)}
  # if(i==4){mtext(side=3, adj=0,"SRS Total Annual Discharge")}
  if(i==4){mtext(side=3,adj=0,"S-335 HW",cex=0.9)}
  box(lwd=1)
}
mtext(side=2,line=1.5,outer=T,"Stage (FT, NGVD29)")
mtext(side=1,line=1,outer=T,"Proportion of Time \u2265 Stage")
dev.off()



# ENP Stages --------------------------------------------------------------

dss_out=opendss(paste0(data.path,"Round2_RSMGL/",alts[3],"/RSMGL_output.dss"))
test=pathsToDataFrame(getCatalogedPathnames(dss_out), simplify=T)

RSM.stage.sites=c("ENP_NESRS1","ENP_NESRS2","ENP_NP-33","ENP_NP-35","ENP_NP-201","ENP_NP-203",
                  "ENP_NP-206","ENP_G3273","ENP_NP-EVER7","ENP_NP-EVER6","ENP_EVER4","G211U",
                  "S18CU","S177D","S178D")

subset(test,PARAMETER=="STAGE"&grepl("G211",LOCATION)==T)

subset(test,LOCATION%in%RSM.stage.sites)

stg.dat=data.frame()
for(j in 1:n.alts){
  dss_out=opendss(paste0(data.path,"Round2_RSMGL/",alts[j],"/RSMGL_output.dss"))  
  test=pathsToDataFrame(getCatalogedPathnames(dss_out), simplify=T)
  for(i in 1:length(RSM.stage.sites)){
    
    if(nrow(subset(test,LOCATION==RSM.stage.sites[i]))==1){
      paths=paste0("/RSMGL/",RSM.stage.sites[i],"/STAGE//1DAY/SIMULATED/")  
      tmp=data.frame(getFullTSC(dss_out,paths))
      tmp$Date=date.fun(rownames(tmp))
      rownames(tmp)<-NULL
      tmp$SITE=RSM.stage.sites[i]
      tmp$Alt=alts[j]
      stg.dat=rbind(tmp,stg.dat)
      print(i)
    }else{
      next
    }
  }
}

RSM.stage.sites.labs=c("NESRS1","NESRS2","NP-33","NP-35","NP-201","NP-203","NP-206","G3273","EVER7","EVER6","EVER4","G211 HW",
                       "S18C HW","S177 TW","S178 TW")

for(j in 1:length(RSM.stage.sites)){
png(filename=paste0(plot.path,"Round2/",RSM.stage.sites[j],"SDC.png"),width=7,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,0.75,0.5,1),oma=c(2.5,3.75,1,0.25));
layout(matrix(1:4,1,4,byrow=T))

  tmp=subset(stg.dat,SITE==RSM.stage.sites[j])
  ylim.val=c(floor(min(tmp$STAGE)),round(max(tmp$STAGE),2)+max(tmp$STAGE)*0.025)
  by.y=1;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
  
xlim.val=c(0,1);by.x=0.5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
for(i in 4:n.alts){
  # ylim.val=c(-0.15,9.5);by.y=1;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
  plot(0:1,0:1,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
  abline(h=ymaj,v=xmaj,lty=3,col="grey")
  with(ecdf_fun(subset(tmp,Alt==alts[1])$STAGE),lines(1-proportion,value,col=cols.alts[1],lty=1,lwd=1.5))
  with(ecdf_fun(subset(tmp,Alt==alts[2])$STAGE),lines(1-proportion,value,col=cols.alts[2],lty=2,lwd=1.5))
  with(ecdf_fun(subset(tmp,Alt==alts[3])$STAGE),lines(1-proportion,value,col=cols.alts[3],lty=1,lwd=1.5))
  if(nrow(subset(tmp,Alt==alts[i]))>0){
    with(ecdf_fun(subset(tmp,Alt==alts[i])$STAGE),lines(1-proportion,value,col=adjustcolor(cols.alts[i],0.5),lwd=2))
  }
  
  legend("topright",legend=c(alts[c(1:3,i)]),
         lty=c(1,2,1),lwd=c(1.5,1.5,1.5),col=c(cols.alts[1:3],adjustcolor(as.character(cols.alts[i]),0.5)),
         ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
  axis_fun(1,xmaj,xmin,format(xmaj))
  if(i==4){axis_fun(2,ymaj,ymin,ymaj)}else{axis_fun(2,ymaj,ymin,NA)}
  # if(i==4){mtext(side=3, adj=0,"SRS Total Annual Discharge")}
  if(i==4){mtext(side=3,adj=0,RSM.stage.sites.labs[j],cex=0.9)}
  box(lwd=1)
}
mtext(side=2,line=1.5,outer=T,"Stage (FT, NGVD29)")
mtext(side=1,line=1,outer=T,"Proportion of Time \u2265 Stage")
dev.off()
print(j)
}


# G211 --------------------------------------------------------------------



G211.q=subset(alt.srs,SITE%in%c("G211"))
G211.q$hydro.season=FL.Hydroseason(G211.q$Date)
G211.q$month=as.numeric(format(G211.q$Date,"%m"))
# wca3b.flow.xtab=dcast(wca3b.q,Alt+CY~SITE,value.var = "Q",sum)

G211.flow.ann=ddply(G211.q,c("Alt","CY","SITE"),summarise,TFlow.kAcftY=sum(Q/1000,na.rm=T))
G211.flow.ann=dcast(G211.flow.ann,Alt~SITE,value.var="TFlow.kAcftY",mean)

G211.flow.seasonal=ddply(G211.q,c("Alt","CY","hydro.season","SITE"),summarise,TFlow.kAcftY=sum(Q/1000,na.rm=T))
G211.flow.seasonal=dcast(G211.flow.seasonal,Alt+SITE~hydro.season,value.var="TFlow.kAcftY",mean)

G211.flow.month=ddply(G211.q,c("Alt","CY","month","SITE"),summarise,TFlow.kAcftY=sum(Q/1000,na.rm=T))
G211.flow.month=ddply(G211.flow.month,c("Alt","month","SITE"),summarise,mean.val=mean(TFlow.kAcftY,na.rm=T))

G211.flow.month.avg=ddply(G211.q,c("Alt","month","SITE"),summarise,mean.cfs=mean(FLOW,na.rm=T))

G211.HW.stg=subset(stg.dat,SITE=="G211U")

cols2=c("lightblue","khaki3")
# png(filename=paste0(plot.path,"Round2/G211_ann_month_sdf.png"),width=8,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.5,1),oma=c(3,2,0.75,0.25),lwd=0.5);
layout(matrix(1:3,1,3))

ylim.val=c(0,200);by.y=50;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
tmp=t(G211.flow.seasonal[,3:4])
x=barplot(t(G211.flow.seasonal[,3:4]),col=cols2,border="grey",
          space=0.05,
          ylim=ylim.val,axes=F,ann=F)
text(x,colSums(tmp),format(round(colSums(tmp),1)),offset=0.1,pos=3,cex=0.75)
text(x,tmp[1,]+tmp[2,]/2,format(round(tmp[2,],1)),font=3,cex=0.75)
text(x,tmp[1,]/2,format(round(tmp[1,],1)),font=3,cex=0.75)
axis_fun(1,x,x,alts,line=-0.5,cex=0.8,las=2)
axis_fun(2,ymaj,ymin,format(ymaj));
box(lwd=1)
mtext(side=3,adj=0," G-211")
mtext(side=2,outer=F,line=2.25,"Avg. Annual Discharge Volume (kAcFt Y\u207B\u00B9)",cex=0.8)
mtext(side=1,outer=F,line=2.6,"Alternatives")
legend("topright",legend=c("Wet (May - Oct)","Dry (Nov - April)"),
       lty=c(0),lwd=c(0.1),col="grey",pch=22,pt.bg=cols2,pt.cex=1.25,
       ncol=1,cex=0.75,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)


xlim.val=c(1,12);by.x=2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0,300);by.y=100;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
plot(0:1,0:1,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
# ylim.val=c(-0.15,9.5);by.y=1;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
with(subset(G211.flow.month.avg,Alt==alts[1]),lines(month,mean.cfs,col=cols.alts[1],lty=1,lwd=1.5))
with(subset(G211.flow.month.avg,Alt==alts[2]),lines(month,mean.cfs,col=cols.alts[2],lty=2,lwd=1.5))
with(subset(G211.flow.month.avg,Alt==alts[3]),lines(month,mean.cfs,col=cols.alts[3],lty=1,lwd=1.5))
for(i in 4:n.alts){
  with(subset(G211.flow.month.avg,Alt==alts[i]),lines(month,mean.cfs,col=cols.alts[i],lty=1,lwd=1.5))
}
axis_fun(1,xmaj,xmin,month.abb[xmaj])
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,outer=F,line=2.25,"Avg. Daily Discharge (ft\u00B3 sec\u207B\u00B9)",cex=0.8)
mtext(side=1,line=2.6,"Month")
legend("topleft",legend=c(alts),
       lty=c(1,2,1,1,1,1,1),lwd=c(1.5),col=c(cols.alts),
       ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)

xlim.val=c(0,1);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(3,8);by.y=1;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
# ylim.val=c(-0.15,9.5);by.y=1;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
plot(0:1,0:1,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
lty.val=c(1,2,1)
for(i in 1:3){
  tmp=ecdf_fun(subset(G211.HW.stg,Alt==alts[i])$STAGE)
  tmp$proportion=1-tmp$proportion
  # tmp$value=exp(tmp$value)
  with(tmp,lines(proportion,value,col=cols.alts[i],lty=lty.val[i],lwd=1.5))
}
for(i in 4:n.alts){
  with(ecdf_fun(subset(G211.HW.stg,Alt==alts[i])$STAGE),lines(1-proportion,value,col=adjustcolor(cols.alts[i],0.5),lwd=2))
}
mtext(side=3,adj=0," G-211 Headwater Stage")
legend("topright",legend=c(alts),
       lty=c(1,2,1,1,1,1,1),lwd=c(1.5),col=c(cols.alts),
       ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
axis_fun(1,xmaj,xmin,format(xmaj))
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,outer=F,line=2.5,"Headwater Stage (ft NGVD29)",cex=0.8)
mtext(side=1,line=2.6,"Proportion of Time \u2265 Stage")
dev.off()