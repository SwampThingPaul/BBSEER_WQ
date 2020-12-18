## STA WQBEL

###
library(openxlsx)
STAQBEL=read.xlsx("C:/Julian_LaCie/Work/Everglades_STAs/STA WQBEL Derivation.xlsx",sheet=3,startRow = 2)[,1:5]
colnames(STAQBEL)=c("STA", "WY", "GM_ugL", "FWM_ugL", "FWM.LT50")
STAQBEL=subset(STAQBEL,FWM.LT50=="Y")
STAQBEL.mean=ddply(STAQBEL,"STA",summarise,mean.FWM=mean(FWM_ugL),mean.GM=mean(GM_ugL))


#png(filename=paste0(plot.path,"STA_WQBEL.png"),width=4.25,height=4.75,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.75,0.25),oma=c(3,2,0.5,0.5));

xlim.val=c(0,55);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(FWM_ugL~GM_ugL,STAQBEL,xlim=xlim.val,ylim=xlim.val,ann=F,axes=F,type="n",xaxs="i",yaxs="i")
abline(h=xmaj,v=xmaj,lty=3,col="grey")
with(STAQBEL,points(GM_ugL,FWM_ugL,pch=21,bg="grey",lwd=0.1))
with(STAQBEL.mean,points(mean.GM,mean.FWM,pch=22,bg="dodgerblue1",lwd=0.1))
abline(0,1,lty=2)
abline(lm(FWM_ugL~GM_ugL+0,STAQBEL),col="indianred1")
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,xmaj,xmin,xmaj);box(lwd=1)
mtext(side=1,line=2,"TP AGM (\u03BCg L\u207B\u00B9)")
mtext(side=2,line=2,"TP FWM (\u03BCg L\u207B\u00B9)")
legend("bottomright",legend=c("Annual Value (STA)","Overall STA Mean Value","Regression Line","1:1 Line"),
       pch=c(21,22,NA,NA),
       lty=c(NA,NA,1,2),
       lwd=c(0.01,0.01,1,1),
       col=c("black","black","indianred1","black"),
       pt.bg=c("grey","dodgerblue1",NA,NA),
       pt.cex=1.5,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)  
dev.off()