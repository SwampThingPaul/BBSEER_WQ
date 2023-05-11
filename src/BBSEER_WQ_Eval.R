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

#Libraries
# devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);
library(plyr)
library(reshape2)

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

utm17=CRS("+init=epsg:26917")
wgs84=CRS("+init=epsg:4326")

## Function

DBHYDRO.meta.bysite=function(site,type,cat="SW"){
  # c("FLOW","STG","GATE")
  # cat = c("SW","RAIN","ETP")
  site.vals=paste0("v_site=",site)
  link=paste0("https://my.sfwmd.gov/dbhydroplsql/show_dbkey_info.show_dbkeys_matched?v_js_flag=Y&v_category=",cat,"&",
              site.vals,
              "&v_data_type=",type,
              "&v_dbkey_list_flag=Y&v_order_by=STATION")
  rslt.table=rvest::read_html(link)
  rslt.table=data.frame(rvest::html_table(rslt.table,fill=T)[[5]])
  rslt.table=rslt.table[,2:ncol(rslt.table)]
  colnames(rslt.table)=toupper(names(rslt.table))
  return(rslt.table)
  
}

DBHYDRO.meta.byDBKEY=function(DBKEY){
  link=paste0("https://my.sfwmd.gov/dbhydroplsql/show_dbkey_info.show_dbkeys_matched?v_js_flag=Y&v_dbkey=",
              paste(DBKEY,collapse="/"),
              "&v_category=SW")
  rslt.table=rvest::read_html(link)
  rslt.table=data.frame(rvest::html_table(rslt.table,fill=T)[[5]])
  rslt.table=rslt.table[,2:ncol(rslt.table)]
  colnames(rslt.table)=toupper(names(rslt.table))
  return(rslt.table)
  
}
leg.fun=function(b,pal,leg.title,
                 top.val=0.8,bot.val=0.2,mid.v.val=NULL,
                 x.max=0.3,x.min=0.1,mid.val=NULL,
                 txt.offset.val=-0.01,txt.y=NULL,leg.txt=NULL,
                 txt.cex=0.75,txt.adj=0,txt.pos=4,txt.offset=0.5,
                 title.cex=0.8,title.pos=3,title.adj=0,
                 title.x=NULL,title.y=NULL,
                 leg.type=c("continuous","categorical"), ...){
  l.b=length(b)
  labs=c(paste0("< ",b[2]),paste(b[2:(l.b-2)],b[3:(l.b-1)],sep=" - "),paste(paste0(">",b[(l.b-1)])))
  n.bks=length(b)-1
  mid.v.val=if(is.null(mid.v.val)==T){bot.val+(top.val-bot.val)/2}else{mid.v.val}
  
  mid.val=if(is.null(mid.val)==T){x.min+(x.max-x.min)/2}else{mid.val}
  if(leg.type=="continuous"){
    legend_image=as.raster(matrix(rev(pal),ncol=1))
    rasterImage(legend_image,x.min,bot.val,x.max,top.val)
    txt.y=if(is.null(txt.y)==T){c(bot.val,top.val)}else(txt.y)
    leg.txt=if(is.null(leg.txt)==T){format(c(min(b),max(b)))}else(leg.txt)
    text(x=x.max, y = txt.y, labels =leg.txt,cex=txt.cex,adj=txt.adj,pos=txt.pos,offset=txt.offset, ...)
  }
  if(leg.type=="categorical"){
    bx.val= seq(bot.val,top.val,(top.val-bot.val)/n.bks)
    rect(x.min,bx.val[1:n.bks],x.max,bx.val[2:(n.bks+1)],col=rev(pal),lty=0)
    text(y=bx.val[2:(n.bks+1)]-c(mean(diff(bx.val[2:(n.bks+1)]))/2), x = x.max, 
         labels = rev(labs),cex=txt.cex,xpd=NA,pos=txt.pos,adj=txt.adj)
  }
  
  title.x=if(is.null(title.x)==T){mid.val}else{title.x}
  title.y=if(is.null(title.y)==T){top.val}else{title.y}
  text(x=title.x,y=title.y,leg.title,adj=title.adj,cex=title.cex,pos=title.pos,xpd=NA)
}
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

NPS=spTransform(readOGR(paste0(gen.GIS,"/DOI/NPSBoundary"),"nps_boundary"),utm17)
BNP=subset(NPS,UNIT_CODE=="BISC")

FDEP_AP=spTransform(readOGR(paste0(gen.GIS,"/FDEP"),"AquaticPreserves"),utm17)

bbseer=spTransform(readOGR(GIS.path,"BBSEER_PRJBND_09092020"),utm17)

nnc=spTransform(readOGR(paste0(gen.GIS,"/FDEP"),"Estuary_NNC"),utm17)

# STORET/WIN locs ---------------------------------------------------------
storet.sites=read.table(paste0(data.path,"BBSEER_WQ_Eval/Station_Results.txt"),sep="|",header=T)
spl=sapply(strsplit(as.character(storet.sites$Latitude),split="\\s+"),as.numeric)
storet.sites$Latitude= spl[1,]+spl[2,]/60+spl[3,]/3600
spl=sapply(strsplit(as.character(storet.sites$Longitude),split="\\s+"),as.numeric)
storet.sites$Longitude= (spl[1,]+spl[2,]/60+spl[3,]/3600)*-1
storet.sites.shp=SpatialPointsDataFrame(storet.sites[,c("Longitude","Latitude")],data=storet.sites,proj4string=wgs84)
storet.sites.shp=spTransform(storet.sites.shp,utm17)

plot(storet.sites.shp)
plot(canals,add=T,col="blue")

win.sites=read.table(paste0(data.path,"BBSEER_WQ_Eval/_WIN_WAVES_JULIAN_P_4_20230502085216_64785.txt"),sep="|",header=T,skip=5)
head(win.sites)
unique(win.sites$Organization.ID)

win.sites.shp=SpatialPointsDataFrame(win.sites[,c("DEP.Longitude","DEP.Latitude")],data=win.sites,proj4string=wgs84)
win.sites.shp=spTransform(win.sites.shp,utm17)
plot(win.sites.shp,add=T,pch=21,bg="red")
plot(storet.sites.shp,add=T,col="white")

library(tmap)
tmap_mode("view")

tm_shape(storet.sites.shp)+tm_dots("red")+
  tm_shape(win.sites.shp)+tm_dots("blue",alpha=0.5)


head(storet.sites.shp@data)
head(win.sites.shp@data)

STORET.locs=cbind(data.frame(Station.ID=storet.sites.shp$Station.ID),coordinates(storet.sites.shp))
colnames(STORET.locs)=c("Station.ID","UTMX","UTMY")
WIN.locs=cbind(data.frame(Station.ID=win.sites.shp$Monitoring.Location.ID),coordinates(win.sites.shp))
colnames(WIN.locs)=c("Station.ID","UTMX","UTMY")

STORET.WIN.locs=rbind(STORET.locs,WIN.locs)
STORET.WIN.locs[duplicated(STORET.WIN.locs[,1])==F,]

plot(STORET.WIN.locs[duplicated(STORET.WIN.locs[,1])==F,c("UTMX","UTMY")])
STORET.WIN.locs=STORET.WIN.locs[duplicated(STORET.WIN.locs[,1])==F,]
STORET.WIN.locs.shp=SpatialPointsDataFrame(STORET.WIN.locs[,c("UTMX","UTMY")],
                                           data=STORET.WIN.locs,
                                           proj4string = utm17)
plot(STORET.WIN.locs.shp)
# -------------------------------------------------------------------------
dat.qual=data.frame(QUALIFIER=c(NA,"!","A","D","E","F","I","R","T","U","*","?","B","H","J","K","L","M","N","O","Q","V","Y","Z"),
                    FATALYN=c("N",rep("N",9),rep("Y",14)))


# dates=date.fun(c("2011-05-01","2022-04-30"))
dates=date.fun(c("2006-05-01","2022-04-30"))

# FDEP IWR database -------------------------------------------------------
# https://ryanpeek.org/2019-09-17-reading-databases-in-r/
#  library(odbc)
# # Microsoft Access Database Engine 2016 Redistributable
# # make sure this is install if first time running ODBC
# # https://www.microsoft.com/en-us/download/details.aspx?id=54920
# odbcListDrivers()
# file_path="C:/Julian_LaCie/_Data/FDEP/IWR/run64/rawDataDB4.accdb"
# accdb_con <- dbConnect(drv = odbc(),
#                        .connection_string = paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=",file_path,";"))
# dbListTables(accdb_con)
# RODBC::odbcConnectAccess2007(file_path)
# 
# dbSendQuery(accdb_con,"SELECT * FROM rawData4 WHERE year = 1994")
# 
# chan.path="C:/Julian_LaCie/_Data/FDEP/IWR/run64/iwr2002_run64_11292002.accdb"
# chan=RODBC::odbcConnectAccess2007(chan.path)
# tab.dat=RODBC::sqlTables(chan)
# subset(tab.dat,TABLE_TYPE=="TABLE")
# read.access(chan.path,"masterCodes_forPickList")
# 
# IWRparams=read.access(chan.path,"Master parameter codes")
# colnames(IWRparams)=gsub(" ","",names(IWRparams))
# unique(IWRparams$Category)
# subset(IWRparams,Category%in%c("Standard","Nutrient"))
# subset(IWRparams,Category%in%c("Nutrient"))
# 
# IWRparams[grepl("ortho",IWRparams$MasterparameterName),]
# 
# ##
# 
# iwr_1=read.access("C:/Julian_LaCie/_Data/FDEP/IWR/run64/rawDataDB1.accdb","rawData1")
# unique(iwr_1$sta)[grep('21FLDADE',unique(iwr_1$sta))]
# unique(iwr_1$sta)[grep('21FLBBAP',unique(iwr_1$sta))]
# range(iwr_1$year)
# 
# iwr_2=read.access("C:/Julian_LaCie/_Data/FDEP/IWR/run64/rawDataDB2.accdb","rawData2")
# unique(iwr_2$sta)[grep('21FLDADE',unique(iwr_2$sta))]
# unique(iwr_2$sta)[grep('21FLBBAP',unique(iwr_2$sta))]
# range(iwr_2$year)
# subset(iwr_2,year==1901)
# range(subset(iwr_2,year!=1901)$year)
# 
# iwr_3=read.access("C:/Julian_LaCie/_Data/FDEP/IWR/run64/rawDataDB3.accdb","rawData3")
# unique(iwr_3$sta)[grep('21FLDADE',unique(iwr_3$sta))]
# unique(iwr_3$sta)[grep('21FLBBAP',unique(iwr_3$sta))]
# range(iwr_3$year)
# 
# IWR.dat=rbind(
#   subset(iwr_1,grepl("21FLDADE",sta)==T|grepl("21FLBBAP",sta)==T),
#   subset(iwr_2,grepl("21FLDADE",sta)==T|grepl("21FLBBAP",sta)==T)
# )
# head(IWR.dat)
# IWR.dat$org=with(IWR.dat,substr(sta,1,8))
# IWR.dat$Station.ID=with(IWR.dat,substr(sta,9,35))
# 
# # unique(IWR.dat$rCode)
# # unique(IWR.dat$xCode)
# # test=IWR.dat[grepl(paste(subset(dat.qual,FATALYN=="Y")$QUALIFIER,collapse="|"),IWR.dat$xCode)==T,]
# # unique(test$rCode)
# # unique(test$xCode)
# 
# # QA/QC qualifiers 
# quals=as.character(unique(IWR.dat$xCode))
# spl=strsplit(quals,split="")
# quals=data.frame(xCode=quals,q1=sapply(spl,"[",1),q2=sapply(spl,"[",2),q3=sapply(spl,"[",3),q4=sapply(spl,"[",4))
# quals$Fatal=with(quals,ifelse(q1%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
#                                 q2%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
#                                 q3%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
#                                 q4%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER,"Y","N"))
# IWR.dat=merge(IWR.dat,quals[,c("xCode","Fatal")],"xCode",all.x=T)
# IWR.dat$date=with(IWR.dat,date.fun(paste(year,month,day,sep="-")))
# 
# unique(IWR.dat$masterCode)[order(unique(IWR.dat$masterCode))]
# ddply(IWR.dat,c("masterCode","param"),summarise,count=N.obs(result))
# 
# data.fram(masterCode=c("TP"))

# STORET/WIN --------------------------------------------------------------
win.dat1=read.table(paste0(data.path,"BBSEER_WQ_Eval/STORET_data/_WIN_WAVES_GUEST_20230502105337_64796__2017_2018.txt"),skip=9,sep="|",header=T)
win.dat2=read.table(paste0(data.path,"BBSEER_WQ_Eval/STORET_data/_WIN_WAVES_GUEST_20230502105337_64796__2018_2019.txt"),skip=9,sep="|",header=T)
win.dat3=read.table(paste0(data.path,"BBSEER_WQ_Eval/STORET_data/_WIN_WAVES_GUEST_20230502105337_64796__2019_2020.txt"),skip=9,sep="|",header=T)
win.dat4=read.table(paste0(data.path,"BBSEER_WQ_Eval/STORET_data/_WIN_WAVES_GUEST_20230502105337_64796__2020_2021.txt"),skip=9,sep="|",header=T)
win.dat5=read.table(paste0(data.path,"BBSEER_WQ_Eval/STORET_data/_WIN_WAVES_GUEST_20230502105337_64796__2021_2022.txt"),skip=9,sep="|",header=T)
win.dat=rbind(win.dat1,win.dat2,win.dat3,win.dat4,win.dat5)
rm(win.dat1,win.dat2,win.dat3,win.dat4,win.dat5)

win.dat$Activity.Start.Date.Time=date.fun(as.character(win.dat$Activity.Start.Date.Time),form="%m/%d/%Y %H:%M:%S")
win.dat$Date.EST=date.fun(win.dat$Activity.Start.Date.Time)
range(win.dat$Date.EST)
unique(win.dat$Monitoring.Location.ID)

win.dat$Station.ID=win.dat$Monitoring.Location.ID

# QA/QC qualifiers 
quals=as.character(unique(win.dat$Value.Qualifier))
spl=strsplit(quals,split="")
quals=data.frame(Value.Qualifier=quals,q1=sapply(spl,"[",1),q2=sapply(spl,"[",2),q3=sapply(spl,"[",3),q4=sapply(spl,"[",4))
quals$Fatal=with(quals,ifelse(q1%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q2%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q3%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q4%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER,"Y","N"))

# unique(win.dat$Matrix)
win.dat.clean=subset(win.dat,Matrix=="AQUEOUS-Surface Water")
win.dat.clean$HalfMDL=with(win.dat.clean,ifelse(Sample.Collection.Type=="Field Testing-Discrete",DEP.Result.Value.Number,
                                                ifelse(DEP.Result.Value.Number<=DEP.MDL,DEP.MDL/2,DEP.Result.Value.Number)))
win.dat.clean=merge(win.dat.clean,quals[,c("Value.Qualifier","Fatal")],"Value.Qualifier",all.x=T)
win.dat.clean=subset(win.dat.clean,Fatal=="N")

unique(win.dat.clean$DEP.Analyte.Name)
win.parameters=data.frame(DEP.Analyte.Name=c("Ammonia (N)", 
                                             "Chlorophyll a- uncorrected", 
                                             "Chlorophyll a, free of pheophytin","Chlorophyll a- corrected", "Dissolved Oxygen","Dissolved Oxygen Saturation", "Nitrate-Nitrite (N)", 
                                             "Nitrogen- Total Kjeldahl", "Orthophosphate (P)", "Phosphorus- Total", 
                                             "Specific Conductance", "Temperature, Water","pH","Salinity","Turbidity","Residues- Nonfilterable (TSS)","Color- Apparent","Carbon- Organic"),
                          param=c("NH4","Chla","Chla_c","Chla_c","DO","DOSat","NOx","TKN","SRP","TP","SPC","Temp","pH","Sal","Turb","TSS","Color","TOC"))
win.dat.clean=merge(win.dat.clean,win.parameters,"DEP.Analyte.Name")

win.xtab=dcast(win.dat.clean,Station.ID+Date.EST~param,value.var="HalfMDL",mean)
win.xtab$WY=WY(win.xtab$Date.EST)
win.xtab$TN=NA
win.xtab$TN=with(win.xtab,TN_Combine(NOx,TKN,TN))
win.xtab$DIN=with(win.xtab,NH4+NOx)

win.xtab$TPReversal=with(win.xtab,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>(TP*1.3),1,0)));# Reversals identified as 1 reversals consistent with TP rule evaluation
win.xtab$TNReversal=with(win.xtab,ifelse(is.na(DIN)==T|is.na(TN)==T,0,ifelse(DIN>(TN*1.3),1,0)));

# sum(win.xtab$TNReversal,na.rm=T)
# sum(win.xtab$TPReversal,na.rm=T)
# subset(win.xtab,TPReversal==1)

par(family="serif",oma=c(1,1,1,1),mar=c(4,4,1,1))
layout(matrix(1:2,1,2,byrow=F))
plot(TN~DIN,win.xtab,ylab="TN (mg L\u207B\u00b9)",xlab="DIN (mg L\u207B\u00b9)",pch=21,bg=ifelse(TNReversal==1,"dodgerblue1",NA),col=adjustcolor("grey",0.8));abline(0,1,col="dodgerblue1")
plot(TP~SRP,win.xtab,ylab="TP (mg L\u207B\u00b9)",xlab="SRP (mg L\u207B\u00b9)",pch=21,bg=ifelse(TPReversal==1,"red",NA),col=adjustcolor("grey",0.8));abline(0,1,col="red")

win.xtab$TP=with(win.xtab,ifelse(TPReversal==1,NA,TP))
win.xtab$SRP=with(win.xtab,ifelse(TPReversal==1,NA,SRP))
win.xtab$TN=with(win.xtab,ifelse(TNReversal==1,NA,TN))
win.xtab$DIN=with(win.xtab,ifelse(TNReversal==1,NA,DIN))


## STORET
storet.dat1=read.table(paste0(data.path,"BBSEER_WQ_Eval/STORET_data/Water_Quality_Results_2005_2008.txt"),sep="|",header=T)
storet.dat2=read.table(paste0(data.path,"BBSEER_WQ_Eval/STORET_data/Water_Quality_Results_2008_2011.txt"),sep="|",header=T)
storet.dat3=read.table(paste0(data.path,"BBSEER_WQ_Eval/STORET_data/Water_Quality_Results_2011_2013.txt"),sep="|",header=T)
storet.dat4=read.table(paste0(data.path,"BBSEER_WQ_Eval/STORET_data/Water_Quality_Results_2011_2013.txt"),sep="|",header=T)
storet.dat5=read.table(paste0(data.path,"BBSEER_WQ_Eval/STORET_data/Water_Quality_Results_2013_2015.txt"),sep="|",header=T)
storet.dat6=read.table(paste0(data.path,"BBSEER_WQ_Eval/STORET_data/Water_Quality_Results_2015_2017.txt"),sep="|",header=T)

storet.dat=rbind(storet.dat1,storet.dat2,storet.dat3,storet.dat4,storet.dat5,storet.dat6)
rm(storet.dat1,storet.dat2,storet.dat3)
storet.dat=storet.dat[!duplicated(storet.dat),]
storet.dat$Date.EST=date.fun(storet.dat$Act.Date,form="%m/%d/%Y")

range(storet.dat$Date.EST)


unique(storet.dat$Medium)
unique(storet.dat$Matrix)
unique(storet.dat$Act.Type)

# QA/QC qualifiers 
quals=as.character(unique(storet.dat$VQ))
spl=strsplit(quals,split="")
quals=data.frame(VQ=quals,q1=sapply(spl,"[",1),q2=sapply(spl,"[",2),q3=sapply(spl,"[",3),q4=sapply(spl,"[",4))
quals$Fatal=with(quals,ifelse(q1%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q2%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q3%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q4%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER,"Y","N"))
storet.dat$Result.Value=with(storet.dat,ifelse(Result.Value=="*Non-detect",MDL*-1,as.numeric(as.character(Result.Value))))

storet.dat$HalfMDL=with(storet.dat,ifelse(Act.Type=="Field Msr/Obs",Result.Value,
                                          ifelse(abs(Result.Value)<=MDL,MDL/2,abs(Result.Value))))
storet.dat=merge(storet.dat,quals[,c("VQ","Fatal")],"VQ",all.x=T)
storet.dat=subset(storet.dat,Fatal=="N")

ddply(storet.dat,c("Characteristic","Result.Units"),summarise,n.val=N.obs(Station.ID))

storet.parameters=data.frame(Characteristic=c("Turbidity", "Total Suspended Solids (TSS)", "Nitrogen, ammonia (NH3) as NH3", 
                                              "Phosphorus as PO4", "Nitrogen, Nitrite (NO2) + Nitrate (NO3) as N", 
                                              "Phosphorus, phosphate (PO4) as PO4","Chlorophyll a, free of pheophytin","Chlorophyll a, corrected for pheophytin", 
                                              "Apparent Color","Nitrogen, Kjeldahl","Phosphorus, orthophosphate as P","Phosphorus, orthophosphate as PO4","pH", "Dissolved oxygen saturation", 
                                              "Dissolved oxygen (DO)", "Salinity", "Specific conductance","Temperature, water", "Nitrogen, ammonia as N","Total Organic Carbon (TOC)"),
                             param=c("Turb","TSS","NH4","TP","NOx","TP","Chla_c","Chla_c","Color","TKN","SRP","SRP","pH","DOSat","DO","Sal","SPC","Temp","NH4","TOC"))

storet.dat=merge(storet.dat,storet.parameters,"Characteristic")

storet.xtab=dcast(storet.dat,Station.ID+Date.EST~param,value.var="HalfMDL",mean)
storet.xtab$WY=WY(storet.xtab$Date.EST)
storet.xtab$TN=NA
storet.xtab$TN=with(storet.xtab,TN_Combine(NOx,TKN,TN))
storet.xtab$DIN=with(storet.xtab,NH4+NOx)

storet.xtab$TPReversal=with(storet.xtab,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>(TP*1.3),1,0)));# Reversals identified as 1 reversals consistent with TP rule evaluation
storet.xtab$TNReversal=with(storet.xtab,ifelse(is.na(DIN)==T|is.na(TN)==T,0,ifelse(DIN>(TN*1.3),1,0)));

# sum(storet.xtab$TNReversal,na.rm=T)
# sum(storet.xtab$TPReversal,na.rm=T)
# subset(storet.xtab,TPReversal==1)

par(family="serif",oma=c(1,1,1,1),mar=c(4,4,1,1))
layout(matrix(1:2,1,2,byrow=F))
plot(TN~DIN,storet.xtab,ylab="TN (mg L\u207B\u00b9)",xlab="DIN (mg L\u207B\u00b9)",pch=21,bg=ifelse(TNReversal==1,"dodgerblue1",NA),col=adjustcolor("grey",0.8));abline(0,1,col="dodgerblue1")
plot(TP~SRP,storet.xtab,ylab="TP (mg L\u207B\u00b9)",xlab="SRP (mg L\u207B\u00b9)",pch=21,bg=ifelse(TPReversal==1,"red",NA),col=adjustcolor("grey",0.8));abline(0,1,col="red")

storet.xtab$TP=with(storet.xtab,ifelse(TPReversal==1,NA,TP))
storet.xtab$SRP=with(storet.xtab,ifelse(TPReversal==1,NA,SRP))
storet.xtab$TN=with(storet.xtab,ifelse(TNReversal==1,NA,TN))
storet.xtab$DIN=with(storet.xtab,ifelse(TNReversal==1,NA,DIN))

head(storet.xtab)
head(win.xtab)
STORET_WIN=rbind(storet.xtab[,names(win.xtab)],
                 win.xtab)
range(STORET_WIN$Date.EST)
plot(TP~Date.EST,subset(STORET_WIN,Station.ID=="AC01"))
plot(result~date,subset(IWR.dat,Station.ID=="AC01"&masterCode=="TP"&year>=2011))

plot(TP~Date.EST,subset(STORET_WIN,Station.ID=="AC01"),ylim=c(0,0.05))
# points(result~date,subset(IWR.dat,Station.ID=="AC01"&masterCode=="TP"&year>=2011),pch=21,bg="red")

range(STORET_WIN$Date.EST)
STORET_WIN=subset(STORET_WIN,Date.EST%in%seq(dates[1],dates[2],"1 days"))

# SFWMD Data --------------------------------------------------------------
## S332B flow from PREF key check metadata
## S332C WQ eliminated in 2020, uses S332B as surrogate
struc.sites=data.frame(
  STRUCT=c("S173","S331","S332B","S332C","S332D","S176","S200","G737","S177","S199","S178",'S18C',"S197",
           "S20F","S20G","S21A","S21","S123","S700","S22"),
  DBKEY=c(91419,91479,"TB064","91483","91485","91422","91437","AN674","91423","91436","91424","91427","91435",
          "91441","91442","91445","91447","91367","91665","91449"),
  WQ.site=c("S331-173","S331","S332B","S332B","S332DX","S332DX","S200","G737","S177","S177","S178","S18C","S197",
            "S20F","S20G","PR01","BL02","CD02","S700","SP01")
)

struc.sites$WQ.site=with(struc.sites,ifelse(STRUCT=="S20F","MW01",WQ.site))
struc.sites$WQ.site=with(struc.sites,ifelse(STRUCT=="S20G","MI01",WQ.site))
subset(struc.sites,STRUCT=="S20F")

struc.sites$SITE=struc.sites$STRUCT

other.WQ=data.frame(
  WQ.site=c("AJC1","MW01","MI01","GL03","EP","BB53","BB52",
            paste0("FLAB",c("07","08","10","11","12","06","09","23","13","04","05","24","21"))),
  STRUCT=c("S199","S20F","S20G",NA,NA,NA,NA,
              rep(NA,13))
)
other.WQ$SITE=other.WQ$WQ.site

## MW01, MI01, CD02/CD01A, SP01 current data not in DBHYDRO
## FLAB07 flow site is HIGHWAY_CR - DBKEY 39203
## FLAB10 flow site is TROUTCR_B - DBKEY AN866

## Flow data
DBHYDRO.meta.byDBKEY(struc.sites$DBKEY)

flow.dat=data.frame()
for(i in 1:nrow(struc.sites)){
  tmp=DBHYDRO_daily(dates[1],dates[2],struc.sites$DBKEY[i])
  tmp$DKEY=as.character(struc.sites$DBKEY[i])
  flow.dat=rbind(tmp,flow.dat)
  print(i)
}

flow.dat$Date.EST=date.fun(flow.dat$Date)
flow.dat=merge(flow.dat,struc.sites,"DBKEY")

plot(Data.Value~Date.EST,subset(flow.dat,STRUCT=="S173"))
plot(Data.Value~Date.EST,subset(flow.dat,STRUCT=="S331"))

flow.dat$data.value.14=with(flow.dat,ave(Data.Value,STRUCT,FUN=function(x) zoo::rollapply(x>0,14,sum,na.rm=T,align="right",fill=NA)))
flow.dat$Q.14dmean=with(flow.dat,ave(Data.Value,STRUCT,FUN=function(x) zoo::rollapply(x,14,mean,na.rm=T,align="right",fill=NA)))

## WQ Data
params=data.frame(Test.Number=c(18,21,80,20,25,23,61,179,7,16,8,9,98,12,13,10),
                  param=c("NOx","TKN","TN","NH4","TP","SRP","Chla","Chla","Temp","TSS","DO","SPC","Sal","Turb","Color","pH"))

WQ.sites=c(unique(struc.sites$WQ.site),unique(other.WQ$WQ.site))
WQ.sites=WQ.sites[!(WQ.sites%in%c("MW01", "MI01", "CD02","CD01A", "SP01"))]

wq.dat=data.frame()
for(i in 1:length(WQ.sites)){
  tmp=DBHYDRO_WQ(dates[1],dates[2],WQ.sites[i],params$Test.Number)
  wq.dat=rbind(wq.dat,tmp)
  print(i)
}

wmd.dat=merge(wq.dat,params,"Test.Number")
unique(wmd.dat$Collection.Method)
wmd.dat=subset(wmd.dat,Collection.Method=="G")

wmd.dat.xtab=dcast(wmd.dat,Station.ID+Date.EST~param,value.var="HalfMDL",mean)
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

wmd.dat.xtab$TP=with(wmd.dat.xtab,ifelse(TPReversal==1,NA,TP))
wmd.dat.xtab$SRP=with(wmd.dat.xtab,ifelse(TPReversal==1,NA,SRP))
wmd.dat.xtab$TN=with(wmd.dat.xtab,ifelse(TNReversal==1,NA,TN))
wmd.dat.xtab$DIN=with(wmd.dat.xtab,ifelse(TNReversal==1,NA,DIN))
wmd.dat.xtab$Sal=with(wmd.dat.xtab,SalinityCalc(SPC,Temp))

wmd.dat.xtab=wmd.dat.xtab[,names(STORET_WIN)]
unique(STORET_WIN$Station.ID)[unique(STORET_WIN$Station.ID)%in%unique(wmd.dat.xtab$Station.ID)]

wmd.dat.xtab$source="DBHYDRO"

## data from SFWMD to fill early 2000s gap
N.mw=14.0067
P.mw=30.973762
C.mw=12.0107

data.path.trend="C:/Julian_LaCie/_GitHub/EVER_FKNMS_WQTrend/Data/"
wmd.dat2=read.csv(paste0(data.path.trend,"SFWMD/CWQMN_thru_Apr2021_For_Analysis.csv"))
names(wmd.dat2)
wmd.names=c("DATE", "month", "Yr", "WY.MON", "WY", "Season", 
            "Station.Num", "Station.ID", "Station.ID.2", "Region", "ZSI", 
            "Temp", "spc_uScm", "Sal", "pH", "DO.mgl", 
            "DO_persat", "Turb", "Chla", "TN.mgL", "NOx.uM", 
            "NO2.uM", "NH4.uM", "TN.uM", "DIN.uM", "TON.uM", 
            "TP.uM", "SRP.uM", "TOC.uM")
colnames(wmd.dat2)=wmd.names
head(wmd.dat2)
summary(wmd.dat2)
wmd.names=c("DATE", "month", "Yr", "WY.MON", "WY", "Season", 
            "Station.Num", "Station.ID", "Station.ID.2", "Region", "ZSI", 
            "Temp", "spc_uScm", "Sal", "pH", "DO.mgl", 
            "DO_persat", "Turb", "Chla", "TN.mgL", "NOx.uM", 
            "NO2.uM", "NH4.uM", "TN.uM", "DIN.uM", "TON.uM", 
            "TP.uM", "SRP.uM", "TOC.uM")
colnames(wmd.dat2)=wmd.names
head(wmd.dat2)
summary(wmd.dat2)

wmd.dat2$Station.ID
wmd.dat2=subset(wmd.dat2,Station.ID%in%paste0("FLAB",c("07","08","10","11","12","06","09","23","13","04","05","24","21")))

subset(wmd.dat2,Turb<=0)
subset(wmd.dat2,Chla<=0)
subset(wmd.dat2,NOx.uM<=0)
subset(wmd.dat2,NH4.uM<=0)
subset(wmd.dat2,TOC.uM<=0)
subset(wmd.dat2,SRP.uM<=0)
subset(wmd.dat2,TP.uM<=0)

wmd.dat2$Turb[wmd.dat2$Turb<=0]=NA
wmd.dat2$Chla[wmd.dat2$Chla<=0]=NA
wmd.dat2$NOx.uM[wmd.dat2$NOx.uM==0]=NA
wmd.dat2$NH4.uM[wmd.dat2$NH4.uM==0]=NA
wmd.dat2$TOC.uM[wmd.dat2$TOC.uM<=0]=NA
wmd.dat2$TP.uM[wmd.dat2$TP.uM==0]=NA
wmd.dat2$SRP.uM[wmd.dat2$SRP.uM==0]=NA

wmd.dat2$DATE=date.fun(wmd.dat2$DATE,form="%m/%d/%Y")
wmd.dat2$Date.EST=wmd.dat2$DATE#date.fun(wmd.dat2$DATE,form="%m/%d/%Y")
wmd.dat2$DO=wmd.dat2$DO.mgl
wmd.dat2$NH4=with(wmd.dat2,ifelse(NH4.uM*N.mw/1000<0.0008,0.0008/2,round(NH4.uM*N.mw/1000,4)))
wmd.dat2$NOx=with(wmd.dat2,ifelse(NOx.uM*N.mw/1000<0.0098,0.0098/2,round(NOx.uM*N.mw/1000,4)))
wmd.dat2$DIN=with(wmd.dat2,NOx+NH4)
wmd.dat2$TKN=NA
wmd.dat2$TN=with(wmd.dat2,ifelse(TN.uM*N.mw/1000<0.02,0.02/2,round(TN.uM*N.mw/1000,4)))
wmd.dat2$TOC=with(wmd.dat2,ifelse(TOC.uM*C.mw/1000<0.05,0.05/2,round(TOC.uM*C.mw/1000,4)))
wmd.dat2$TP=with(wmd.dat2,ifelse(TP.uM*P.mw/1000<0.002,0.002/2,round(TP.uM*P.mw/1000,4)))
wmd.dat2$SRP=with(wmd.dat2,ifelse(SRP.uM*P.mw/1000<0.002,0.002/2,round(SRP.uM*P.mw/1000,4)))
wmd.dat2$DOC=NA
wmd.dat2$WY=WY(wmd.dat2$Date.EST)
wmd.dat2$source="CWQMN"
wmd.dat2$SALINITY=with(wmd.dat2,SalinityCalc(spc_uScm,Temp))
wmd.dat2$SPC=wmd.dat2$spc_uScm

wmd.dat2$season=FL.Hydroseason(wmd.dat2$Date.EST)
wmd.dat2$STATION=wmd.dat2$Station.ID

wmd.dat2$TPReversal=with(wmd.dat2,ifelse(is.na(SRP)==T|is.na(TP)==T,0,ifelse(SRP>(TP*1.3),1,0)));# Reversals identified as 1 reversals consistent with TP rule evaluation
wmd.dat2$TNReversal=with(wmd.dat2,ifelse(is.na(DIN)==T|is.na(TN)==T,0,ifelse(DIN>(TN*1.3),1,0)));

par(family="serif",oma=c(1,1,1,1),mar=c(4,4,1,1))
layout(matrix(1:2,1,2,byrow=F))
plot(TN~DIN,wmd.dat2,ylab="TN (mg L\u207B\u00b9)",xlab="DIN (mg L\u207B\u00b9)",pch=21,bg=ifelse(TNReversal==1,"dodgerblue1",NA),col=adjustcolor("grey",0.8));abline(0,1,col="dodgerblue1")
plot(TP~SRP,wmd.dat2,ylab="TP (mg L\u207B\u00b9)",xlab="SRP (mg L\u207B\u00b9)",pch=21,bg=ifelse(TPReversal==1,"red",NA),col=adjustcolor("grey",0.8));abline(0,1,col="red")

wmd.dat2$TP=with(wmd.dat2,ifelse(TPReversal==1,NA,TP))
wmd.dat2$SRP=with(wmd.dat2,ifelse(TPReversal==1,NA,SRP))
wmd.dat2$TN=with(wmd.dat2,ifelse(TNReversal==1,NA,TN))
wmd.dat2$DIN=with(wmd.dat2,ifelse(TNReversal==1,NA,DIN))

wmd.dat2$source="SFWMD_ERDB"
vars=c("Station.ID", "Date.EST", "DO", "NH4", "NOx", "pH", "SPC", 
  "SRP", "Temp", "TKN", "TP", "WY", "TN", "DIN", "TPReversal", 
  "TNReversal", "source")
names(wmd.dat2)
match(vars,names(wmd.dat2))
wmd.dat2=wmd.dat2[,vars]

dat1=subset(wmd.dat.xtab,Station.ID%in%paste0("FLAB",c("07","08","10","11","12","06","09","23","13","04","05","24","21")))

vars=c("Station.ID","Date.EST")
sum(duplicated(rbind(wmd.dat2[,vars],dat1[,vars])))

wmd.dat2=rbind(wmd.dat2,dat1)[duplicated(rbind(wmd.dat2[,vars],dat1[,vars]))==F,]
layout(matrix(1:15,3,5))

sites=paste0("FLAB",c("07","08","10","11","12","06","09","23","13","04","05","24","21"))
for(i in 1:length(sites)){
plot(TN~Date.EST,subset(wmd.dat2,Station.ID=="FLAB05"))
points(TN~Date.EST,subset(dat1,Station.ID=="FLAB05"),pch=21,bg="red")
points(TN~Date.EST,subset(test,Station.ID=="FLAB05"),pch=21,bg="green")
}
# WQ Data Combined --------------------------------------------------------
STORET_WIN$source="STORET"
FLAB.sites=paste0("FLAB",c("07","08","10","11","12","06","09","23","13","04","05","24","21"))

WQ.dat.all=rbind(subset(wmd.dat.xtab,!(Station.ID%in%FLAB.sites)),
                 wmd.dat2,
                 STORET_WIN)
tmp=ddply(WQ.dat.all,"Station.ID",summarise,min.date=min(Date.EST),max.date=max(Date.EST))
subset(tmp,substr(Station.ID,1,2)=="BI")
subset(tmp,substr(Station.ID,1,2)=="L3")
subset(tmp,substr(Station.ID,1,2)=="S2")

# write.csv(WQ.dat.all,paste0(export.path,"/20230504_BBSEERWQ_data.csv"),row.names = F)
idvars=c("Station.ID","WY","Date.EST")
vars=c("TP","SRP","TN","DIN")
WQ.dat.all=melt(WQ.dat.all[,c(idvars,vars)],id.vars=idvars)

WQ.dat.all$season=FL.Hydroseason(WQ.dat.all$Date.EST)
# WQ.dat.all=ddply(WQ.dat.all,c("Station.ID","Date.EST","WY","season","variable"),summarise,value=mean(value,na.rm=T),N.val=N.obs(value))
# subset(WQ.dat.all,N.val>1)

season.screen=dcast(WQ.dat.all,Station.ID+variable+WY~season,value.var="value",N.obs)
season.screen$totalN=rowSums(season.screen[,c("A_Wet","B_Dry")],na.rm=T)
season.screen$screen=with(season.screen, ifelse(A_Wet>0&B_Dry>0&totalN>=4,1,0))

WQ.dat.all=merge(WQ.dat.all,season.screen[,c("Station.ID","WY","variable","screen")],c("Station.ID","WY","variable"))

head(flow.dat)
q.vars=c("Date.EST","WQ.site","SITE","Data.Value","data.value.14","Q.14dmean")
head(flow.dat[,q.vars])
struct.wq=merge(WQ.dat.all,flow.dat[,q.vars],by.x=c("Station.ID","Date.EST"),by.y=c("WQ.site","Date.EST"),all.x=T)
struct.wq$SITE2=with(struct.wq,ifelse(is.na(SITE)==T,Station.ID,SITE))
struct.wq$type=with(struct.wq,ifelse(is.na(SITE)==T,"ambient","struct"))


subset(struct.wq,Station.ID=="S178"&variable=="TP")

## Structure Data ----------------------------------------------------------
## Annual means
WY_FWM_GM.struct=ddply(subset(struct.wq,screen==1&type=="struct"),c("SITE2","WY","variable"),summarise,
           N.wq=N.obs(value),
           N.Q=N.obs(ifelse(is.na(Data.Value)==T|Data.Value==0,NA,value)),
           FWM=Hmisc::wtd.mean(value,Data.Value,na.rm=T),
           GM=exp(mean(log(ifelse(is.na(Data.Value)==T|Data.Value==0,NA,value)),na.rm=T)))
unique(WY_FWM_GM.struct$SITE2)

## Trends
N.trend.screen=ddply(WY_FWM_GM.struct,c("SITE2","variable"),summarise,N.val=N.obs(GM))
N.trend.screen$trend.screen=with(N.trend.screen,ifelse(N.val>3,1,0))
WY_FWM_GM.struct=merge(WY_FWM_GM.struct,N.trend.screen[,c("SITE2","variable","trend.screen")],c("SITE2","variable"),all.x=T)
WY_GM_trend.struct=ddply(subset(WY_FWM_GM.struct,trend.screen==1),c("SITE2","variable"),summarise,
      tau=as.numeric(cor.test(GM,WY,method="kendall")$estimate),
      pval=as.numeric(cor.test(GM,WY,method="kendall")$p.value),
      sen.slope=as.numeric(zyp::zyp.sen(GM~WY)$coefficients[2]),
      N.WY=N.obs(WY))
subset(WY_GM_trend.struct,pval<0.05)

##
wmd.struct2=subset(wmd.struct,NAME%in%unique(WY_FWM_GM.struct$SITE2))
wmd.struct2=cbind(data.frame(SITE2=wmd.struct2$NAME),coordinates(wmd.struct2))
colnames(wmd.struct2)=c("SITE2","UTMX","UTMY")

## Mean WY TP SHP
mean.WY_FWM_GM.struct=ddply(WY_FWM_GM.struct,c("SITE2","variable"),summarise,
                     mean.FWM=mean(FWM,na.rm=T),
                     mean.GM=mean(GM,na.rm=T),
                     sd.GM=sd(GM,na.rm=T),
                     N.yrs=N.obs(GM))
mean.WY_FWM_GM.struct$DOF=with(mean.WY_FWM_GM.struct,N.yrs-1)
mean.WY_FWM_GM.struct$Tp=with(mean.WY_FWM_GM.struct,abs(qt(0.05,DOF)))
mean.WY_FWM_GM.struct$UCI=with(mean.WY_FWM_GM.struct, mean.GM+(sd.GM*Tp)/sqrt(N.yrs))

mean.WY_FWM_GM.struct=merge(mean.WY_FWM_GM.struct,N.trend.screen[,c("SITE2","variable","trend.screen")],c("SITE2","variable"),all.x=T)
mean.WY_FWM_GM.struct=subset(mean.WY_FWM_GM.struct,trend.screen==1)

mean.WY_FWM_GM.struct.TP=merge(subset(mean.WY_FWM_GM.struct,variable=="TP"),
                               wmd.struct2,"SITE2")
mean.WY_FWM_GM.struct.TP$mean.GM=mean.WY_FWM_GM.struct.TP$mean.GM*1000
mean.WY_FWM_GM.struct.TP$UCI=mean.WY_FWM_GM.struct.TP$UCI*1000
# mean.WY_FWM_GM.struct.TP.shp=SpatialPointsDataFrame(mean.WY_FWM_GM.struct.TP[,c("UTMX","UTMY")],
#                                                     data=mean.WY_FWM_GM.struct.TP,proj4string=utm17)
# ## Trend TP SHP
WY_GM_trend.struct.TP=merge(subset(WY_GM_trend.struct,variable=="TP"),
                               wmd.struct2,"SITE2")
# WY_GM_trend.struct.TP.shp=SpatialPointsDataFrame(WY_GM_trend.struct.TP[,c("UTMX","UTMY")],
#                                                     data=WY_GM_trend.struct.TP,proj4string=utm17)

## Mean WY TN SHP
mean.WY_FWM_GM.struct.TN=merge(subset(mean.WY_FWM_GM.struct,variable=="TN"),
                               wmd.struct2,"SITE2")
# mean.WY_FWM_GM.struct.TN.shp=SpatialPointsDataFrame(mean.WY_FWM_GM.struct.TN[,c("UTMX","UTMY")],
#                                                     data=mean.WY_FWM_GM.struct.TN,proj4string=utm17)

## Trend TN SHP
WY_GM_trend.struct.TN=merge(subset(WY_GM_trend.struct,variable=="TN"),
                            wmd.struct2,"SITE2")
# WY_GM_trend.struct.TN.shp=SpatialPointsDataFrame(WY_GM_trend.struct.TN[,c("UTMX","UTMY")],
#                                              data=WY_GM_trend.struct.TN,proj4string=utm17)



## "Ambient Data" ------------------------------------------------------------
## Annual means
WY_FWM_GM.ambient=ddply(subset(struct.wq,screen==1&type=="ambient"),c("SITE2","WY","variable"),summarise,
                        N.wq=N.obs(value),
                        N.Q=NA,
                        FWM=NA,
                        GM=exp(mean(log(value),na.rm=T)))
mean.WY_FWM_GM.ambient=ddply(WY_FWM_GM.ambient,c("SITE2","variable"),summarise,
                            mean.FWM=mean(FWM,na.rm=T),
                            mean.GM=mean(GM,na.rm=T),
                            sd.GM=sd(GM,na.rm=T),
                            N.yrs=N.obs(GM))
mean.WY_FWM_GM.ambient$DOF=with(mean.WY_FWM_GM.ambient,N.yrs-1)
mean.WY_FWM_GM.ambient$Tp=with(mean.WY_FWM_GM.ambient,abs(qt(0.05,DOF)))
mean.WY_FWM_GM.ambient$UCI=with(mean.WY_FWM_GM.ambient, mean.GM+(sd.GM*Tp)/sqrt(N.yrs))

unique(WY_FWM_GM.ambient$SITE2)

# plot(GM~WY,subset(WY_FWM_GM.ambient,SITE2=='FLAB04'))

## Trends
N.trend.screen=ddply(WY_FWM_GM.ambient,c("SITE2","variable"),summarise,N.val=N.obs(GM))
N.trend.screen$trend.screen=with(N.trend.screen,ifelse(N.val>3,1,0))
WY_FWM_GM.ambient=merge(WY_FWM_GM.ambient,N.trend.screen[,c("SITE2","variable","trend.screen")],c("SITE2","variable"),all.x=T)
WY_GM_trend.ambient=ddply(subset(WY_FWM_GM.ambient,trend.screen==1),c("SITE2","variable"),summarise,
                         tau=as.numeric(cor.test(GM,WY,method="kendall")$estimate),
                         pval=as.numeric(cor.test(GM,WY,method="kendall")$p.value),
                         sen.slope=as.numeric(zyp::zyp.sen(GM~WY)$coefficients[2]),
                         N.WY=N.obs(WY))
subset(WY_GM_trend.ambient,pval<0.05)
plot(GM~WY,subset(WY_FWM_GM.ambient,SITE2=="FLAB04"&variable=="TP"))

mean.WY_FWM_GM.ambient=merge(mean.WY_FWM_GM.ambient,N.trend.screen[,c("SITE2","variable","trend.screen")],c("SITE2","variable"),all.x=T)
mean.WY_FWM_GM.ambient=subset(mean.WY_FWM_GM.ambient,trend.screen==1)
## 
win.sites$Station.ID=win.sites$Monitoring.Location.ID
win.sites$Latitude=win.sites$DEP.Latitude
win.sites$Longitude=win.sites$DEP.Longitude

STORET.WIN.sites=rbind(
  win.sites[,c("Station.ID","Latitude","Longitude")],
  storet.sites[,c("Station.ID","Latitude","Longitude")]
)
sum(duplicated(STORET.WIN.sites$Station.ID))

STORET.WIN.sites=STORET.WIN.sites[duplicated(STORET.WIN.sites$Station.ID)==F,]
STORET.WIN.sites=SpatialPointsDataFrame(STORET.WIN.sites[,c("Longitude","Latitude")],
                                        data=STORET.WIN.sites,proj4string=wgs84)
STORET.WIN.sites=spTransform(STORET.WIN.sites,utm17)

## Mean WY TP SHP
wmd.sites=c("EP",paste0("FLAB",c("07","08","10","11","12","06","09","23","13","04","05","24","21")),"AJC1")
WY_FWM_GM.ambient.TP1=subset(mean.WY_FWM_GM.ambient,!(SITE2%in%wmd.sites)&variable=="TP")
WY_FWM_GM.ambient.TP2=subset(mean.WY_FWM_GM.ambient,SITE2%in%wmd.sites&variable=="TP")

STORET.WIN.sites1=cbind(data.frame(SITE2=subset(STORET.WIN.sites,!(Station.ID%in%wmd.sites))$Station.ID),
                        coordinates(subset(STORET.WIN.sites,!(Station.ID%in%wmd.sites))))
colnames(STORET.WIN.sites1)=c("SITE2","UTMX","UTMY")
WY_FWM_GM.ambient.TP1=merge(WY_FWM_GM.ambient.TP1,STORET.WIN.sites1,"SITE2")

wmd.mon1=cbind(data.frame(SITE2=subset(sfwmd.mon,ACTIVITY_S=="Surface Water Grab"&STATION%in%wmd.sites)$STATION),
               coordinates(subset(sfwmd.mon,ACTIVITY_S=="Surface Water Grab"&STATION%in%wmd.sites)))
colnames(wmd.mon1)=c("SITE2","UTMX","UTMY")
WY_FWM_GM.ambient.TP2=merge(WY_FWM_GM.ambient.TP2,wmd.mon1,"SITE2")

WY_FWM_GM.ambient.TP=rbind(WY_FWM_GM.ambient.TP1,WY_FWM_GM.ambient.TP2)
WY_FWM_GM.ambient.TP$mean.GM=WY_FWM_GM.ambient.TP$mean.GM*1000
WY_FWM_GM.ambient.TP$UCI=WY_FWM_GM.ambient.TP$UCI*1000
## All data AGM TP Combined
mean.WY_FWM_GM.TP.all=rbind(mean.WY_FWM_GM.struct.TP,WY_FWM_GM.ambient.TP)
mean.WY_FWM_GM.TP.all.shp=SpatialPointsDataFrame(mean.WY_FWM_GM.TP.all[,c("UTMX","UTMY")],
                                                    data=mean.WY_FWM_GM.TP.all,proj4string=utm17)

## Trend TP SHP
WY_GM_trend.ambient.TP1=merge(subset(WY_GM_trend.ambient,!(SITE2%in%wmd.sites)&variable=="TP"),
                            STORET.WIN.sites1,"SITE2")
WY_GM_trend.ambient.TP2=merge(subset(WY_GM_trend.ambient,(SITE2%in%wmd.sites)&variable=="TP"),
                            wmd.mon1,"SITE2")
WY_GM_trend.ambient.TP=rbind(WY_GM_trend.ambient.TP1,WY_GM_trend.ambient.TP2)

## All data trend TP Combined
WY_GM_trend.TP.all=rbind(WY_GM_trend.struct.TP,WY_GM_trend.ambient.TP)
WY_GM_trend.TP.all.shp=SpatialPointsDataFrame(WY_GM_trend.TP.all[,c("UTMX","UTMY")],
                                                 data=WY_GM_trend.TP.all,proj4string=utm17)

## Mean WY TN SHP
WY_FWM_GM.ambient.TN1=merge(subset(mean.WY_FWM_GM.ambient,!(SITE2%in%wmd.sites)&variable=="TN"),
                            STORET.WIN.sites1,"SITE2")
WY_FWM_GM.ambient.TN2=merge(subset(mean.WY_FWM_GM.ambient,(SITE2%in%wmd.sites)&variable=="TN"),
                            wmd.mon1,"SITE2")
WY_FWM_GM.ambient.TN=rbind(WY_FWM_GM.ambient.TN1,WY_FWM_GM.ambient.TN2)

## All data AGM TN Combined
mean.WY_FWM_GM.TN.all=rbind(mean.WY_FWM_GM.struct.TN,WY_FWM_GM.ambient.TN)
mean.WY_FWM_GM.TN.all.shp=SpatialPointsDataFrame(mean.WY_FWM_GM.TN.all[,c("UTMX","UTMY")],
                                                 data=mean.WY_FWM_GM.TN.all,proj4string=utm17)


## Trend TN SHP
WY_GM_trend.ambient.TN1=merge(subset(WY_GM_trend.ambient,!(SITE2%in%wmd.sites)&variable=="TN"),
                              STORET.WIN.sites1,"SITE2")
WY_GM_trend.ambient.TN2=merge(subset(WY_GM_trend.ambient,(SITE2%in%wmd.sites)&variable=="TN"),
                              wmd.mon1,"SITE2")
WY_GM_trend.ambient.TN=rbind(WY_GM_trend.ambient.TN1,WY_GM_trend.ambient.TN2)

## All data trend TN Combined
WY_GM_trend.TN.all=rbind(WY_GM_trend.struct.TN,WY_GM_trend.ambient.TN)
WY_GM_trend.TN.all.shp=SpatialPointsDataFrame(WY_GM_trend.TN.all[,c("UTMX","UTMY")],
                                              data=WY_GM_trend.TN.all,proj4string=utm17)




# png(filename=paste0(plot.path,"BBSEER_WQ_AvgGM_TP.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(1,1,1,1));
layout(matrix(1:2,1,2),widths=c(1,0.5))

bks=c(0,5,10,25,50,100,150)
cols.vals=findInterval(mean.WY_FWM_GM.TP.all.shp$mean.GM,bks)
cols=viridis::inferno(length(bks)-1,alpha=0.75,direction=-1)
# cols=wesanderson::wes_palette("Zissou1",length(bks)-1,"continuous")

bbox.lims=bbox(mean.WY_FWM_GM.TP.all.shp)
plot(shore,col="cornsilk",border="grey",bg="lightblue",
     ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(wcas,add=T,col="darkseagreen2",border=NA)
plot(enp.shore,add=T,col="darkseagreen2",border="grey",lwd=0.1)
plot(ENP,add=T,col=NA,border="white",lty=2,lwd=2)
plot(subset(FDEP_AP,SHORT_NAME=="Biscayne Bay"),add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5)
plot(BNP,add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5,lty=2)
plot(canals,add=T,col="lightblue")
plot(mean.WY_FWM_GM.TP.all.shp,pch=21,bg=cols[cols.vals],cex=0.8,add=T,lwd=0.01)
plot(bbseer,add=T,border="red",lty=2)
box(lwd=1)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=0.75,seg.len=4,outer=F);
plot(0:1,0:1,ann=F,axes=F,type="n")
leg.fun(bks,cols,"Avg TP GM\n(\u03BCg L\u207B\u00B9)",
        leg.type="categorical")
mtext(side=1,line=-1.25,adj=1,
      paste0("WY",WY(dates[1])," - WY",WY(dates[2])))

dev.off()

range(mean.WY_FWM_GM.TP.all.shp$UCI)

# png(filename=paste0(plot.path,"BBSEER_WQ_GM_UCI_TP.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(1,1,1,1));
layout(matrix(1:2,1,2),widths=c(1,0.5))

bks=c(0,2,10,20,50,100,150)
cols.vals=findInterval(mean.WY_FWM_GM.TP.all.shp$UCI,bks)
cols=viridis::inferno(length(bks)-1,alpha=0.75,direction=-1)
# cols=wesanderson::wes_palette("Zissou1",length(bks)-1,"continuous")

bbox.lims=bbox(mean.WY_FWM_GM.TP.all.shp)
plot(shore,col="cornsilk",border="grey",bg="lightblue",
     ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(wcas,add=T,col="darkseagreen2",border=NA)
plot(enp.shore,add=T,col="darkseagreen2",border="grey",lwd=0.1)
plot(ENP,add=T,col=NA,border="white",lty=2,lwd=2)
plot(subset(FDEP_AP,SHORT_NAME=="Biscayne Bay"),add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5)
plot(BNP,add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5,lty=2)
plot(canals,add=T,col="lightblue")
plot(mean.WY_FWM_GM.TP.all.shp,pch=21,bg=cols[cols.vals],cex=0.8,add=T,lwd=0.01)
plot(bbseer,add=T,border="red",lty=2)
box(lwd=1)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=0.75,seg.len=4,outer=F);
plot(0:1,0:1,ann=F,axes=F,type="n")
leg.fun(bks,cols,"Upper 95% CI\nTP GM\n(\u03BCg L\u207B\u00B9)",
        leg.type="categorical")
mtext(side=1,line=-1.25,adj=1,
      paste0("WY",WY(dates[1])," - WY",WY(dates[2])))

dev.off()

# png(filename=paste0(plot.path,"BBSEER_WQ_AvgGM_TN.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(1,1,1,1));
layout(matrix(1:2,1,2),widths=c(1,0.5))

bks=c(0,0.25,0.5,1,2,5)
cols.vals=findInterval(mean.WY_FWM_GM.TN.all.shp$mean.GM,bks)
cols=viridis::inferno(length(bks)-1,alpha=0.75,direction=-1)
# cols=wesanderson::wes_palette("Zissou1",length(bks)-1,"continuous")

bbox.lims=bbox(mean.WY_FWM_GM.TN.all.shp)
plot(shore,col="cornsilk",border="grey",bg="lightblue",
     ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(wcas,add=T,col="darkseagreen2",border=NA)
plot(enp.shore,add=T,col="darkseagreen2",border="grey",lwd=0.1)
plot(ENP,add=T,col=NA,border="white",lty=2,lwd=2)
plot(subset(FDEP_AP,SHORT_NAME=="Biscayne Bay"),add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5)
plot(BNP,add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5,lty=2)
plot(canals,add=T,col="lightblue")
plot(mean.WY_FWM_GM.TN.all.shp,pch=21,bg=cols[cols.vals],cex=0.8,add=T,lwd=0.01)
plot(bbseer,add=T,border="red",lty=2)
box(lwd=1)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=0.75,seg.len=4,outer=F);
plot(0:1,0:1,ann=F,axes=F,type="n")
leg.fun(format(bks),cols,"Avg TN GM\n(mg L\u207B\u00B9)",
        leg.type="categorical")
mtext(side=1,line=-1.25,adj=1,
      paste0("WY",WY(dates[1])," - WY",WY(dates[2])))
dev.off()

# png(filename=paste0(plot.path,"BBSEER_WQ_GM_UCI_TN.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(1,1,1,1));
layout(matrix(1:2,1,2),widths=c(1,0.5))

bks=c(0,0.25,0.5,1,2,5)
cols.vals=findInterval(mean.WY_FWM_GM.TN.all.shp$UCI,bks)
cols=viridis::inferno(length(bks)-1,alpha=0.75,direction=-1)
# cols=wesanderson::wes_palette("Zissou1",length(bks)-1,"continuous")

bbox.lims=bbox(mean.WY_FWM_GM.TN.all.shp)
plot(shore,col="cornsilk",border="grey",bg="lightblue",
     ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(wcas,add=T,col="darkseagreen2",border=NA)
plot(enp.shore,add=T,col="darkseagreen2",border="grey",lwd=0.1)
plot(ENP,add=T,col=NA,border="white",lty=2,lwd=2)
plot(subset(FDEP_AP,SHORT_NAME=="Biscayne Bay"),add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5)
plot(BNP,add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5,lty=2)
plot(canals,add=T,col="lightblue")
plot(mean.WY_FWM_GM.TN.all.shp,pch=21,bg=cols[cols.vals],cex=0.8,add=T,lwd=0.01)
plot(bbseer,add=T,border="red",lty=2)
box(lwd=1)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=0.75,seg.len=4,outer=F);
plot(0:1,0:1,ann=F,axes=F,type="n")
leg.fun(format(bks),cols,"Upper 95% CI\nTN GM\n(mg L\u207B\u00B9)",
        leg.type="categorical")
mtext(side=1,line=-1.25,adj=1,
      paste0("WY",WY(dates[1])," - WY",WY(dates[2])))
dev.off()

range(WY_GM_trend.TP.all.shp$sen.slope)
subset(WY_GM_trend.TP.all.shp@data,pval<0.05&tau<0)

nrow(subset(WY_GM_trend.TP.all.shp@data,pval<0.05&tau>0))

# png(filename=paste0(plot.path,"BBSEER_WQ_AvgGM_TPTrend.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(1,1,1,1));
layout(matrix(1:2,1,2),widths=c(1,0.5))

bks=seq(-0.04,0.04,0.01)
pal2=colorRampPalette(c("blue","grey90","red"))(length(bks)-1)
cols.vals=findInterval(WY_GM_trend.TP.all.shp$sen.slope,bks)
pch.vals=with(WY_GM_trend.TP.all.shp@data,ifelse(pval>0.05,21,ifelse(pval<0.05&tau<0,25,24)))

bbox.lims=bbox(WY_GM_trend.TP.all.shp)
plot(shore,col="cornsilk",border="grey",bg="lightblue",
     ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(wcas,add=T,col="darkseagreen2",border=NA)
plot(enp.shore,add=T,col="darkseagreen2",border="grey",lwd=0.1)
plot(ENP,add=T,col=NA,border="white",lty=2,lwd=2)
plot(subset(FDEP_AP,SHORT_NAME=="Biscayne Bay"),add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5)
plot(BNP,add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5,lty=2)
plot(canals,add=T,col="lightblue")
plot(WY_GM_trend.TP.all.shp,pch=pch.vals,bg=pal2[cols.vals],cex=0.8,add=T,lwd=0.01)
plot(bbseer,add=T,border="red",lty=2)
box(lwd=1)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=0.75,seg.len=4,outer=F);

plot(0:1,0:1,ann=F,axes=F,type="n")
leg.fun(format(bks),pal2,"TP GM\nThiel-Sen Slope\n(\u03BCg L\u207B\u00B9 Y\u207B\u00B9)",
        leg.type="categorical",
        top.val=0.9,bot.val=0.4,
        x.max=0.55,x.min=0.40)
mtext(side=1,line=-1.25,adj=1,
      paste0("WY",WY(dates[1])," - WY",WY(dates[2])))
legend(0.5,0.25,legend=c("Sig. Increasing","No Trend","Sig. Decreasing"),
       pch=c(24,21,25),lty=0,lwd=c(0.1),
       pt.bg="grey",col="black",
       pt.cex=1.25,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

range(WY_GM_trend.TN.all.shp$sen.slope)

# png(filename=paste0(plot.path,"BBSEER_WQ_AvgGM_TNTrend.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(1,1,1,1));
layout(matrix(1:2,1,2),widths=c(1,0.5))

bks=seq(-0.3,0.3,0.1)
bks[4]=0;# interpolation error
pal2=colorRampPalette(c("blue","grey90","red"))(length(bks)-1)
cols.vals=findInterval(WY_GM_trend.TN.all.shp$sen.slope,bks)
pch.vals=with(WY_GM_trend.TN.all.shp@data,ifelse(pval>0.05,21,ifelse(pval<0.05&tau<0,25,24)))

bbox.lims=bbox(WY_GM_trend.TN.all.shp)
plot(shore,col="cornsilk",border="grey",bg="lightblue",
     ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(wcas,add=T,col="darkseagreen2",border=NA)
plot(enp.shore,add=T,col="darkseagreen2",border="grey",lwd=0.1)
plot(ENP,add=T,col=NA,border="white",lty=2,lwd=2)
plot(subset(FDEP_AP,SHORT_NAME=="Biscayne Bay"),add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5)
plot(BNP,add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5,lty=2)
plot(canals,add=T,col="lightblue")
plot(WY_GM_trend.TN.all.shp,pch=pch.vals,bg=pal2[cols.vals],cex=0.8,add=T,lwd=0.01)
plot(bbseer,add=T,border="red",lty=2)
box(lwd=1)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=0.75,seg.len=4,outer=F);

plot(0:1,0:1,ann=F,axes=F,type="n")
leg.fun(format(bks),pal2,"TN GM\nThiel-Sen Slope\n(mg L\u207B\u00B9 Y\u207B\u00B9)",
        leg.type="categorical",
        top.val=0.9,bot.val=0.4,
        x.max=0.55,x.min=0.40)
mtext(side=1,line=-1.25,adj=1,
      paste0("WY",WY(dates[1])," - WY",WY(dates[2])))
legend(0.5,0.25,legend=c("Sig. Increasing","No Trend","Sig. Decreasing"),
       pch=c(24,21,25),lty=0,lwd=c(0.1),
       pt.bg="grey",col="black",
       pt.cex=1.25,ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

tm_shape(WY_GM_trend.TN.all.shp)+tm_dots(col="sen.slope")
tm_shape(WY_GM_trend.TP.all.shp)+tm_dots(col="sen.slope")

subset(WY_GM_trend.TP.all.shp,SITE2=="AR03")

##
bbseer.r=raster(bbseer)
res(bbseer.r)=1000

# Spatial Trend -----------------------------------------------------------
library(fields)
library(spdep)
# TN
m=Tps(coordinates(WY_GM_trend.TN.all.shp),WY_GM_trend.TN.all.shp$sen.slope)
tps=interpolate(bbseer.r,m)
tps.TN.trend=raster::mask(tps,bbseer)
plot(tps.TN.trend)
plot(WY_GM_trend.TN.all.shp,add=T)
plot(rasterToContour(tps.TN.trend,levels=c(0.01),nlevels=2),col="yellow",lwd=2,add=T)

# TP
m=Tps(coordinates(WY_GM_trend.TP.all.shp),WY_GM_trend.TP.all.shp$sen.slope)
tps=interpolate(bbseer.r,m)
tps.TP.trend=raster::mask(tps,bbseer)
plot(tps.TP.trend)
plot(WY_GM_trend.TP.all.shp,add=T)
plot(rasterToContour(tps.TP.trend,levels=c(0.005),nlevels=2),col="yellow",lwd=2,add=T)

plot(GM~WY,subset(WY_FWM_GM.struct,variable=="TP"&SITE2=="S178"))



# NNC eval ----------------------------------------------------------------

STORET.WIN.locs.shp
locs.NNC=data.frame(sf::st_intersection(sf::st_as_sf(STORET.WIN.locs.shp),sf::st_as_sf(nnc)))

unique(locs.NNC$ESTUARY)
subset(nnc,ESTUARY=="Biscayne Bay")@data
subset(nnc,ESTUARY=="Biscayne Bay")

struct.wq.nnc=merge(struct.wq,locs.NNC[,c("Station.ID","ESTUARY_SE","ESTUARY")],by.x="Station.ID",by.y="Station.ID")
      

struct.wq.nnc.GM=ddply(subset(struct.wq.nnc,variable%in%c("TN","TP")&screen==1),c("Station.ID","WY","variable","ESTUARY_SE","ESTUARY"),summarise,GM=exp(mean(log(value),na.rm=T)))
NNC.meanGM=dcast(struct.wq.nnc.GM,ESTUARY_SE+ESTUARY+WY~variable,value.var="GM",mean,na.rm=T)

fill=data.frame(expand.grid(ESTUARY_SE=paste0("ENRH",1:9),
                 ESTUARY="Biscayne Bay",
                 WY=seq(2007,2022,1)))
NNC.meanGM=merge(NNC.meanGM,fill,c("ESTUARY_SE","ESTUARY","WY"),all.y=T)

NNC.meanGM=merge(NNC.meanGM,subset(nnc,ESTUARY=="Biscayne Bay")@data[,c("ESTUARY_SE","TN_CRITERI","TP_CRITERI")],"ESTUARY_SE")
NNC.meanGM$TP.ann.exceed=with(NNC.meanGM,
                          ifelse(TP>TP_CRITERI,1,0))
NNC.meanGM$TP.exceed.3yr=with(NNC.meanGM,ave(TP.ann.exceed,ESTUARY_SE,FUN=function(x) zoo::rollapply(x,3,sum,na.rm=T,align="right",fill=NA)))
NNC.meanGM$TP.LT.exceed=with(NNC.meanGM,ifelse(TP.exceed.3yr==3,1,0))

NNC.meanGM$TN.ann.exceed=with(NNC.meanGM,
                          ifelse(TN>TN_CRITERI,1,0))
NNC.meanGM$TN.exceed.3yr=with(NNC.meanGM,ave(TN.ann.exceed,ESTUARY_SE,FUN=function(x) zoo::rollapply(x,3,sum,na.rm=T,align="right",fill=NA)))
NNC.meanGM$TN.LT.exceed=with(NNC.meanGM,ifelse(TN.exceed.3yr==3,1,0))

NNC.seg=paste0("ENRH",1:9)
# png(filename=paste0(plot.path,"BBSEER_WQ_WQSEval.png"),width=8,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,1,0.5,0.5),oma=c(3,2.5,0.75,0.25),lwd=0.5);
layout(matrix(1:18,2,9,byrow = T))

xlim.val=c(2007,2022);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0,0.015);by.y=0.005;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
for(i in 1:length(NNC.seg)){
  tmp=subset(NNC.meanGM,ESTUARY_SE==NNC.seg[i])
  nnc.tmp=subset(nnc@data,ESTUARY_SE==NNC.seg[i])
  plot(TP~WY,tmp,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
  abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
  with(tmp,pt_line(WY,TP,2,"dodgerblue1",1,21,
                   ifelse(TP.exceed==1,"red","dodgerblue1"),pt.lwd = 0.01,cex=1.25))
  abline(h=nnc.tmp$TP_CRITERI,lty=2,col="red",lwd=1.5)
  # text(xlim.val[1]+0.5,nnc.tmp$TP_CRITERI,paste0(nnc.tmp$TP_CRITERI*1000," \u03BCg L\u207B\u00B9"),font=3,pos=3,offset=0.1,col="red")
  axis_fun(1,xmaj,xmin,NA,line=-0.5)
  if(i==1){axis_fun(2,ymaj,ymin,format(ymaj*1000));}else{axis_fun(2,ymaj,ymin,NA);}
  box(lwd=1)
  if(i==1){mtext(side=2,outer=F,line=2,"Avg. TP GM (\u03BCg L\u207B\u00B9)",cex=0.8)}
  # mtext(side=1,line=2.6,"WY")
  # mtext(side=3,adj=0,line=-0.5,paste0(" ",nnc.tmp$SEGMENT_NA,"\n(",nnc.tmp$ESTUARY_SE,")"),cex=0.6,padj=1)
  mtext(side=3,adj=0,paste0( nnc.tmp$ESTUARY_SE),cex=0.6)
}

ylim.val=c(0,0.9);by.y=0.3;ymaj=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(max(c(0,ylim.val[1])),ylim.val[2],by.y/2)
for(i in 1:length(NNC.seg)){
  tmp=subset(NNC.meanGM,ESTUARY_SE==NNC.seg[i])
  nnc.tmp=subset(nnc@data,ESTUARY_SE==NNC.seg[i])
  plot(TN~WY,tmp,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
  abline(h=ymaj,v=xmaj,lty=3,col="grey",lwd=0.5)
  with(tmp,pt_line(WY,TN,2,"darkolivegreen3",1,21,
                   ifelse(TN.exceed==1,"red","darkolivegreen3"),pt.lwd = 0.01,cex=1.25))
  abline(h=nnc.tmp$TN_CRITERI,lty=2,col="red",lwd=1.5)
  # text(xlim.val[1]+0.75,nnc.tmp$TN_CRITERI,paste0(nnc.tmp$TN_CRITERI," mg L\u207B\u00B9"),font=3,pos=3,offset=0.2,col="red")
  axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
  if(i==1){axis_fun(2,ymaj,ymin,format(ymaj));}else{axis_fun(2,ymaj,ymin,NA);}
  box(lwd=1)
  if(i==1){mtext(side=2,outer=F,line=2,"Avg. TN GM (mg L\u207B\u00B9)",cex=0.8)}
  if(i==9){mtext(side=1,line=2.6,"WY")}
  # mtext(side=3,adj=0,line=-1.25,paste0(" ",nnc.tmp$SEGMENT_NA," (",nnc.tmp$ESTUARY_SE,")"),cex=0.75)
}
mtext(side=1,line=1,outer=T,"Water Year")
dev.off()


exceed.prop=ddply(NNC.meanGM,"ESTUARY_SE",summarise,
                  TP.exceed=(sum(TP.LT.exceed,na.rm=T)/N.obs(TP.LT.exceed))*100,
                  TN.exceed=(sum(TN.LT.exceed,na.rm=T)/N.obs(TN.LT.exceed))*100)

nnc.exceed=merge(subset(nnc,ESTUARY=="Biscayne Bay"),exceed.prop,"ESTUARY_SE")


plot(nnc.exceed)

bks=seq(0,100,20)
cols.vals=findInterval(nnc.exceed$TN.exceed,bks)
cols=colorRampPalette(c("lightgreen","yellow","indianred"))(length(bks)-1)

plot(nnc.exceed,col=cols[cols.vals])


cols.vals=findInterval(nnc.exceed$TP.exceed,bks)
cols=colorRampPalette(c("lightgreen","yellow","indianred"))(length(bks)-1)

plot(nnc.exceed,col=cols[cols.vals])


# png(filename=paste0(plot.path,"BBSEER_WQ_NNCExceed.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(1,1,2,1));
layout(matrix(1:3,1,3),widths=c(1,1,0.5))

bks=seq(0,100,20)
cols.vals=findInterval(nnc.exceed$TP.exceed,bks)
cols=colorRampPalette(c("lightgreen","yellow","indianred"))(length(bks)-1)

bbox.lims=bbox(nnc.exceed)
plot(shore,col="cornsilk",border="grey",bg="lightblue",
     ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(wcas,add=T,col="grey",border=NA)
plot(enp.shore,add=T,col="grey",border="grey",lwd=0.1)
plot(ENP,add=T,col=NA,border="white",lty=2,lwd=2)
plot(canals,add=T,col="lightblue")
plot(nnc.exceed,col=cols[cols.vals],add=T,lwd=0.1)
raster::text(nnc.exceed,nnc.exceed$ESTUARY_SE,halo=T,col="red",font=2,cex=0.75)
plot(bbseer,add=T,border="red",lty=2)
box(lwd=1)
mtext(side=3,"Total Phosphorous")
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=0.75,seg.len=4,outer=F);

cols.vals=findInterval(nnc.exceed$TN.exceed,bks)
cols=colorRampPalette(c("lightgreen","yellow","indianred"))(length(bks)-1)

plot(shore,col="cornsilk",border="grey",bg="lightblue",
     ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(wcas,add=T,col="grey",border=NA)
plot(enp.shore,add=T,col="grey",border="grey",lwd=0.1)
plot(ENP,add=T,col=NA,border="white",lty=2,lwd=2)
plot(canals,add=T,col="lightblue")
plot(nnc.exceed,col=cols[cols.vals],add=T,lwd=0.1)
raster::text(nnc.exceed,nnc.exceed$ESTUARY_SE,halo=T,col="red",font=2,cex=0.75)
plot(bbseer,add=T,border="red",lty=2)
box(lwd=1)
mtext(side=3,"Total Nitrogen")

plot(0:1,0:1,ann=F,axes=F,type="n")
leg.fun(format(bks),cols,paste0("Percent NNC Exceedance\n WY",WY(dates[1])," - WY",WY(dates[2])),
        leg.type="categorical",
        top.val=0.7,bot.val=0.3,
        x.max=0.55,x.min=0.40)

dev.off()



# Cluster Analysis --------------------------------------------------------
## See if we can idetifiy WQ regions. 
# https://www.r-bloggers.com/2021/04/cluster-analysis-in-r/

struct.wq2=ddply(subset(struct.wq,variable%in%c("TN","TP")&screen==1),c("Station.ID","WY","variable"),summarise,GM=exp(mean(log(value),na.rm=T)))
struct.wq2=dcast(struct.wq2,Station.ID~variable,value.var="GM",mean,na.rm=T)

z <- struct.wq2[,c("TP","TN")]
means <- apply(z,2,mean,na.rm=T)
sds <- apply(z,2,sd,na.rm=T)
nor <- scale(z,center=means,scale=sds)

distance = dist(nor)

mydata.hclust = hclust(distance)
plot(mydata.hclust)
plot(mydata.hclust,labels=struct.wq2$Station.ID,main='Default from hclust')
plot(mydata.hclust,hang=-1, labels=struct.wq2$Station.ID,main='Default from hclust')

mydata.hclust<-hclust(distance,method="average") 
plot(mydata.hclust,hang=-1) 

member = cutree(mydata.hclust,6)
table(member)
member

struct.wq2$clust1=member

boxplot(TP~clust1,struct.wq2)
boxplot(TN~clust1,struct.wq2)

STORET.WIN.sites2=cbind(data.frame(Station.ID=STORET.WIN.sites@data$Station.ID),coordinates(STORET.WIN.sites))
colnames(STORET.WIN.sites2)=c("Station.ID","UTMX","UTMY")

colnames(wmd.struct2)=c("Station.ID","UTMX","UTMY")

sites.all=rbind(STORET.WIN.sites2,wmd.struct2)
struct.wq2=merge(struct.wq2,sites.all,"Station.ID")
struct.wq2.shp=SpatialPointsDataFrame(struct.wq2[,c("UTMX","UTMY")],
                                      data=struct.wq2,
                                      proj4string = utm17)

# png(filename=paste0(plot.path,"BBSEER_WQ_AvgGM_TNTrend.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(1,1,1,1));
layout(matrix(1:2,1,2),widths=c(1,0.5))

pal2=viridis::cividis(max(member))

bbox.lims=bbox(struct.wq2.shp)
plot(shore,col="cornsilk",border="grey",bg="lightblue",
     ylim=bbox.lims[c(2,4)],xlim=bbox.lims[c(1,3)],lwd=0.1)
plot(wcas,add=T,col="darkseagreen2",border=NA)
plot(enp.shore,add=T,col="darkseagreen2",border="grey",lwd=0.1)
plot(ENP,add=T,col=NA,border="white",lty=2,lwd=2)
plot(subset(FDEP_AP,SHORT_NAME=="Biscayne Bay"),add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5)
plot(BNP,add=T,col=adjustcolor("darkseagreen2",0.5),border="white",lwd=0.5,lty=2)
plot(canals,add=T,col="lightblue")
plot(struct.wq2.shp,pch=21,bg=pal2[as.numeric(struct.wq2.shp$clust1)],cex=0.8,add=T,lwd=0.01)
plot(bbseer,add=T,border="red",lty=2)
box(lwd=1)
mapmisc::scaleBar(utm17,"bottomleft",bty="n",cex=0.75,seg.len=4,outer=F);
