

## Download BBSEER Round 2 Modeling
##
## Code was compiled by Paul Julian
## contact info: pjulian@sccf.org

## BAD 
## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
## Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

## Libraries
#devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);
library(RCurl)


# -------------------------------------------------------------------------
## ECB and FWO
link="ftp://ftppub.sfwmd.gov/outgoing/BBSEER/20230201/PM_ECB22_FWO_FWOI/Model_Output/"

result=getURL(link)
result2=strsplit(result, "\r*\n")[[1]]
result2=strsplit(result2,"\\s+")

alt=sapply(result2,"[",9)
alt=alt[3:length(alt)]
alt=alt[!(alt%in%c("FWOI"))]

data.path1="C:/Julian_LaCie/_GitHub/BBSEER_WQ/Data/Round2_RSMGL"
Folder.Maker(paste(data.path1,alt,sep="/"))

for(i in 1:length(alt)){
  url=paste0(link,alt[i],"/RSMGL/")
  dest=paste0(data.path1,"/",alt[i],"/RSMGL_output.dss")
  download.file(paste0(url,"RSMGL_output.dss"),dest,mode="wb",cacheOK = F)
  print(i)
}


####
link="ftp://ftppub.sfwmd.gov/outgoing/BBSEER/20230314/PM_FWOi_ALT21_ALT22_ALT23_ALT24/Model_Output/"

result=getURL(link)
result2=strsplit(result, "\r*\n")[[1]]
result2=strsplit(result2,"\\s+")

alt=sapply(result2,"[",9)
alt=alt[3:length(alt)]

Folder.Maker(paste(data.path1,alt,sep="/"))
for(i in 1:length(alt)){
  url=paste0(link,alt[i],"/RSMGL/")
  dest=paste0(data.path1,"/",alt[i],"/RSMGL_output.dss")
  download.file(paste0(url,"RSMGL_output.dss"),dest,mode="wb",cacheOK = F)
  print(i)
}
i=1

url=paste0(link,alt[i],"/RSMGL/")
dest=paste0(data.path1,"/",alt[i],"/RSMGL_output.dss")
download.file(paste0(url,"RSMGL_output.dss"),dest,mode="wb",cacheOK = F)
print(i)



## -------------------------------------------------------------------------
## ECB and FWO
link="ftp://ftppub.sfwmd.gov/outgoing/BBSEER/20230201/PM_ECB22_FWO_FWOI/Model_Output/"

result=getURL(link)
result2=strsplit(result, "\r*\n")[[1]]
result2=strsplit(result2,"\\s+")

alt=sapply(result2,"[",9)
alt=alt[3:length(alt)]
alt=alt[!(alt%in%c("FWOI"))]

data.path1="C:/Julian_LaCie/_GitHub/BBSEER_WQ/Data/Round2_RSMGL"
Folder.Maker(paste(data.path1,alt,sep="/"))

for(i in 1:length(alt)){
  url=paste0(link,alt[i],"/RSMGL/")
  dest=paste0(data.path1,"/",alt[i],"/transect_flows.dss")
  download.file(paste0(url,"transect_flows.dss"),dest,mode="wb",cacheOK = F)
  print(i)
}


####
link="ftp://ftppub.sfwmd.gov/outgoing/BBSEER/20230314/PM_FWOi_ALT21_ALT22_ALT23_ALT24/Model_Output/"

result=getURL(link)
result2=strsplit(result, "\r*\n")[[1]]
result2=strsplit(result2,"\\s+")

alt=sapply(result2,"[",9)
alt=alt[3:length(alt)]
alt=alt[!(alt%in%c("WaterBudgets"))]

Folder.Maker(paste(data.path1,alt,sep="/"))
for(i in 1:length(alt)){
  url=paste0(link,alt[i],"/RSMGL/")
  dest=paste0(data.path1,"/",alt[i],"/transect_flows.dss")
  download.file(paste0(url,"transect_flows.dss"),dest,mode="wb",cacheOK = F)
  print(i)
}


# i=1
# url=paste0(link,alt[i],"/RSMGL/")
# dest=paste0(data.path1,"/",alt[i],"/transect_flows.dss")
# download.file(paste0(url,"transect_flows.dss"),dest,mode="wb",cacheOK = F)
# print(i)

