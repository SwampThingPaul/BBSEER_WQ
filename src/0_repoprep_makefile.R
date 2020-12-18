## 
## BBSEER Repo set-up and presentation makefile
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

wd="C:/Julian_LaCie/_GitHub/BBSEER_WQ"

paths=paste0(wd,c("/Plots/","/Exports/","/Data/","/src/","/GIS"))
# Folder.Maker(paths);#One and done. Creates folders in working DIR


# -------------------------------------------------------------------------
# Presentation makefile

# Render slides
# rmarkdown::render("Julian_BBSEER_OFW.Rmd")

list.files(full.names=T)
files=c("./libs/","./Plots/","./resources/","./Julian_BBSEER_OFW.html","./BBSEER_OFW_map.html")

# local webpage
webpage.loc="c:/Julian_LaCie/_GitHub/SwampThingPaul.github.io/slides/BBSEER"
# Folder.Maker(webpage.loc)
file.copy(files,webpage.loc,overwrite=T,recursive=T)


## BBSEER WQ planning targets
files=c("./libs/","./Plots/","./resources/","./BBSEER_OFW_map.html","./BBSEER_WQPlanTarget.html")
# local webpage
webpage.loc="c:/Julian_LaCie/_GitHub/owper-tech.github.io/slides/BBSEER"
# Folder.Maker(webpage.loc)
file.copy(files,webpage.loc,overwrite=T,recursive=T)


