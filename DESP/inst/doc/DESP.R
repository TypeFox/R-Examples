## ----VersionDate,echo=FALSE,message=FALSE,results='hide'--------------
options(width=72)
knitr::opts_knit$set(width=72)
knitr::opts_chunk$set(fig.pos = "", fig.align = "center")
set.seed(0)
library(DESP, quietly=TRUE)
DESPversion <- packageDescription("DESP")$Version
DESPdateRaw <- packageDescription("DESP")$Date
DESPdateYear <- as.numeric(substr(DESPdateRaw, 1, 4))
DESPdateMonth <- as.numeric(substr(DESPdateRaw, 6, 7))
DESPdateDay <- as.numeric(substr(DESPdateRaw, 9, 10))
DESPdate <- paste0(month.name[DESPdateMonth], " ",
                     DESPdateDay, ", ",
                     DESPdateYear)

## ----Install,eval=FALSE-----------------------------------------------
#  install.packages("DESP")

## ----Load,eval=FALSE--------------------------------------------------
#  library(DESP)
#  data(iris3)

## ----ExFun------------------------------------------------------------
settings <- list(diagElem='AD')
v <- 5

## ----ExRun,message=FALSE,warning=FALSE--------------------------------
set.seed(1)
categories <- colnames(iris3[1,,])
params <- vector(mode="list", length=length(categories))
for(c in 1:length(categories)){
  obs <- 1:25
  lr <- (9/10)^(0:9)
  gr <- (1/sqrt(2))^(0:4) * sqrt(2*log(ncol(iris3[,,c])))
  params[[c]] <- desp.cv(iris3[obs,,c], v=v, lambda.range=lr, 
    gamma.range=gr, settings=settings)
}

## ----ExRes------------------------------------------------------------
params[[1]]

