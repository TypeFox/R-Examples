source("Methods.R")
library(flare)
library(lasso2)

## Set seed
Seed  <- 1234
nfold <-    5
set.seed(Seed)

## Prostate data
data(Prostate)
ProsDataY    <- Prostate[,9]
ProsDataX    <- as.matrix(Prostate[,-9])
nP           <- nrow(Prostate)
dp           <- round(nP/nfold)
CVpr         <- NULL

for(ifold in 1:nfold){
  lfold  <- (1+(ifold-1)*dp):(ifold*dp)
  xt     <- ProsDataX[-lfold,]; yt <- ProsDataY[-lfold]
  xv     <- ProsDataX[+lfold,]; yv <- ProsDataY[+lfold]
  cv     <- compare(xt,yt,xv,yv,Seed)
  CVpr   <- rbind(CVpr,cv)
}

## Eye Data
data(eyedata)
EyeDataX <- x
EyeDataY <- y
nE       <- nrow(EyeDataX)
dp       <- round(nE/nfold)
CVed     <- NULL

for(ifold in 1:nfold){
  lfold  <- (1+(ifold-1)*dp):(ifold*dp)
  xt     <- EyeDataX[-lfold,]; yt <- EyeDataY[-lfold]
  xv     <- EyeDataX[+lfold,]; yv <- EyeDataY[+lfold]
  cv     <- compare(xt,yt,xv,yv,Seed)
  CVed   <- rbind(CVed,cv)
}


numExpRealData <- list(ProstateData=CVpr,EyeData=CVed)
save(list=c("numExpRealData"),file="numExpRealData.RData")
save.image("RealDataExample.RData")


## Create table
createTable2 <- function(numExpRealData){
  pData <- numExpRealData[["ProstateData"]][,1:27]
  eData <- numExpRealData[["EyeData"]][,1:27]
  basenames <- colnames(pData)[1:9]
  colnames(pData) <- paste(basenames,rep(c("cv","df","tt"),each=9),sep="_")
  colnames(eData) <- paste(basenames,rep(c("cv","df","tt"),each=9),sep="_")
  methname <- basenames; names(methname) <- methname
  methname["ELNET"]   <- "Elastic net"
  methname["CLERE_s"] <- "CLERE$_0$"
  methname["SS"]      <- "Spike and Slab"
  
  
  genLine <- function(refData,meth){
    cv <- 100*refData[,paste(meth,"_cv",sep="")]
    df <- refData[,paste(meth,"_df",sep="")]
    tt <- refData[,paste(meth,"_tt",sep="")]
    
    line <- paste(
                  methname[meth],"      &  ",format(mean(cv),digit=2,nsmall=1)," (",format(sd(cv)/sqrt(5),digit=1)," )",
                  "                     &  ",format(mean(df),digit=2,nsmall=1)," (",format(sd(df)/sqrt(5),digit=1)," )",
                  "                     &  ",format(mean(tt),digit=2,nsmall=1)," (",format(sd(tt)/sqrt(5),digit=1)," ) \\\\\n")
    return(line)
  }
  tab <- paste("\\tiny\n",
               "\\begin{center}\n",
               "\\begin{table}[h!]\n",
               "\\begin{tabular}{lccc}\n",
               "\\toprule\n",
               "                    & 100$\\times$Averaged CV-statistic & Number of parameters & CT (seconds)\\\\\n",
               "                    & (Std. Error)                     & (Std. Error)         & (Std. Error) \\\\\n",
               "\\midrule\n",
               "\\code{Prostate dataset}                   &                       &                               &  \\\\\n",
               "\\midrule\n",
               paste(sapply(basenames,function(meth) genLine(pData,meth)),collapse=""),
               "\\\\\n",
               "\\midrule\n",
               "\\code{eyedata}                    &                       &                               &  \\\\\n",
               "\\midrule\n",
               paste(sapply(basenames,function(meth) genLine(eData,meth)),collapse=""),
               "\\bottomrule\n",
               "\\end{tabular}\n",
               "\\caption{\\label{tab:realdata} Real data analysis. Out-of-sample prediction error (averaged CV-statistic) was estimated using cross-validation in 100 splitted datasets. The number of parameters reported for CLERE/CLERE$_0$ was selected using AIC. CT stands for the average Computational Time.}\n",
               "\\end{table}\n",
               "\\end{center}\n",
               "\\normalsize\n")
  
  return(tab)  
}
cat(createTable2(numExpRealData))
