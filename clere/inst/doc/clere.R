### R code from vignette source 'clere.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: clere.Rnw:591-655
###################################################
library(clere)
data(algoComp)
  nsim   <- 200
  lnames <- names(algoComp)
  for(i in 1:length(lnames)){
    Tmp <- algoComp[[i]][,1:4]
    colnames(Tmp) <- c("MCEM5","MCEM25","MCEM125","SEM")
    if(lnames[i]=="Liks"){
      tmp <- format(100*table(factor(apply(Tmp,1,which.max),levels=1:4))/nrow(Tmp),digit=2,nsmall=1)
      names(tmp) <- colnames(Tmp)
    }else{
      tmp <- 100*table(factor(apply(Tmp,1,which.min),levels=1:4))/nrow(Tmp)
      names(tmp) <- colnames(Tmp)
    }
    assign(paste("isBest",lnames[i],sep=""),tmp)
    ## Mean
    avTmp <- format(apply(Tmp,2,median),digit=2,nsmall=1)
    sdTmp <- format(1.253*apply(Tmp,2,sd)/sqrt(nsim),digit=2,nsmall=1)
    assign(paste("av",lnames[i],sep=""),avTmp)
    assign(paste("sd",lnames[i],sep=""),sdTmp)
  }
  avTPpred <- format( median(algoComp$Pred[,"ORACLE"]),digit=2,nsmall=1)
  sdTPpred <- format( 1.253*sd(algoComp$Pred[,"ORACLE"])/sqrt(nsim),digit=2,nsmall=1)
  avTPliks <- format( median(algoComp$Liks[,"ORACLE"]),digit=2,nsmall=1)
  sdTPliks <- format( 1.253*sd(algoComp$Liks[,"ORACLE"])/sqrt(nsim),digit=2,nsmall=1)
  
  tab <- paste(
               "\\begin{center}\n",
               "\\begin{table}[h!]\n",
               "\\begin{tabular}{llrr}",
               "\\toprule\n",
               "              &                  & \\% of times                    & Median       \\\\\n",
               "Performance indicators & Algorithms & the algorithm was best       & (Std. Err.)\\\\\n",
               "\\midrule\n",
               "CT (seconds)  &  SEM             &  ",isBestTime["SEM"],"           &  ",avTime["SEM"]," ( ",sdTime["SEM"]," ) \\\\\n",
               "              &  MCEM$_5$        &  ",isBestTime["MCEM5"],"         &  ",avTime["MCEM5"]," ( ",sdTime["MCEM5"]," ) \\\\\n",
               "              &  MCEM$_{25}$     &  ",isBestTime["MCEM25"],"        &  ",avTime["MCEM25"]," ( ",sdTime["MCEM25"]," ) \\\\\n",
               "              &  MCEM$_{125}$    &  ",isBestTime["MCEM125"],"       &  ",avTime["MCEM125"]," ( ",sdTime["MCEM125"]," ) \\\\\n",
               "\\\\\n",
               "\\midrule\n",
               "MSEE          &  SEM             &  ",isBestBias["SEM"],"           &  ",avBias["SEM"]," ( ",sdBias["SEM"]," ) \\\\\n",
               "              &  MCEM$_5$        &  ",isBestBias["MCEM5"],"         &  ",avBias["MCEM5"]," ( ",sdBias["MCEM5"]," ) \\\\\n",
               "              &  MCEM$_{25}$     &  ",isBestBias["MCEM25"],"        &  ",avBias["MCEM25"]," ( ",sdBias["MCEM25"]," ) \\\\\n",
               "              &  MCEM$_{125}$    &  ",isBestBias["MCEM125"],"       &  ",avBias["MCEM125"]," ( ",sdBias["MCEM125"]," ) \\\\\n",
               "\\\\\n",
               "\\midrule\n",
               "MSPE          &  SEM             &  ",isBestPred["SEM"],"           &  ",avPred["SEM"]," ( ",sdPred["SEM"]," ) \\\\\n",
               "              &  MCEM$_5$        &  ",isBestPred["MCEM5"],"         &  ",avPred["MCEM5"]," ( ",sdPred["MCEM5"]," ) \\\\\n",
               "              &  MCEM$_{25}$     &  ",isBestPred["MCEM25"],"        &  ",avPred["MCEM25"]," ( ",sdPred["MCEM25"]," ) \\\\\n",
               "              &  MCEM$_{125}$    &  ",isBestPred["MCEM125"],"       &  ",avPred["MCEM125"]," ( ",sdPred["MCEM125"]," ) \\\\\n",                             
               "              &  True parameter  &  ---                             &  ",avTPpred," (",sdTPpred," )\\\\\n",
               "\\\\\n",
               "\\midrule\n",
               "ML            &  SEM             &  ",isBestLiks["SEM"],"           &  ",avLiks["SEM"]," ( ",sdLiks["SEM"]," ) \\\\\n",
               "              &  MCEM$_5$        &  ",isBestLiks["MCEM5"],"         &  ",avLiks["MCEM5"]," ( ",sdLiks["MCEM5"]," ) \\\\\n",
               "              &  MCEM$_{25}$     &  ",isBestLiks["MCEM25"],"        &  ",avLiks["MCEM25"]," ( ",sdLiks["MCEM25"]," ) \\\\\n",
               "              &  MCEM$_{125}$    &  ",isBestLiks["MCEM125"],"       &  ",avLiks["MCEM125"]," ( ",sdLiks["MCEM125"]," ) \\\\\n",        
               "              &  True parameter  &  ---                             &  ",avTPliks," (",sdTPliks," )\\\\\n",
               "\\bottomrule\n",
               "\\end{tabular}\n",
               "\\caption{\\label{tab:simulations} Performance indicators used to compare SEM and MCEM algorithms. Computational Time (CT) was measured on a Intel(R) Xeon(R) CPU E7- 4870  @ 2.40GHz processor. The best algorithm is defined as the one that either reached the largest log-likelihood (ML) or the lowest CT, Mean Squared Prediction Error (MSPE) and Mean Squared Estimation Error (MSEE).}\n",
               "\\end{table}\n",
               "\\end{center}\n")
cat(tab)  


###################################################
### code chunk number 2: clere.Rnw:752-770
###################################################
library(clere)
data(numExpSimData)
meths <- c("CLERE0","CLERE","PACS","LASSO",
           "AVG","Ridge","Elastic net",
           "Spike and Slab")
o    <- order(apply(numExpSimData[,1:9],2,median))
dfs  <- round(apply(numExpSimData[,10:18][,o[1:8]],2,mean),1)
sdf  <- round(apply(numExpSimData[,10:18][,o[1:8]],2,sd),1)
tts  <- round(apply(numExpSimData[,19:27][,o[1:8]],2,mean),1)
stt  <- round(apply(numExpSimData[,19:27][,o[1:8]],2,sd),2)
cols <- rainbow(9)
par(mar=c(5, 2, 4, 7)+0.1)
boxplot(numExpSimData[,1:9][,o[1:8]],horizontal=TRUE,log="x",col=cols,axes=FALSE,pch=18,xlab="Mean Squared Prediction Error")
axis(1)
labs <- paste(meths,"\ndf: ",dfs," (",sdf,")",sep="")
axis(4,at=1:8,labels=labs,las=2)
cts <- paste(tts,"s (",stt,")",sep="")
legend("topleft",legend=cts,box.lty=0,lwd=2,lty=1,col=cols,title="Computational time")


###################################################
### code chunk number 3: clere.Rnw:814-818
###################################################
library(clere)
data(Prostate)
y <- Prostate[,"lpsa"]
x <- as.matrix(Prostate[,-which(colnames(Prostate)=="lpsa")]) 


###################################################
### code chunk number 4: clere.Rnw:821-824
###################################################
itraining <- 1:(0.8*nrow(x))
xt <- x[ itraining,] ; yt <- y[ itraining]
xv <- x[-itraining,] ; yv <- y[-itraining]


###################################################
### code chunk number 5: clere.Rnw:827-832
###################################################
Seed <- 1234
mod  <- fitClere(y=yt,x=xt,g=5,analysis="aic",parallel=FALSE,nstart=5,
                 sparse=TRUE,nItEM=2000,nBurn=1000,nItMC=10,dp=5,nsamp=1000,
                 seed=Seed)
summary(mod)


###################################################
### code chunk number 6: clere.Rnw:836-837
###################################################
plot(mod)


###################################################
### code chunk number 7: clere.Rnw:844-845
###################################################
clusters(mod,thresold=0.7)


###################################################
### code chunk number 8: clere.Rnw:848-849
###################################################
mod@P


###################################################
### code chunk number 9: clere.Rnw:852-854
###################################################
error <- mean( (yv - predict(mod,xv))^2 )
error


###################################################
### code chunk number 10: clere.Rnw:899-943
###################################################
library(clere)
data(numExpRealData)
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
cat(tab)


