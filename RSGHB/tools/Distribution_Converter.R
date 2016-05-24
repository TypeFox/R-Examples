library(RSGHB)

#################################################################
# This script calculates the moments and median and mode for the
# random parameters in your mixed logit model produced by RSGHB
# Current only works for normals and positive/negative log-normals
# Johnson SB and truncated normals are not supported
#
# The user of the script needs to:
# 1. set the working directory
# 2. read in the A and D files from RSGHB
# 3. copy in the gDist vector from model control file.
#################################################################


setwd("C:\\Work\\Code\\HB\\Tools\\Log-Normal Converter")

Afile <- read.table("MNL_withLogNormals_A.csv",sep=",",header=T)
Dfile <- read.table("MNL_withLogNormals_D.csv",sep=",",header=T)

gDIST <- c(2,3)


#################################################################
#################################################################
########### nothing below this line should be modified ##########
#################################################################
#################################################################

numIters <- dim(Afile)[1]
numParms <- dim(Afile)[2]-1

LNdists <- (1:numParms)[(gDIST==2|gDIST==3)*(1:numParms)]
Ndists  <- (1:numParms)[(gDIST==1)*(1:numParms)]

meansMat     <- matrix(0,nrow=numIters,numParms+1)
mediansMat   <- matrix(0,nrow=numIters,numParms+1)
modesMat     <- matrix(0,nrow=numIters,numParms+1)
variancesMat <- matrix(0,nrow=numIters,numParms+1)
skewnessMat  <- matrix(0,nrow=numIters,numParms+1)
kurtosisMat  <- matrix(0,nrow=numIters,numParms+1)

colnames(meansMat)     <- c("iteration",names(Afile)[-1])
colnames(mediansMat)   <- c("iteration",names(Afile)[-1])
colnames(modesMat)     <- c("iteration",names(Afile)[-1])
colnames(variancesMat) <- c("iteration",names(Afile)[-1])
colnames(skewnessMat)  <- c("iteration",names(Afile)[-1])
colnames(kurtosisMat)  <- c("iteration",names(Afile)[-1])

LNmetrics <- function(iter,LNdists)
{

     A.i <- as.matrix(Afile[iter,c(1+LNdists)])
     D.i <- as.matrix(Dfile[iter,-1])

     meansMat[iter,c(1,1+LNdists)]     <<- c(iter,exp(A.i+diag(xpnd(D.i))/2))
     mediansMat[iter,c(1,1+LNdists)]   <<- c(iter,exp(A.i))
     modesMat[iter,c(1,1+LNdists)]     <<- c(iter,exp(A.i-diag(xpnd(D.i))))
     variancesMat[iter,c(1,1+LNdists)] <<- c(iter,(exp(diag(xpnd(D.i))[LNdists])-1)*exp(2*A.i+diag(xpnd(D.i))[LNdists]))
     skewnessMat[iter,c(1,1+LNdists)]  <<- c(iter,(exp(diag(xpnd(D.i))[LNdists])+2)*sqrt(exp(diag(xpnd(D.i))[LNdists])-1))
     kurtosisMat[iter,c(1,1+LNdists)]  <<- c(iter,exp(4*diag(xpnd(D.i))[LNdists])+2*exp(3*diag(xpnd(D.i))[LNdists])+3*exp(2*diag(xpnd(D.i))[LNdists])-6)

}

Nmetrics <- function(iter,Ndists)
{
     
     A.i <- as.matrix(Afile[iter,c(1+Ndists)])
     D.i <- as.matrix(Dfile[iter,-1])
     
     meansMat[iter,c(1,1+Ndists)]     <<- c(iter,A.i)
     mediansMat[iter,c(1,1+Ndists)]   <<- c(iter,A.i)
     modesMat[iter,c(1,1+Ndists)]     <<- c(iter,A.i)
     variancesMat[iter,c(1,1+Ndists)] <<- c(iter,diag(xpnd(D.i))[Ndists])
     skewnessMat[iter,c(1,1+Ndists)]  <<- c(iter,0)
     kurtosisMat[iter,c(1,1+Ndists)]  <<- c(iter,0)
     
}

for(iter in 1:numIters)
{
     if(!is.na(LNdists[1]))
     {
          LNmetrics(iter,LNdists)
     }
     if(!is.na(Ndists[1]))
     {
          Nmetrics(iter,Ndists)
     }
}

calcStats <- function(mat)
{
     c(colMeans(mat[,-1]),apply(mat[,-1],2,sd))
}

resultsMat <- matrix(0,nrow=6,ncol=2*numParms)

resultsMat[1,] <- calcStats(meansMat)*c(((gDIST!=3)*1 + (gDIST==3)*-1),rep(1,length(gDIST)))
resultsMat[2,] <- calcStats(mediansMat)*c(((gDIST!=3)*1 + (gDIST==3)*-1),rep(1,length(gDIST)))
resultsMat[3,] <- calcStats(modesMat)*c(((gDIST!=3)*1 + (gDIST==3)*-1),rep(1,length(gDIST)))
resultsMat[4,] <- calcStats(variancesMat)
resultsMat[5,] <- calcStats(skewnessMat)*c(((gDIST!=3)*1 + (gDIST==3)*-1),rep(1,length(gDIST)))
resultsMat[6,] <- calcStats(kurtosisMat)*c(((gDIST!=3)*1 + (gDIST==3)*-1),rep(1,length(gDIST)))

rownames(resultsMat) <- c("mean","median","mode","variance","skewness","kurtosis")
colnames(resultsMat) <- c(names(Afile)[-1],
                          sapply(names(Afile)[-1],paste,"(sd)"))

resultsMat
