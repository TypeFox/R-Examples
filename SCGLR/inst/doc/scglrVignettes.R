## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
  concordance=FALSE, fig.path='scglr-',tidy=FALSE,size="small"
)

## ----echo=FALSE-------------------------------------------
options(width=60,prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)

## ----echo=FALSE-------------------------------------------
suppressPackageStartupMessages(library("SCGLR"))

## ----eval=FALSE-------------------------------------------
#  results.scglr <- scglr(formula,data,family,K,size,offset,subset,
#                         na.action,crit,method)

## ----results='hide',echo=FALSE----------------------------
ny <- paste("y",1:2,sep="")
nx <-paste("x",1:5,sep="")
nz <- paste("z",1:3,sep="")

## ----eval=TRUE--------------------------------------------
myformula <- multivariateFormula(ny,nx,nz)
myformula

## ----eval=FALSE-------------------------------------------
#  scglrCrossVal(formula,data,family,K,nfolds,types,size,
#                  offset,subset,na.action,crit,method)

## ---------------------------------------------------------
data("genus")
dim(genus)

## ---------------------------------------------------------
names(genus)

## ----eval=TRUE--------------------------------------------
ny <- names(genus)[1:27]
sx <- which(names(genus) %in% c("geology","surface"))
nx <- names(genus)[-c(1:27,sx)] 
family <- rep("poisson",length(ny))
formula <- multivariateFormula(ny,c(nx,"I(lon*lat)"),"geology")
formula
offset <- genus$surface

## ----tempo,eval=TRUE--------------------------------------
criterion <- t(apply(genus.cv,1,function(x) x/mean(x)))
criterion.mean <- apply(criterion,2,mean) 
K.cv <- which.min(criterion.mean)-1

## ----plotCv, fig.align='center', fig.pos='!ht', fig.cap='Mean Squared Prediction Error (MSPE) as a function of the number of components.'----
plot(0:K,criterion.mean, type="l",
     xlab="K, number of components", ylab="Criterion (MSPE)")
Axis(side=1,at=0:K)
abline(v=K.cv,col=2)

## ----eval=TRUE--------------------------------------------
print(genus.scglr)

## ----barplotScglr, fig.align='center', fig.pos='!ht', fig.cap='Barplot of inertia per component'----
barplot(genus.scglr)

## ----samplePlots, fig.align='center',fig.pos='!ht',fig.cap='Two sample plots',fig.subcap=c('Simple correlation plot','Correlation plot with linear predictors and covariates passing a threshold of $0.8$'),fig.show='hold',out.width='0.49\\linewidth'----
plot(genus.scglr)
plot(genus.scglr, threshold=0.8, predictors=TRUE)

## ----pairsScglr,fig.align='center', fig.pos='!ht', fig.cap='Correlation plots on planes spanned by components'----
pairs(genus.scglr,ncol=2,label.size=0.5) 

