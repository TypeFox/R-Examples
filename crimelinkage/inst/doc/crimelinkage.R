## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = TRUE, echo=TRUE,comment = "#>")
library(crimelinkage)
CACHE = FALSE # Set TRUE for development, then set to FALSE for release

## ------------------------------------------------------------------------
#- If the crimelinkage package is not installed, type: install.packages{"crimelinkage"}
library(crimelinkage)

## ------------------------------------------------------------------------
data(crimes)
head(crimes)

## ------------------------------------------------------------------------
data(offenders)
head(offenders)

## ------------------------------------------------------------------------
cs = getCrimeSeries("O:40",offenders)
cs

## ------------------------------------------------------------------------
getCrimes("O:40",crimes,offenders)

## ------------------------------------------------------------------------
getCriminals(cs$crimeID,offenders)

## ------------------------------------------------------------------------
getCriminals("C:78",offenders)

## ----seriesData----------------------------------------------------------
seriesData = makeSeriesData(crimedata=crimes,offenderTable=offenders)

## ------------------------------------------------------------------------
nCrimes = table(seriesData$offenderID)  # length of each crime series
table(nCrimes)                          # distribution of crime series length
mean(nCrimes>1)                         # proportion of offenders with multiple crimes

## ------------------------------------------------------------------------
nCO = table(seriesData$crimeID) # number of co-offenders per crime
table(nCO)                      # distribution of number of co-offenders
mean(nCO>1)                     # proportion of crimes with multiple co-offenders

## ---- get.indices,cache=CACHE--------------------------------------------
set.seed(1)         # set random seed for replication
allPairs = makePairs(seriesData,thres=365,m=40)

## ----make.linkage.data, cache=CACHE, results='hide'----------------------
varlist = list( spatial = c("X", "Y"), 
                temporal = c("DT.FROM","DT.TO"), 
                categorical = c("MO1",  "MO2", "MO3"))    # crime variables list
X = compareCrimes(allPairs,crimedata=crimes,varlist=varlist,binary=TRUE) # Evidence data
Y = ifelse(allPairs$type=='linked',1,0)      # Linkage indicator. 1=linkage, 0=unlinked

## ------------------------------------------------------------------------
head(X)
table(Y)         

## ------------------------------------------------------------------------
set.seed(3)                                        # set random seed for replication
train = sample(c(TRUE,FALSE),nrow(X),replace=TRUE,prob=c(.7,.3))  # assign pairs to training set
test = !train
D.train = data.frame(X[train,],Y=Y[train])          # training data

## ------------------------------------------------------------------------
vars = c("spatial","temporal","tod","dow","MO1","MO2","MO3") 
fmla.all = as.formula(paste("Y ~ ", paste(vars, collapse= "+")))
fmla.all

## ----logistic.regression, message=FALSE,warning=FALSE--------------------
fit.logistic = glm(fmla.all,data=D.train,family=binomial,weights=weight) 
fit.logisticstep = step(fit.logistic,trace=0)
summary(fit.logisticstep)

## ---- results='hold',fig.keep='high',fig.show='hold',fig.width=7,fig.height=4, out.width="100%"----
NB = naiveBayes(fmla.all,data=D.train,weights=weight,df=10,nbins=15,partition='quantile')

#- Component Plots
plot(NB,ylim=c(-2,2))       # ylim sets the limits of the y-axis

## ----gbm, cache=FALSE, results='hide',fig.keep='none',fig.show='hold',fig.width=6,fig.height=6,message=FALSE,warning=FALSE----
library(gbm) # if not installed type: install.packages("gbm")
set.seed(4)
fit.gbm = gbm(fmla.all,distribution="bernoulli",data=D.train,weights=weight,
              shrinkage=.01,n.trees=500,interaction.depth=3)
nt.opt = fit.gbm$n.trees     # see gbm.perf() for better options

#- Relative influence and plot
print(relative.influence(fit.gbm,n.trees=nt.opt))
par(mfrow=c(2,4))
for(i in 1:7) plot(fit.gbm,i)

## ------------------------------------------------------------------------
D.test = data.frame(Y=Y[test],X[test,])
X.test = D.test[,vars]
Y.test = D.test[,"Y"]

## ------------------------------------------------------------------------
#- Predict the Bayes factor for each model
BF = data.frame(
  logisticstep = predict(fit.logisticstep,newdata=X.test,type='link'),
  nb = predict(NB,newdata=X.test),
  gbm = predict(fit.gbm,newdata=X.test,n.trees=nt.opt)
  )

#- Evaluation via ROC like metrics
roc = apply(BF,2,function(x) getROC(x,Y.test))

## ----fig.show='hold',fig.width=10, out.width="100%"----------------------
with(roc[[1]], plot(Total,TPR,typ='l',xlim=c(0,1000),ylim=c(0,1),las=1,
                    xlab='Total Cases Examined',ylab='True Positive Rate'))
grid()
for(i in 2:length(roc)){
  with(roc[[i]], lines(Total,TPR,col=i))
}
legend('topleft',names(roc),col=1:length(roc),lty=1, cex=.75,bty="n")
title("Proportion of overall linked crimes captured")

## ----fig.show='hold',fig.width=10, out.width="100%"----------------------
with(roc[[1]], plot(Total,PPV,typ='l',xlim=c(0,1000),ylim=c(0,1),las=1,
                    xlab='Total Cases Examined',ylab='Precision'))
grid()
for(i in 2:length(roc)){
  with(roc[[i]], lines(Total,PPV,col=i))
}
legend('topright',names(roc),col=1:length(roc),lty=1, cex=.75,bty="n")
title("Proportion of identified crimes that are actually linked")

