## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = TRUE, echo=TRUE,comment = "#>",
                      cache=FALSE)
library(crimelinkage)

## ----loaddata------------------------------------------------------------
#-- Load the package and get the example crime data
library(crimelinkage)
data(crimes)               # some example crime incident data
data(offenders)            # some example crime offender data
seriesData = makeSeriesData(crimedata=crimes,offenderTable=offenders)

## ----caselinkage, results='hide'-----------------------------------------
#-- Make Crime Pairs for Case Linkage
set.seed(1)         # set random seed for replication
allPairs = makePairs(seriesData,thres=365,m=40)

#-- Make Evidence Variables for Case Linkage
varlist = list( spatial = c("X", "Y"), 
                temporal = c("DT.FROM","DT.TO"), 
                categorical = c("MO1",  "MO2", "MO3"))    # crime variables list
X = compareCrimes(allPairs,crimedata=crimes,varlist=varlist,binary=TRUE) # Evidence data
Y = ifelse(allPairs$type=='linked',1,0)      # Linkage indicator. 1=linkage, 0=unlinked


#-- Get Training Data
set.seed(3)                                        # set random seed for replication
train = sample(c(TRUE,FALSE),nrow(X),replace=TRUE,prob=c(.7,.3))  # assign pairs to training set
test = !train
D.train = data.frame(X[train,],Y=Y[train])          # training data

#-- Fit naive Bayes model and make estimateBF() function 
vars = c("spatial","temporal","tod","dow","MO1","MO2","MO3") 
fmla.all = as.formula(paste("Y ~ ", paste(vars, collapse= "+")))
NB = naiveBayes(fmla.all,data=D.train,weights=weight,df=10,nbins=15,partition='quantile')

estimateBF <- function(X){           # estimateBF() returns the estimated log Bayes factor
  predict(NB,newdata=X)
}

## ----hierclustering, results='hide', out.width="100%",fig.width=12,fig.height=7----
#-- Get unsolved crimes
unsolved = subset(crimes, !crimeID %in% seriesData$crimeID)

#-- Run agglomerative hierarchical crime clustering
tree = crimeClust_hier(unsolved,varlist,estimateBF,linkage='average', binary=TRUE)

#-- Plot results in dendrogram using plot_hcc()
plot_hcc(tree,yticks=seq(-2,6,by=2),type="triangle",hang=.05,main="Average Linkage") 

## ------------------------------------------------------------------------
#-- Examine crimes C:431 and C:460 
subset(crimes,crimeID %in% c('C:431','C:460'))

## ------------------------------------------------------------------------
#-- Find path info for crime C:429
cp = clusterPath('C:429',tree)
cp[cp$logBF>0,]                 # only return path for scores > 0

## ------------------------------------------------------------------------
solved = subset(crimes, crimeID %in% seriesData$crimeID)
unsolved = subset(crimes, !crimeID %in% seriesData$crimeID)

## ------------------------------------------------------------------------
crime = unsolved[2,]             # use the 2nd unsolved crime C:392
crime
results = seriesID(crime,solved,seriesData,varlist,estimateBF)
head(results$score)

## ------------------------------------------------------------------------
subset(results$groups,group=='12')      # most similar crime series
subset(results$groups,group=='154')     # 2nd most similar series
subset(results$groups,group=='9')       # a series with multiple crimes

## ------------------------------------------------------------------------
crime4 = unsolved[4,]             # use the 4th unsolved crime
results4 = seriesID(crime4,solved,seriesData,varlist,estimateBF)
head(results4$score)

## ------------------------------------------------------------------------
#- using crime C:394 (the 4th unsolved crime)
pairs = data.frame(i1=unsolved$crimeID[4],i2=unique(unsolved$crimeID[-4]))  
X = compareCrimes(pairs,unsolved,varlist,binary=TRUE)     # Evidence data
score = data.frame(pairs,logBF=estimateBF(X))  
head(score[order(-score$logBF),])

## ------------------------------------------------------------------------
C429 = which(unsolved$crimeID %in% 'C:429')       # now use crime C:429
pairs = data.frame(i1=unsolved$crimeID[C429],i2=unique(unsolved$crimeID[-C429]))  
X = compareCrimes(pairs,unsolved,varlist,binary=TRUE)     # Evidence data
score = data.frame(pairs,logBF=estimateBF(X))  
head(score[order(-score$logBF),])             

#-- results from hierarchical clustering
cp = clusterPath('C:429',tree)
cp[cp$logBF>0,]   

## ----MCMC,eval=FALSE,message=FALSE,results='hide',fig.keep='none'--------
#  #-- Make the crime group labels for each crime (NA for unsolved crimes)
#  seriesData$CG = makeGroups(seriesData,method=2)      # method=2 uses unique co-offenders
#  group_info = unique(seriesData[,c('crimeID','CG')])  # extract the group info
#  A = merge(crimes,group_info,by="crimeID",all.x=TRUE) # add group info to crimes
#  A = A[order(A$CG),]                                  # order by crime group
#  
#  #-- Run MCMC
#  fit = crimeClust_bayes(A$CG, spatial=A[,c('X','Y')],t1=A$DT.FROM,t2=A$DT.TO,
#                         Xcat=A[,c("MO1","MO2","MO3")],maxcriminals=1000,
#                         iters=3000,burn=1000,update=100,seed=5)
#  
#  #-- Extract pairwise probabilities
#  pp = fit$p.equal    # probability that crime i is linked to crime j
#  diag(pp) = NA

## ----echo=FALSE,eval=FALSE-----------------------------------------------
#  #-- Save the results that take a long time to run
#  save(A,fit,pp,file="vignettes/MCMC-results.RData")

## ----echo=FALSE,eval=TRUE------------------------------------------------
load("MCMC-results.RData")

## ----message=FALSE,out.width="95%",fig.width=10,fig.height=8-------------
library(fields) # if not installed, type: install.packages("fields")

#-- Get index of unsolved crimes
ind.unsolved = which(is.na(A$CG))          # index of unsolved crimes
n = nrow(A)                                # number of crimes

#-- Image plot of linkage probabilities
fields::image.plot(1:n,ind.unsolved,pp[1:n,ind.unsolved],
           xlab="Crime",ylab="Unsolved Crime",
           main="Probability crimes are linked")

## ----fig.height=6,fig.width=8,out.width="70%"----------------------------
#-- Find strongest linkages
unsolved.probs = apply(pp[ind.unsolved,],1,max,na.rm=TRUE)  # maximum probability
plot(ind.unsolved,unsolved.probs,xlab="unsolved crime",ylab='maximum probability of linkage')
abline(h=0.25)
ind = ind.unsolved[unsolved.probs > 0.25]
investigate = as.character(A$crimeID[ind])       # crimeIDs for crimes with strongest linkage
investigate

## ------------------------------------------------------------------------
bp = bayesProb(pp[A$crimeID %in% "C:417"])
bp$crimeID = A$crimeID[bp$index]
bp$CG = A$CG[bp$index]
head(bp)

