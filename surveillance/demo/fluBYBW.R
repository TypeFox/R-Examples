################################################################################
### Demo of hhh4() modelling of influenza in Southern Germany - data("fluBYBW")
### RUNNING THE WHOLE SCRIPT TAKES ~20 MINUTES!
###
### Copyright (C) 2009-2012 Michaela Paul, 2012-2013,2016 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

library("surveillance")

## Weekly counts of influenza in 140 districts of Bavaria and Baden-Wuerttemberg
data("fluBYBW")


##################################################
# Fit the models from the Paul & Held (2011) paper
##################################################

## generate formula for temporal and seasonal trends
f.end <- addSeason2formula(f = ~ -1 + ri(type="iid", corr="all") + 
                               I((t-208)/100), S=3, period=52)

## settings for the optimizer
opt <- list(stop = list(tol=1e-5, niter=200),
            regression = list(method="nlminb"),
            variance = list(method="nlminb"))

## models
# A0
cntrl_A0 <- list(ar = list(f = ~ -1),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose = 1)
summary(res_A0 <- hhh4(fluBYBW,cntrl_A0))

# B0
cntrl_B0 <- list(ar = list(f = ~ 1),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1)
res_B0 <- hhh4(fluBYBW,cntrl_B0)               
 

# C0
cntrl_C0 <- list(ar = list(f = ~ -1 + ri(type="iid", corr="all")),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1)
res_C0 <- hhh4(fluBYBW,cntrl_C0)               


#A1

# weight matrix w_ji = 1/(No. neighbors of j) if j ~ i, and 0 otherwise
wji <- neighbourhood(fluBYBW)/rowSums(neighbourhood(fluBYBW))

cntrl_A1 <- list(ar = list(f = ~ -1),
                 ne = list(f = ~ 1, weights = wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1)
res_A1 <- hhh4(fluBYBW,cntrl_A1)               


# B1
cntrl_B1 <- list(ar = list(f = ~ 1),
                 ne = list(f = ~ 1, weights = wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1)
res_B1 <- hhh4(fluBYBW,cntrl_B1)               


# C1
cntrl_C1 <- list(ar = list(f = ~ -1 + ri(type="iid", corr="all")),
                 ne = list(f = ~ 1, weights = wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1)
res_C1 <- hhh4(fluBYBW,cntrl_C1)               


#A2
cntrl_A2 <- list(ar = list(f = ~ -1),
                 ne = list(f = ~ -1 + ri(type="iid",corr="all"), weights=wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1)
res_A2 <- hhh4(fluBYBW,cntrl_A2)               


# B2
cntrl_B2 <- list(ar = list(f = ~ 1),
                 ne = list(f = ~ -1 + ri(type="iid",corr="all"), weights =wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1)
res_B2 <- hhh4(fluBYBW,cntrl_B2)               

# C2
cntrl_C2 <- list(ar = list(f = ~ -1 + ri(type="iid", corr="all")),
                 ne = list(f = ~ -1 + ri(type="iid",corr="all"), weights =wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1,
                 start=list(fixed=fixef(res_B0),random=c(rep(0,140),
                         ranef(res_B0)), sd.corr=c(-.5,res_B0$Sigma.orig,0)))
res_C2 <- hhh4(fluBYBW,cntrl_C2)               


# D
cntrl_D <- list(ar = list(f = ~ 1),
                ne = list(f = ~ -1 + ri(type="iid"), weights = wji),
                end = list(f =addSeason2formula(f = ~ -1 + ri(type="car") + 
                                             I((t-208)/100), S=3, period=52), 
                          offset = population(fluBYBW)),
                family = "NegBin1", optimizer = opt, verbose=1)
res_D <- hhh4(fluBYBW,cntrl_D)               


######################################################################
# Compare the predictive performance of the models by computing
# one-step-ahead predictions to be assessed by proper scoring rules
######################################################################

## do 1-step ahead predictions for the last two years

tp <- nrow(fluBYBW)-2*52

## for this demo: only calculate pseudo-predictions based on the final fit
## to avoid the time-consuming sequential refitting at each step.
TYPE <- "final"
## use "rolling" for true one-step-ahead predictions => TAKES ~8 HOURS!

val_A0 <- oneStepAhead(res_A0, tp=tp, type=TYPE)
val_B0 <- oneStepAhead(res_B0, tp=tp, type=TYPE)
val_C0 <- oneStepAhead(res_C0, tp=tp, type=TYPE)

val_A1 <- oneStepAhead(res_A1, tp=tp, type=TYPE)
val_B1 <- oneStepAhead(res_B1, tp=tp, type=TYPE)
val_C1 <- oneStepAhead(res_C1, tp=tp, type=TYPE)

val_A2 <- oneStepAhead(res_A2, tp=tp, type=TYPE)
val_B2 <- oneStepAhead(res_B2, tp=tp, type=TYPE)
val_C2 <- oneStepAhead(res_C2, tp=tp, type=TYPE)

val_D <- oneStepAhead(res_D, tp=tp, type=TYPE)


## compute scores

vals <- ls(pattern="val_")
nam <- substring(vals,first=5,last=6)

whichScores <- c("logs", "rps", "ses")
scores_i <- vector(mode="list", length=length(vals))
meanScores <- NULL
for(i in seq_along(vals)){
  sc <- scores(get(vals[i]), which=whichScores, individual=TRUE)
  scores_i[[i]] <- sc
  meanScores <- rbind(meanScores,colMeans(sc, dims=2))
}

names(scores_i) <- nam
rownames(meanScores) <- nam

print(meanScores)
## Note that the above use of "final" fitted values instead of "rolling"
## one-step-ahead predictions leads to different mean scores than reported
## in Paul & Held (2011, Table IV).


## comparison with best model B2 

compareWithBest <- function(best, whichModels, nPermut=9999, seed=1234){
  set.seed(seed)
  pVals <- NULL
  for(score in seq_along(whichScores)){
    p <- c()
    for(model in whichModels){
      if(model==best) p <- c(p,NA)
      else p <- c(p,permutationTest(scores_i[[model]][,,score],scores_i[[best]][,,score],
     plot=TRUE,nPermutation=nPermut, verbose=TRUE)$pVal.permut)
    }  
    pVals <- cbind(pVals,p)
  }
  return(pVals)
}

pVals_flu <- compareWithBest(best=9, whichModels=1:10,
                             nPermut=999, # reduced for this demo
                             seed=2059710987)
rownames(pVals_flu) <- nam
colnames(pVals_flu) <- whichScores

print(pVals_flu)
