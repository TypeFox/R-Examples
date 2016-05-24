####**********************************************************************
####**********************************************************************
####
####  SPIKE AND SLAB 1.1.5
####
####  Copyright 2013, Cleveland Clinic Foundation
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  ----------------------------------------------------------------
####  Written and Developed by:
####  ----------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Director of Statistical Methodology
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    J. Sunil Rao, Ph. D.
####    Professor and Director of the Division of Biostatistics, 
####    Department of Epidemiology & Public Health
####    Clinical Research Bldg, R-669
####    1120 NW 14th Street, Room 1056
####    Miami, FL 33136
####    email:  rao.jsunil@gmail.com
####    URL:    http://biostat.med.miami.edu/people/primary-faculty/sunil-rao
####  ----------------------------------------------------------------
####  Maintained by:
####    Udaya B. Kogalur, Ph.D.
####    Adjunct Staff
####    Dept of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####    
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  ubk@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************

cv.spikeslab <- function(
 x=NULL,              #x matrix
 y=NULL,              #y response
 K=10,                #K-fold
 plot.it=TRUE,        #plot cv
 n.iter1=500,         #no. burn-in samples
 n.iter2=500,         #no. Gibbs sampled values (following burn-in)
 mse=TRUE,            #mse estimate (TRUE --> ridge/forest estimate)
 bigp.smalln=FALSE,   #used for p>>n 
 bigp.smalln.factor=1,#p>>n factor (relative to n) used in filtering variables
 screen=(bigp.smalln),#filter variables?
 r.effects=NULL,      #used for grouping variables
 max.var=500,         #max no. vars allowed in final model
 center=TRUE,         #center data (FALSE --> used for array data: !!CAUTION!!)
 intercept=TRUE,      #include an intercept?
 fast=TRUE,           #update beta in blocks (for bigp.small n, controls screening)
 beta.blocks=5,       #no. of beta blocks in beta Gibbs update (fast=TRUE)
 verbose=TRUE,        #verbose details
 save.all=TRUE,       #save spikeslab object for each fold?
 ntree=300,           #number RF trees
 seed=NULL,           #seed
  ...)
{

#tidy up x
x <- as.matrix(x)
if (length(unique(colnames(x))) != ncol(x)) {
    colnames(x) <- paste("x.", 1:ncol(x), sep = "")
}
  
# define the folds
# last fold is the full data and corresponds to the "primary object"
all.folds <-  split(sample(1:nrow(x)), rep(1:K, length = nrow(x)))
K <- length(all.folds)
all.folds[[K+1]] <- nrow(x) + 1

#core cv function
eval.fold <- function(k, ...) {
  if (verbose) {
    if (k <= K) {
      cat("\t K-fold:", k, "\n")
    }
    else {
      cat("\t final analysis (full-data)\n")
    }
  }
  omit <- all.folds[[k]]
  obj <- spikeslab(x = as.matrix(x[-omit,, drop = FALSE]), y = y[-omit],
                     n.iter1 = n.iter1, n.iter2 = n.iter2,
                     mse = mse, bigp.smalln = bigp.smalln,
                     bigp.smalln.factor = bigp.smalln.factor,
                     screen = screen, r.effects = r.effects,
                     max.var = max.var, center = center, intercept = intercept,
                     fast = fast, beta.blocks = beta.blocks, verbose = FALSE,
                     ntree = ntree, seed = seed)
  if (k <= K) {
    #test-set prediction
    pred.obj <- predict(obj, as.matrix(x[omit,, drop = FALSE])) 
    yhat.k <- pred.obj$yhat.gnet
    yhat.path.k <- rbind(pred.obj$yhat.gnet.path)
    gnet.path.k <- obj$gnet.path$path
    #lars should only return steps 0, ..., p; yet there seems to be ties
    #apply an ad-hoc beta-breaker here
    model.size.k <- apply(gnet.path.k, 1,
             function(sbeta){sum(abs(sbeta) > .Machine$double.eps, na.rm = TRUE)})
    beta.break.k <- which(!duplicated(model.size.k))
    if (length(beta.break.k) > 0) {
      gnet.path.k <- as.matrix(gnet.path.k[beta.break.k, ])
      yhat.path.k <- as.matrix(yhat.path.k[, beta.break.k])
      #test-set error
      cv.k <- mean((y[omit] - yhat.k)^2, na.rm = TRUE)
      cv.path.k <- colMeans((y[omit] - yhat.path.k)^2, na.rm = TRUE)
      gnet.k <- gnet.path.k[which(cv.path.k == min(cv.path.k))[1], ]
    }
    else {
      cv.path.k <- cv.k <- mean((y[omit] - mean(y[-omit]))^2, na.rm = TRUE)
      gnet.k <- rep(0, length(obj$names))
    }
    return(list(
              obj=(if (save.all) obj else NULL),
              model.size.k = model.size.k,
              cv.k = cv.k,
              cv.path.k = cv.path.k,
              gnet.k = gnet.k,
              bma.k = obj$bma, 
              names = obj$names))
  }
  else {
    return(list(obj = obj))
  }
}

eval.fold.obj <- mclapply(1:(K+1), eval.fold, ...)

## extract the primary object
primary.obj <- eval.fold.obj[[K+1]]$obj
eval.fold.obj[[K+1]] <- NULL
all.folds[[K+1]] <- NULL

## parse the eval.fold.obj
varnames <- primary.obj$names
p <- length(varnames)
cv <- max.model.size <- rep(NA, K)
cv.path <- model.size <- list(length = K)
bma.path <- gnet.path <- stability <- matrix(0, K, p)
cv.spikeslab.obj <- gnet.model <- vector("list", K)
cv.plot.path <-  matrix(NA, p + 1, K)
for (k in 1:K) {
  #extract the cv obj
  if (save.all) cv.spikeslab.obj[[k]] <- eval.fold.obj[[k]]$obj
  #extract cv, cv.path, gnet, gnet.path, stability, model size
  cv[k] <- eval.fold.obj[[k]]$cv.k
  cv.path[[k]] <- eval.fold.obj[[k]]$cv.path.k
  gnet.k.pt <- (abs(eval.fold.obj[[k]]$gnet.k) >  .Machine$double.eps)
  bma.path[k, is.element(varnames, eval.fold.obj[[k]]$names)] <- eval.fold.obj[[k]]$bma.k
  gnet.path[k, is.element(varnames, eval.fold.obj[[k]]$names)] <- eval.fold.obj[[k]]$gnet.k
  gnet.model[[k]] <- match(eval.fold.obj[[k]]$names[gnet.k.pt], varnames)
  stability[k, is.element(varnames, eval.fold.obj[[k]]$names)] <- 1 * gnet.k.pt
  model.size[[k]] <- eval.fold.obj[[k]]$model.size.k
  max.model.size[k] <- max(model.size[[k]], na.rm = TRUE) 
  #convert mse into matrix format more conducive for plotting/printing
  cv.plot.path[1:length(cv.path[[k]]), k] <- cv.path[[k]]
}
if (!save.all) cv.spikeslab.obj <- NULL
cv.plot.mean <- apply(cv.plot.path, 1, mean, na.rm = TRUE)
cv.plot.se <- apply(cv.plot.path, 1, SD)/sqrt(K)

## parse the primary obj
## apply an ad-hoc beta-breaker for the full path
## in extreme cases "pad" the full path with NA's
## determine the gnet as the cv optimized one
## extract other primary objects to be passed outside the wrapper
primary.gnet.path <- primary.obj$gnet.path$path
primary.model.size <- apply(primary.gnet.path, 1,
             function(sbeta){sum(abs(sbeta) > .Machine$double.eps, na.rm = TRUE)})
beta.break <- which(!duplicated(primary.model.size))
if (length(beta.break) > 0) {
    primary.gnet.path <- as.matrix(primary.gnet.path[beta.break, ])
    primary.model.size <- apply(primary.gnet.path, 1,
             function(sbeta){sum(abs(sbeta) > .Machine$double.eps, na.rm = TRUE)})
}
if (nrow(primary.gnet.path) < (p + 1)) {
  primary.gnet.path <- rbind(primary.gnet.path, matrix(NA, (p + 1) - nrow(primary.gnet.path), p))
}

primary.gnet.cv <- cv.plot.mean
#in rare cases the full data gnet may have a smaller path than the cv-gnet paths
which.best.path <- min(which(primary.gnet.cv == min(primary.gnet.cv, na.rm = TRUE))[1],
                   1 + sum(abs(primary.obj$bma) > .Machine$double.eps))
primary.gnet.scale <- primary.gnet.path[which.best.path, ]
primary.gnet <- primary.gnet.scale * primary.obj$x.scale
primary.gnet.path <- list(path = primary.gnet.path, cv = primary.gnet.cv, model.size = model.size)
primary.obj$gnet <- primary.gnet
primary.obj$gnet.scale <- primary.gnet.scale


## plot it
if (plot.it) {
  matplot(0:p, cv.plot.path, type = c("l", "n")[1 + 1 * (K > 20)], lty = 3, col = "gray", 
         xlim = range(c(0, max.model.size), na.rm = TRUE),
         ylim = range(c(cv.plot.path, cv.plot.mean + cv.plot.se, cv.plot.mean - cv.plot.se), na.rm = TRUE),
         xlab="Model Size", ylab="Cross-Validated MSE")
  lines(0:p, cv.plot.mean, lty = 1, lwd = 2, col = 4)
  error.bars(0:p, cv.plot.mean + cv.plot.se, cv.plot.mean - cv.plot.se, width = 0.0025, col = 2)
}

# stability analysis; make it pretty for the return
tally.stability <- cbind(primary.obj$bma,
                         apply(bma.path, 2, mean, na.rm = TRUE),
                         primary.obj$gnet,
                         apply(gnet.path, 2, mean, na.rm = TRUE) * primary.obj$x.scale,
                         100 * apply(stability, 2, mean, na.rm = TRUE))
colnames(tally.stability) <- c("bma", "bma.cv", "gnet", "gnet.cv", "stability")
rownames(tally.stability) <- varnames
tally.stability <- tally.stability[order(tally.stability[, 5], abs(tally.stability[, 1]),
                          decreasing = TRUE),, drop = FALSE]


### --------------------------------------------------------------
###	Terminal Output
###     Save as list
### --------------------------------------------------------------	

# functions for determining model sizes
get.model.size <- function(mn, se) {
  ms.upper <- min(which(mn == min(mn, na.rm = TRUE)))
  mn.upper <- mn[ms.upper]
  ms.range <- which((mn + se >= mn.upper) & (mn - se <= mn.upper))
  ms.lower <- min(ms.range, na.rm = TRUE)
  if (is.infinite(ms.lower)) ms.lower <- ms.upper
  ms.all <- unique(c(ms.lower, ms.upper))
  if (length(ms.all) == 1) return(ms.all)
  paste("[", ms.lower, ",", ms.upper, "]", sep = "")
}

verbose.list <- list(
  c("cross-validation"),
  c(bigp.smalln),
  c(screen),
  c(fast),
  c(nrow(x)),
  c(ncol(x)),
  c(n.iter1),
  c(n.iter2),
  c(K),
  c(paste(round(mean(cv, na.rm = TRUE), 3) , "+/-",
    round(sd(cv, na.rm = TRUE)/sqrt(K), 3))),
  c(get.model.size(cv.plot.mean, cv.plot.se)),
  c(sum(abs(primary.gnet) > .Machine$double.eps, na.rm = TRUE))
)

if (verbose) {
     cat("-------------------------------------------------------------------","\n")
    cat("Variable selection method     :",verbose.list[[1]],"\n")
    cat("Big p small n                 :",verbose.list[[2]],"\n")
    cat("Screen variables              :",verbose.list[[3]],"\n")
    cat("Fast processing               :",verbose.list[[4]],"\n")
    cat("Sample size                   :",verbose.list[[5]],"\n")
    cat("No. predictors                :",verbose.list[[6]],"\n")
    cat("No. burn-in values            :",verbose.list[[7]],"\n")
    cat("No. sampled values            :",verbose.list[[8]],"\n")
    cat("K-fold                        :",verbose.list[[9]],"\n")
    cat("CV mean-squared error         :",verbose.list[[10]],"\n")
    cat("Model size                    :",verbose.list[[11]],"\n")
    cat("\n\nTop variables in terms of stability:\n")
    print(head(round(tally.stability, 3), verbose.list[[12]]))
    cat("-------------------------------------------------------------------","\n")
}

#return the goodies
object <- list(
               spikeslab.obj = primary.obj,          #spikeslab obj from full data
               cv.spikeslab.obj=cv.spikeslab.obj,    #cv spikeslab obj for each fold
               cv.folds=all.folds,                   #cv folds
               cv = cv,                              #mean cv for each fold
               cv.path = cv.plot.path,               #cv path (p+1) x K
               stability = tally.stability,          #stability values
               bma = primary.obj$bma,                #bma coefficients, for scaled x   
               bma.scale = primary.obj$bma.scale,    #bma coefficiebnts, for original x   
               gnet = primary.obj$gnet,              #cv-optimized gnet, for scaled x
               gnet.scale = primary.obj$gnet.scale,  #cv-optimized gnet, for original x   
               gnet.model = gnet.model,              #models used to determine gnet
               gnet.path = primary.gnet.path,        #gnet path (full data)
               gnet.obj = primary.obj$gnet.obj,      #gnet object (full data)
               gnet.obj.vars =  primary.obj$gnet.obj.vars,#variables used to define gnet
               verbose = verbose.list)               #verbose (for print)
class(object) <- c("spikeslab", "cv")
invisible(object)

}
