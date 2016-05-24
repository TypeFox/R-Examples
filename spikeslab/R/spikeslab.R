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

spikeslab <- function(
 formula,             #R formula 
 data=NULL,           #data frame
 x=NULL,              #NULL formula --> X,Y data
 y=NULL,              #NULL formula --> X,Y data
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
 verbose=FALSE,       #verbose details?
 ntree=300,           #number RF trees
 seed=NULL,           #seed
  ...)
{


### --------------------------------------------------------------
###
###  Global Parameters
###  Do not touch unless you know what you are doing !!!
###  
### --------------------------------------------------------------

eps                <-  .Machine$double.eps          #effective zero
ridge.tolerance    <- 0.01                          #lower bound ridge
min.df.mse         <- 25                            #min df for ridge mse estimate
prior              <- c("cb","mcb","lasso","bimodal")[1]
                                                    #cb=continuous bimodal
                                                    #mcb=mixed continuous bimodal
                                                    #lasso=double-exponential
                                                    #bimodal=two-atom mixture

###additional priors not currently implemented
###if (prior == "lasso" | prior == "mcb") library(statmod)

### --------------------------------------------------------------
###  Read data & preprocess
### --------------------------------------------------------------

### ---------------------------------
#  Set seed if appropriate
#  Process via non-formula or formula mode
#  Read data + set dimensions
#  Fit intercept by centering Y (unless user requests no intercept)
if (is.null(seed)) seed <- -1.0*abs(round(rnorm(1, sd=1e5)))
mf <- match.call(expand.dots = FALSE)
if (is.na(match("formula",names(mf)))) {
  if (is.null(y) | is.null(x)) stop("x and or y is missing")
  Y.org <- y
  X.org <- x
  beta.names <- colnames(X.org)
  if (length(unique(beta.names)) != ncol(X.org)) {
    colnames(X.org) <- beta.names <- paste("x.", 1:ncol(X.org), sep = "")
  }
  mt <- NULL
}
else {
 m <- match(c("formula", "data"), names(mf), 0)
 mf <- mf[c(1, m)]
 mf$drop.unused.levels <- T
 mf[[1]] <- as.name("model.frame")
 mf <- eval(mf, parent.frame())
 mt <- attr(mf, "terms")
 attr(mt, "intercept") <- 0
 Y.org <- model.response(mf, "numeric")
 X.org <- model.matrix(mt, mf)
 beta.names <- unlist(dimnames(X.org)[2])
}
#double check for intercepts (shouldn't have made it this far)
if(intercept & any(beta.names=="(Intercept)")) {
 int.pt <- (beta.names == "(Intercept)")
 X.org <- X.org[, !int.pt]
 beta.names <- beta.names[!int.pt]
}
X.org <- as.matrix(X.org)
n.data <- nrow(X.org)
n.cov <- length(beta.names)
rm(mf)

### ---------------------------------
### check coherence of n.iter
n.iter1 <- round(max(1, n.iter1))
n.iter2 <- round(max(1, n.iter2))

### ---------------------------------
### check coherence of max.var
max.var <- max(1, min(max.var, n.cov), na.rm = TRUE)

### ---------------------------------
# clean up beta.names
beta.names <- sub('`',"",beta.names)
beta.names <- sub('`',"",beta.names)

### ---------------------------------
# random effects
# non-null fixed effects handled specially
turn.sigma.on <- FALSE
f.eff <- r.eff <- NULL
if (!is.null(r.effects)) {
  r.eff <- vector("list", length(r.effects))
  for (k in 1:length(r.eff)) {
    match.k <- match(r.effects[[k]], beta.names)
    match.k <- match.k[!is.na(match.k)]
    if (length(match.k) == 0) {
      stop("r.effects error: variable grouping contains variable names not found in the data") 
    }
    if (length(unique(match.k)) < length(match.k)) {
        stop("r.effects error: variable names are duplicated within a variable grouping")
    }
    r.eff[[k]] <- match.k
  }
  r.eff.intersect <- NULL
  if (length(r.eff) > 1) {
    for (k in 1:(length(r.eff)-1)) {
      r.eff.intersect <- c(r.eff.intersect, intersect(r.eff[[k]], r.eff[[k+1]]))
      if (length(r.eff.intersect) > 0) {
        stop("r.effects error: variable groups are not distinct") 
      }
    }
  }
  f.eff <- setdiff(1:n.cov, unlist(r.eff))
  if (length(f.eff) == 0) f.eff <- NULL
  r.eff.names <- (1:length(r.effects))
  prior <- "cb"
  turn.sigma.on <- TRUE
}

### ---------------------------------
# variable selection method no longer supported
method <- "AIC"

### ---------------------------------
# check coherence of priors
if (length(intersect(prior, c("cb","mcb","lasso","bimodal"))) !=1)
  stop("Incorrect prior choice: ", prior)

### ---------------------------------
# big p, small n details

if (bigp.smalln) {
  if (n.cov < n.data) stop("n > p: big p small n option makes no sense!!\n")
  forest <- FALSE
  mse <- FALSE
  turn.sigma.on <- TRUE
  screen <- TRUE
  prior <- "cb"
}

### ---------------------------------
#  Save sufficient statistics for original X and Y
Y.org.mean <- mean.center(Y.org, center = intercept)
Y.center <- Y.org - Y.org.mean
X.org.sd <- as.double(c(apply(X.org, 2, sd.center, center = center)))
X.org.mean <- as.double(apply(X.org, 2, mean.center, center = center))

### ---------------------------------
### Define working X matrix
X.wrk <- scale(X.org, center=X.org.mean, scale=X.org.sd)
X.wrk[, X.org.sd==0] <- 0
X.org.mean[X.org.sd==0] <- 0
X.wrk <- as.matrix(X.wrk)
row.names(X.wrk) <- colnames(X.wrk) <- NULL

### ---------------------------------
# External mse estimate 
# Reverts to ridge if df are large enough, otherwise uses RF
if (mse) {
  if (verbose) {
    if (n.cov <= n.data) {
      cat("\n", "\t pre-processing data... \n")
    }
    else {
      cat("\n", "\t pre-processing data (p is bigger than n, *consider* using 'bigp.smalln=TRUE')... \n")
    }
  }
  if ((n.data - n.cov)  > min.df.mse) {
    mse.hat <- sum((Y.center - cbind(X.wrk)%*%Ridge(Y.center, X.wrk, ridge.tolerance))^2)/(n.data-n.cov)
  }
  else {
    rf.fit <- randomForest(cbind(X.wrk), Y.center, ntree=ntree, importance=TRUE)
    rf.res <- rf.fit$y - predict(rf.fit, X.wrk)
    mse.hat <- rf.fit$mse[ntree]-mean(rf.res^2, na.rm = TRUE)
  }
  turn.sigma.on <- FALSE
}
else {
  mse.hat <- 1
  turn.sigma.on <- TRUE
}



### --------------------------------------------------------------
###
###  MAIN FUNCTION
###
### --------------------------------------------------------------

  if (verbose) cat("\t running spike and slab regression...\n")

  ### ---------------------------------
  ### Rescaled Y, XX, XY 
  ### Compute scale factor, sf  
  sf <- sqrt(n.data/mse.hat)
  Y.wrk <- Y.center*sf
  XX.wrk <- XX.multiply(X.wrk, bigp.smalln)
  sum.xy.wrk <- t(X.wrk)%*%Y.wrk

  nozap.pt <- 1:n.cov
  hyperv <- model <- NULL

  ### ---------------------------------
  ### Call core Gibbs routine:
  ### 1. Screen approach for complexity reduction
  ### 2. Straight call
  if (screen) {
    ### ---------------------------------
    ### pre-filter variables
    ###
    ### bigp.smalln=F         : ss{regular    &x-y cor}
    ### fast=T, bigp.smalln=T : ss{orthogonal &x-y cor} 
    ### fast=F, bigp.smalln=T : ss{regular    &x-y cor}
    ###
    ### core Gibbs call with noise
    gibbs <- spikeslab.GibbsCore(
               n.iter1=(if (bigp.smalln & !fast) min(400, n.iter1) else n.iter1),
               n.iter2=(if (bigp.smalln & !fast) min(400, n.iter2) else n.iter2),
               orthogonal=(bigp.smalln & fast),
               prior,
               fast, beta.blocks,
               X=X.wrk, Y=Y.wrk, XX=XX.wrk, sum.xy=sum.xy.wrk,
               seed=seed, verbose=verbose,
               bigp.smalln=bigp.smalln,
               turn.sigma.on=turn.sigma.on,
               correlation.filter=TRUE,
               r.eff=r.eff)
    complexity.vec <-  gibbs$complexity.vec
    nozap.pt <- which(abs(gibbs$b.m) > eps)
    ### make sure p < min(max.var, (factor * n)) 
    ### remove least important variables
    if (bigp.smalln) {
      max.bigp.cov <- min(max.var, bigp.smalln.factor * n.data)
    }
    else {
      max.bigp.cov <- max.var
    }
    if (length(nozap.pt) > max.bigp.cov) {
      nozap.new.pt <- order(abs(gibbs$b.m), decreasing = TRUE)
      nozap.pt <- nozap.new.pt[1:max.bigp.cov]
    }
    ### re-adjust r.eff for variable screening
    if (!is.null(r.eff) & (length(nozap.pt) > 0)) {
      n.reff <- 0
      r.eff <- vector("list", 0)
      r.eff.names <- NULL
      for (k in 1:length(r.effects)) {
        match.k <- match(r.effects[[k]], beta.names[nozap.pt])
        match.k <- match.k[!is.na(match.k)]
        if (length(match.k) > 0) {
          n.reff <- n.reff + 1
          r.eff[[n.reff]] <- match.k
          r.eff.names <- c(r.eff.names, k)
        }
      }
      if (length(r.eff) == 0) r.eff <- NULL
    }
    ### core Gibbs call with complexity reduced space
    ### for big p small n, invoke fast variable screening
    b.m <- rep(0, n.cov)
    if (length(nozap.pt) > 0) {
      if (bigp.smalln)  {
        XX.wrk.reduced <- XX.multiply(as.matrix(X.wrk[, nozap.pt]), FALSE)
      }
      else {
        XX.wrk.reduced <- as.matrix(XX.wrk[nozap.pt, nozap.pt])
      }
      gibbs <- spikeslab.GibbsCore(n.iter1, n.iter2, 
          orthogonal=FALSE, prior,
          fast=(fast | (bigp.smalln)),
          beta.blocks,
          X=as.matrix(X.wrk[, nozap.pt]), Y=Y.wrk,
          XX=XX.wrk.reduced, sum.xy=sum.xy.wrk[nozap.pt],
          seed=seed, verbose=verbose,
          bigp.smalln=FALSE,
          turn.sigma.on=turn.sigma.on,
          r.eff=r.eff)
      b.m[nozap.pt] <-  gibbs$b.m
      complexity.vec <-  gibbs$complexity.vec
      hyperv <-  lapply(1:length(gibbs$hyperv), function(i) {
        hyperv.i <- rep(Inf, n.cov)
        hyperv.i[nozap.pt] <- gibbs$hyperv[[i]]
        hyperv.i
      })
      model <-  lapply(1:length(gibbs$model), function(i) {
        nozap.pt[gibbs$model[[i]]]
      })
      if (turn.sigma.on) mse.hat <- mean(gibbs$sigma.vec, na.rm=T)
    }
    b.m[is.na(b.m)] <- 0  #remove NA's from the bma (rare)
    phat.bma <- sum(abs(b.m) > eps, na.rm = TRUE)
    penal <- rep(Inf, n.cov)
    resid <- rep(0, n.cov)
    if (length(nozap.pt) > 0) {
      resid[nozap.pt] <- (sum.xy.wrk[nozap.pt] - XX.wrk.reduced %*% b.m[nozap.pt])
      penal[nozap.pt] <- resid[nozap.pt] / b.m[nozap.pt]
    }

  }
  
  else {
    
    ### ---------------------------------
    ### Straight core Gibbs call
    gibbs <- spikeslab.GibbsCore(n.iter1, n.iter2, 
                orthogonal=FALSE, prior,
                fast, beta.blocks,
                X=X.wrk, Y=Y.wrk, XX=XX.wrk, sum.xy=sum.xy.wrk,
                seed=seed, verbose=(verbose),
                bigp.smalln=bigp.smalln, turn.sigma.on=turn.sigma.on, r.eff=r.eff)
    b.m <- gibbs$b.m
    b.m[is.na(b.m)] <- 0  #remove NA's from the bma (rare)    
    complexity.vec <- gibbs$complexity.vec
    hyperv <- gibbs$hyperv
    model <- gibbs$model
    if (turn.sigma.on) mse.hat <- mean(gibbs$sigma.vec, na.rm=T)
    phat.bma <- min(sum(abs(b.m) > eps, na.rm = TRUE), max.var, na.rm=T)
    resid <- (sum.xy.wrk - XX.wrk %*% b.m)
    penal <- resid / b.m
     
  }

  ### ---------------------------------
  # Save terms
  # Rescale beta: careful with X.org.sd=0
  # !!!NOTE!!!! bma.scale is used for prediction

  b.m <- b.m/sf
  b.m[is.na(b.m)] <- 0
  bma <- b.m
  bma.scale <- bma/X.org.sd
  bma.scale[X.org.sd==0] <- NA
  if (is.na(phat.bma)) phat.bma <- 0


  ### ---------------------------------
  # Order variables using the bma
  # Track the ordered variables used to fit the gnet

  if (verbose) cat("\t primary loop completed... \r")
  o.r <- order(abs(bma), decreasing = TRUE)
  if (phat.bma > 0) {
    gnet.obj.vars <- o.r[1:phat.bma]
  }
  else {
    gnet.obj.vars <- NULL
  }

### --------------------------------------------------------------
###
###  Variable Selection via generalized elastic net (gnet)
###
### gnet solution corresponds to ellipsoid optimization around GRR estimator
### closest to the bma; model selection based on AIC 
###
### --------------------------------------------------------------

if (verbose) cat("\t generalized elastic net (gnet) variable selection...                \r")
if (length(nozap.pt) > 0) {
  penal.constant <- mean(abs(resid[nozap.pt]), na.rm = TRUE)
  penal <-  (sqrt(n.data) * abs(penal))  / penal.constant
  #adjust penalty constant if random effects present)
  if (is.null(r.eff)) {
    penal.new <- penal
    for (k in 1:length(r.eff)) {
      match.k <- match(r.effects[[k]], beta.names)
      match.k <- match.k[!is.na(match.k)]
      penal.constant.k <- mean(abs(resid[match.k]), na.rm = TRUE)
      penal.new[match.k] <-  (sqrt(n.data) * abs(penal[match.k]))  / penal.constant.k
    }
    penal <- penal.new
  }
  gnet.out <- gnet.get(Y.center, X.wrk, o.r, phat.bma, penal, mse = mse.hat)
}
else {
  gnet.out <- list(gnet = rep(0, n.cov), gnet.path = NULL, gnet.obj = NULL)
}
###bad case when lasso fails
if (is.null(gnet.out$gnet)) {
  gnet.out$gnet <- rep(0, n.cov)
  gnet.out$gnet.path <- NULL
}
###extract objects
gnet <- gnet.out$gnet
gnet.path <- gnet.out$gnet.path
gnet.obj <- gnet.out$gnet.obj
gnet.scale <- gnet/X.org.sd
gnet.scale[X.org.sd==0] <- NA
phat <- sum(abs(gnet) > eps, na.rm = TRUE)
###scale and center gnet object/path information
###used later for prediction and lars-processing
if (!is.null(gnet.path)) {
  gnet.path$path <- t(t(gnet.path$path)/X.org.sd)
}
else {
  gnet.path$aic <- NULL
  gnet.path$path <- matrix(0, 1, n.cov)
}
if (!is.null(gnet.obj)) {
  gnet.obj$mu <- Y.org.mean
  gnet.obj$meanx <- X.org.mean
}
  

  
### --------------------------------------------------------------
###
###  SUMMARY VALUES (MAIN FUNCTION COMPLETED)
###
### --------------------------------------------------------------

ss.summary <- as.data.frame(cbind(bma=bma, gnet=gnet, bma.scale=bma.scale, gnet.scale = gnet.scale))
rownames(ss.summary) <- beta.names
ss.summary <- ss.summary[o.r, ]
names(bma) <- names(bma.scale) <- names(gnet) <- names(gnet.scale) <- beta.names
if (!is.null(r.eff)) {#nice labels for complexity matrix (for random effects)
  complexity.vec <- as.matrix(complexity.vec)
  colnames(complexity.vec) <- paste("f.eff.", ncol(complexity.vec):1, sep = "")
  colnames(complexity.vec)[1:length(r.eff.names)] <- paste("r.eff.", r.eff.names, sep = "")
}

### --------------------------------------------------------------
###	Terminal Output
###     Save as list
### --------------------------------------------------------------	

verbose.list <- list(
  c(method),
  c(bigp.smalln),
  c(screen),
  c(fast),
  c(n.data),
  c(n.cov),
  c(n.iter1),
  c(n.iter2),
  c(round(mean(mse.hat), 4)),
  c(phat)
)
  
if (verbose){
  
    cat("-------------------------------------------------------------------","\n")
    cat("Variable selection method     :",verbose.list[[1]],"\n")
    cat("Big p small n                 :",verbose.list[[2]],"\n")
    cat("Screen variables              :",verbose.list[[3]],"\n")
    cat("Fast processing               :",verbose.list[[4]],"\n")
    cat("Sample size                   :",verbose.list[[5]],"\n")
    cat("No. predictors                :",verbose.list[[6]],"\n")
    cat("No. burn-in values            :",verbose.list[[7]],"\n")
    cat("No. sampled values            :",verbose.list[[8]],"\n")
    cat("Estimated mse                 :",verbose.list[[9]],"\n")
    cat("Model size                    :",verbose.list[[10]],"\n")
    cat("\n\n")
    cat("---> Top variables:\n")
    print(round(ss.summary[ss.summary[, 2] != 0, ], 3))
    cat("-------------------------------------------------------------------","\n")
}

### --------------------------------------------------------------
###
###    Return the goodies
###	
### --------------------------------------------------------------	


out <- list(
       summary=ss.summary,                     #ordered summary object
       verbose=verbose.list,                   #verbose details (for printing)
       terms=mt,                               #terms
       sigma.hat=mse.hat,                      #internal (or external) mse
       y=Y.org,                                #original Y
       xnew=X.wrk,                             #centered rescaled X
       x=X.org,                                #original X            
       y.center=Y.org.mean,                    #centering for Y (0 if center=F)
       x.center=X.org.mean,                    #centering for original X
       x.scale=X.org.sd,                       #scaling for original X
       names=beta.names,                       #variable names 
       bma=bma,                                #bma coefficients for xnew
       bma.scale=bma.scale,                    #bma coefficients rescaled for x
       gnet=gnet,                              #gnet coefficients for xnew
       gnet.scale=gnet.scale,                  #gnet coefficients rescaled for x 
       gnet.path=gnet.path,                    #gnet full solution path
       gnet.obj=gnet.obj,                      #gnet object (lars type)
       gnet.obj.vars=gnet.obj.vars,            #variables (in order) used to define gnet
       gnet.parms=penal,                       #grr parameters used to define gnet
       phat=phat,                              #estimated dimension
       complexity=complexity.vec,              #complexity estimates
       ridge=hyperv,                           #gamma hypervariance
       model=model                             #list of models sampled (can be NULL)
       )

class(out) <- "spikeslab"
return(out)

}
