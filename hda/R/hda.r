#############################
### Definition of method hda
hda <- function (x, ...) 
    UseMethod("hda")

################################################
### formuala interface for hda
hda.formula <- function(formula, data = NULL, ...)
{
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval.parent(m$data))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if ((yvar <- attr(Terms, "response")) > 0) 
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    if (xint > 0) 
        x <- x[, -xint, drop = FALSE]
    res <- hda(x, grouping, ...)
    res$terms <- Terms
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- xlev
    attr(res, "na.message") <- attr(m, "na.message")
    if (!is.null(attr(m, "na.action"))) 
        res$na.action <- attr(m, "na.action")
    res
}

##############################################
### print function
print.hda <- function(x, ...){
  cat("\nHeteroscedastic discriminant analysis\n\n")
  cat("Call:\n")
  print(x$hda.call)

  cat("\nDimension of reduced space: \n", x$reduced.dimension, "\n\n")
  
  cat("\nClass means in reduced space: \n")
  print(x$class.dist$new.classmeans)
  cat("\n\n")

  cat("\nClass covariances in reduced space: \n")
  print(x$class.dist$new.classcovs)
  cat("\n\n")

  cat("\nComponent class discrimination accuracies in reduced space: \n")
  print(x$comp.acc)
  cat("\n\n")
  
  
  invisible(x)
  }

######################
### visualize loadings
showloadings <- function(object, comps = 1:object$reduced.dimension, loadings = TRUE, ...){
  
  if (max(comps) > ncol(object$hda.loadings)) 
    stop("Component ids have to be <= dimension of object$hda.loadings!")
  
  vnames <- rownames(object$hda.loadings)
  d <- length(comps)
  
  # if no loadings should be plotted: replace loadings by lifts
  if(!loadings) {
    object$hda.loadings <- t(object$vlift$total.lift)
    if(d == 1) object$hda.loadings <- t(object$hda.loadings)
  }
  if (max(comps) > ncol(object$hda.loadings)) 
    stop("Component ids for loadings == FALSE have to be <= reduced dimension!")
  
  if(d > 2){
    op <- par(mfrow=c(d,d))
    for(i in comps){
      for(j in comps){
        plot(object$hda.loadings[,c(j,i)], type="n", ...)
        for(k in 1:nrow(object$hda.loadings)) 
          text(object$hda.loadings[k,j], object$hda.loadings[k,i], vnames[k])        
      }
    }
  }
  
  if(d==2){
    op <- par(mfrow=c(1,1))
    plot(object$hda.loadings[,c(comps[1],comps[2])], type="n", ...)
    for(k in 1:nrow(object$hda.loadings)) 
      text(object$hda.loadings[k,comps[1]], object$hda.loadings[k,comps[2]], vnames[k])     
  }
  
  if(d==1){
    op <- par(mfrow=c(1,1))
    plot(object$hda.loadings[,comps], type="n", xlab="Variable index", ylab=colnames(object$hda.loadings)[comps], ...)
    for(k in 1:nrow(object$hda.loadings)) 
      text(k, object$hda.loadings[k,comps], vnames[k])    
  }
  
  par(op)
}


######################
### plot function to visualize scores
plot.hda <- function(x, comps = 1:x$reduced.dimension, scores = TRUE, col = x$grouping, ...){
  if (max(comps) > nrow(x$hda.loadings)) 
    stop("Component ids have to be <= dimension of object$hda.loadings!")
  # loadings plot
  if(scores){
    if (length(comps) > 1) 
      plot(as.data.frame(x$hda.scores[,comps]), col = col, main = "Scores", ...)
    if (length(comps) == 1) 
      plot(x$hda.scores[,comps], col = col, ylab = paste("comp",comps,sep=" "), main = "Scores", ...)
  }
  # density plot
  if(!scores){
    priors <- table(x$grouping) / length(x$grouping)
    ncl <- length(table(x$grouping))
    d1 <- d2 <- ceiling(sqrt(length(comps)))
    for(i in 3:1){
      for(j in 3:1) if ( (i*j) >= length(comps)) {d1 <- i; d2 <- j}
    }
    par(mfrow=c(d1,d2))
    for(d in comps){
      k <- 1
      mean.kd <- x$class.dist$new.classmeans[[k]][d]
      if(is.matrix(x$class.dist$new.classcovs[[k]])) sd.kd <- sqrt(x$class.dist$new.classcovs[[k]][d,d])
      if(!is.matrix(x$class.dist$new.classcovs[[k]])) sd.kd <- sqrt(x$class.dist$new.classcovs[[k]]) # distinction: for newdim=1 classcovs is no matrix
      xpts <- seq(min(x$hda.scores[,d]), max(x$hda.scores[,d]), length.out = 500)
      ymax <- max(priors[k]*dnorm(xpts, mean = mean.kd, sd = sd.kd))
      if(is.matrix(x$class.dist$new.classcovs[[k]])) for (j in 2:ncl) ymax <- max(ymax, priors[j]*dnorm(xpts, mean = x$class.dist$new.classmeans[[j]][d], sd = sqrt(x$class.dist$new.classcovs[[j]][d,d])))
      if(!is.matrix(x$class.dist$new.classcovs[[k]])) for (j in 2:ncl) ymax <- max(ymax, priors[j]*dnorm(xpts, mean = x$class.dist$new.classmeans[[j]][d], sd = sqrt(x$class.dist$new.classcovs[[j]])))
      
      
      plot(xpts, priors[k]*dnorm(xpts, mean = mean.kd, sd = sd.kd), col = k, type = "l",
           ylim = c(0, ymax), xlab = colnames(x$hda.scores)[d], ylab = "prior[k] * f(x|k)")
      lines(rep(mean.kd, 2), c(0, priors[k]*dnorm(mean.kd, mean.kd, sd.kd)), lty = "dotted", col =k)
      for(k in 2:ncl){
        mean.kd <- x$class.dist$new.classmeans[[k]][d]
        if(is.matrix(x$class.dist$new.classcovs[[k]])) sd.kd <- sqrt(x$class.dist$new.classcovs[[k]][d,d])
        if(!is.matrix(x$class.dist$new.classcovs[[k]])) sd.kd <- sqrt(x$class.dist$new.classcovs[[k]])
        xpts <- seq(min(x$hda.scores[,d]), max(x$hda.scores[,d]), length.out = 500)
        lines(xpts, priors[k]*dnorm(xpts, mean = mean.kd, sd = sd.kd), col = k)
        lines(rep(mean.kd, 2), c(0, priors[k]*dnorm(mean.kd, mean.kd, sd.kd)), lty = "dotted", col =k)        
      }
    }        
  }
}

###########################################################################
### predict function for easy transformation of new data with a given model
predict.hda <- function(object, newdata, alldims = FALSE, task = c("dr", "c"), ...){
  if (class(object) != "hda") 
    stop("Object must be of class hda!")
  task <- match.arg(task)
  if(is.data.frame(newdata)) 
    newdata <- as.matrix(newdata)
  if(!is.matrix(newdata)) 
    stop("Newdata must be of type matrix or data frame!")
  if(dim(object$hda.loadings)[2] != dim(newdata)[2]) 
    stop("Newdata must be of same dimension as the explanatory input variables data set!")
  if(task == "dr"){
    new.transformed.data <- newdata %*% object$hda.loadings
    if(!alldims) 
        new.transformed.data <- new.transformed.data[,1:object$reduced.dimension,drop=FALSE]
    return(new.transformed.data)
  }
  if (task == "c"){
    if(class(object$naivebayes) != "naiveBayes") 
        stop("Classification of newdata can only be done id option 'crule = TRUE' has been chosen at hda() call")
    new.transformed.data <- newdata %*% object$hda.loadings
    new.transformed.data <- new.transformed.data[,1:object$reduced.dimension,drop=FALSE]
    prediction <- "Remove argument type = 'raw' in predict() call of naiveBayes() call or set it to 'class'."
    posteriors <- "Remove argument type = 'class' in predict() call of naiveBayes() call or set it to 'raw'."
    if (object$reduced.dimension > 1){
      if(is.null(object$naivebayes$levels))
        warning("Function naiveBayes() {e1071} requires vector of type factor for prediction of type 'class'!")
      if (sum(ls() == "type") == 0 ){
        prediction <- predict(object$naivebayes, new.transformed.data, type = "class", ...)
        posteriors <- predict(object$naivebayes, new.transformed.data, type = "raw", ...)
        }
      if (sum(ls() == "type") > 0){
        nbpred      <- predict(object$naivebayes, new.transformed.data, ...)
        if (length(dim(nbpred))  > 0) posteriors  <- nbpred
        if (length(dim(nbpred)) == 0) prediction  <- nbpred
        remove(nbpred)
        }
      }
    if (object$reduced.dimension == 1){
      warning("Function predict.naiveBayes {e1071} not implemented for dimension 1 of the reduced space! Prediction is done without using naiveBayes().\n")
    
      priors <- object$naivebayes$apriori/ sum(object$naivebayes$apriori)  
      classids <- 1:nrow(object$naivebayes$tables[[1]])  ### muss hier anstatt x der variablenname hin?
    
      priodens <- function(i) return(priors[i]*dnorm(new.transformed.data, mean = object$naivebayes$tables[[1]][i,1], sd = object$naivebayes$tables[[1]][i,2]))
      posteriors <- sapply(classids, priodens)
      posteriors <- t(apply(posteriors, 1, function(x) return(x / sum(x))))
      prediction <- as.factor(rownames(object$naivebayes$tables[[1]])[apply(posteriors, 1, which.max)])
    }
    
    result <- list(prediction, posteriors, new.transformed.data)
    names(result) <- c("prediction", "posteriors", "new.transformed.data")
    return(result)
    }                     
  }

##############################################################################
### calculation of heteroscedastic discriminant analysis loadings matrix 
compute.loadings <- function(covarray, clsizes, newd, initial = NULL, iters = 7, initers = 10, verb = TRUE){
    
    # covarray: array (of dimension oldd, oldd, classnumber) of the estimated class-specific covariance matrices (of size (oldd, oldd)
    # clsizes : number of corresponding observations per class
    # newd    : dimension of the desired subspace
    # iters   : number of optimization iterations
    # initers : iterations inner loop
    
    if (dim(covarray)[1] != dim(covarray)[2]) stop("Covariance matrices must be quadratic")
    if (newd >= dim(covarray)[2]) stop("Dimension of reduced space should be lower than original")
 
    ncl             <- dim(covarray)[3]
    oldd            <- dim(covarray)[1]
    totN            <- sum(clsizes)  # Total number of observations
    commonwithin    <- matrix(rowSums(sapply(1:ncl, function(i) clsizes[i] * covarray[,,i])), ncol=oldd, nrow=oldd) / totN

    # initialization loadings matrix
    Trafo           <- initial  
    invG            <- array(dim=c(oldd,oldd,oldd))

    for (iter in 1:iters){
      Qu  <- totN * log(det(Trafo)^2) - totN * oldd * (log(2*pi)+1)
      for (i in 1:oldd){
          if (i <= newd){
              G <- matrix(0,oldd,oldd)
              for (m in 1:ncl){
                  sigma_i <- as.numeric(t(Trafo[,i]) %*% covarray[,,m] %*% Trafo[,i])
                  G       <- G + clsizes[m] / sigma_i * covarray[,,m] 
                  Qu       <- Qu - clsizes[m] * log(sigma_i)
                  }
              }
          if (i > newd){
              sigma_i <- as.numeric(t(Trafo[,i]) %*% commonwithin %*% Trafo[,i])
              G       <- totN / sigma_i * commonwithin
              Qu       <- Qu - totN * log(sigma_i)
              }
           
          invG[,,i] <- solve(G)
          }
          
        Qu <- Qu/2
        if (verb) cat(paste("Iteration", iter,"Log Likelihood: ",  Qu, "\n"))
        for (initer in 1:initers){
            for (i in 1:oldd){
                Ce          <- t(solve(t(Trafo)) *  det(t(Trafo)))
                ci_invG     <- Ce[i,] %*% invG[,,i]
                Trafo[,i]  <- t(ci_invG * sqrt(totN / as.numeric((ci_invG %*% t(t(Ce[i,]))))))
                }
            }
        }
     
    Trafo = Trafo / (det(Trafo)^(1/oldd))
    return(Trafo)
    } 

#####################################################################################
### regularize covariance matrices according to input parameters as in Friedman (1989)  
regularize <- function(covarray, clsizes, lambd, gamm){

  if(lambd < 0 || lambd > 1) 
    stop("lambd and gamm must be in [0,1]")
  if(gamm < 0 || gamm > 1) 
    stop("lambd and gamm must be in [0,1]")

  # calculation of common covaiance matrix (over al classes)
  commoncov <- (clsizes[1]-1) * covarray[,,1]
  for(i in 2:dim(covarray)[3]) 
    commoncov <- commoncov + (clsizes[i]-1) * covarray[,,i]
  commoncov <- commoncov / (sum(clsizes)-dim(covarray)[3])
  
  # regularization of covaiances - array
  for(i in 1:dim(covarray)[3]){
    # towards equal covariances
    covarray[,,i] <- (lambd * clsizes[i] * covarray[,,i] + (1-lambd) * sum(clsizes) * commoncov) / 
                     (lambd * clsizes[i] + (1-lambd) * sum(clsizes)) 
    # towards diagonality with average variance
    average.variance.i  <- sum(diag(covarray[,,i]))/dim(covarray)[1]
    shrinkcov.i         <- diag(average.variance.i, dim(covarray)[1])
    covarray[,,i]       <- gamm * (covarray[,,i]) + (1-gamm) * shrinkcov.i     
    }  

  return(covarray)
  }

######################################################  
### main call of heteroscedastic discriminant analysis 
# ... calculation of the input parameters from data set and call of compute.loadings        
hda.default <- function(x, grouping, newdim = 1:(ncol(x)-1), crule = FALSE, 
    reg.lamb = NULL, reg.gamm = NULL, initial.loadings = NULL, 
    sig.levs = c(0.05,0.05), noutit = 7, ninit = 10, verbose = TRUE, ...){
  # x:                  data frame containing the (metric) variables  
  # grouping:           vector of class labels
  # newdim:             dimension of reduced space
  # crule:              TRUE if naiveBayes() model shoulf be built using pkg e1071
  # reg.lamb,reg.gamm:  regularization parameters as in rda() [Friedman, 1989]
  # initial.laodings:   (optional) initial guess of the (quadratic) loadings matrix 
  # sig.levs:           significance levels for tests if no unique newdim is specified
  # noutit:             (optional) number of outer iterations of the algorithm
  # ninit:              (optional) number of inner iterations of the algorithm

  # check for possible errors in function call
  if(length(table(grouping)) < 2) 
    stop("Class vector should contain different levels")
  if(dim(x)[2] < 2) 
    stop("Dimensionality reduction only meaningful if x has at least two coloumns")
  
  # class parameters
  clsizes <- table(grouping)
  clnb    <- length(clsizes)
  # array of covariance matrices:
  covlist <- by(as.data.frame(x), grouping, cov)
  # dimension of the original data space
  dms     <- dim(covlist[[1]])[1]
  carray  <- array(0,c(dms,dms,clnb))
  for (i in seq_along(covlist)){
    carray[,,i] <- covlist[[i]]
  }
  # regularization ov covariance estimates if specified by input parameters  
  if (sum(c(is.null(reg.lamb),is.null(reg.gamm))) < 2) 
    carray <- regularize(carray, clsizes, reg.lamb, reg.gamm)

  
  # initialization loadings matrix
  if (is.character(initial.loadings)) {if (initial.loadings == "random"){
    cat("Initialization by random orthonormal matrix.\n")
    initial.loadings   <- qr.Q(qr(matrix(runif(dms^2), dms, dms)))
  }}
  else if(is.matrix(initial.loadings)){
    if (sum(dim(initial.loadings)== (dim(carray)[1:2])) != 2) 
      stop("(Optional) initalization of loading matrix must be quadratic and of size dim(x)[2]")
  }    
  else if (is.null(initial.loadings)){
    cat("Initialization by the identity.\n")
    initial.loadings   <- diag(dms)
  }
  else stop("Incorrect specification of initial.loadings!")
  
  
  hda.loadings <- NULL
  trace.newdim <- NULL
  # recursive call of hda.default if not a single newdim as specified
  if(length(newdim) > 1){
    pvaleqmeans     <- 0
    pvalhomogcovs   <- 0
    d               <- 1 
    while(((d <= length(newdim)) * ((pvaleqmeans <= sig.levs[1]) + (pvalhomogcovs <= sig.levs[2]))) > 0){
      if (verbose) 
        cat("\nnewdim = ", newdim[d], "\n")
      dummy <- hda(x=x, grouping=grouping, newdim=newdim[d], reg.lamb = reg.lamb, 
                   reg.gamm = reg.gamm, initial.loadings = initial.loadings, 
                   sig.levs = sig.levs, noutit = noutit, ninit = ninit, verbose = verbose)
      if((ncol(x)-newdim[d]) > 1) 
        pvaleqmeans <- dummy$eqmean.test[[3]][6]
      if((ncol(x)-newdim[d]) == 1) 
        pvaleqmeans <- dummy$eqmean.test[[3]][1,5]
      pvalhomogcovs <- dummy$homog.test$pValue
      trace.newdim  <- cbind(trace.newdim, c(pvaleqmeans, pvalhomogcovs)) 
      rownames(trace.newdim) <- c("eqmean.tests","homog.tests")
      colnames(trace.newdim) <- newdim[1:d]
      d <- d+1
    }
    if (verbose) cat("\n")
    hda.loadings  <- dummy$hda.loadings
    newdim    <- newdim[d-1]
  }           


  # calculate transformation matrix if not already done (i.e. if 
  if (is.null(hda.loadings)){
    if(round(newdim) != newdim) 
        stop("newdim must be an integer")
    hda.loadings <- compute.loadings(carray, clsizes, newdim, 
                                     initial = initial.loadings, iters = noutit, 
                                     initers = ninit, verb = verbose)
  }
  rownames(hda.loadings) <- colnames(as.data.frame(x))
  compnames  <- NULL
  for(i in 1:dms) compnames <- c(compnames, paste("comp",i,sep=""))
  colnames(hda.loadings) <- compnames 
  
  # calculate the transformed space
  newspace  <- as.data.frame(as.matrix(x) %*% hda.loadings)
  
  # class distribution parameters in new space
  new.classmeans  <- by(as.data.frame(newspace[,1:newdim]),grouping,colMeans)
  new.classcovs   <- by(as.data.frame(newspace[,1:newdim]),grouping,cov)
  new.classdist   <- list(new.classmeans = new.classmeans, new.classcovs = new.classcovs)
  
  ### several tests of assumptions  
  # test on homogenity of the classwise covariance matrices in remaining dimensions 
  newcovs     <- by(newspace, grouping, cov) # ml-estimates of classes in redundant dimensions
  newcommon   <- matrix(0,dms-newdim,dms-newdim) # initialize common covariance matrix in redundant transformed space
  stat        <- 0 # initialize value of statistic
  for(i in seq_along(clsizes)){ 
    newcovs[[i]]  <- (newcovs[[i]])[(newdim+1):dms,(newdim+1):dms] * (clsizes[i]-1)/clsizes[i]
    stat          <- stat - clsizes[i] * log(det(as.matrix(newcovs[[i]])))
    newcommon     <- newcommon + newcovs[[i]] * clsizes[i] 
  }
  newcommon       <- newcommon / sum(clsizes)
  stat           <- as.numeric(stat + sum(clsizes) * log(det(as.matrix(newcommon))))
  dfs            <- (dms-newdim)*(dms-newdim+1)*(length(clsizes)-1)/2
  pval           <- 1 - pchisq(stat, dfs) # corresponding p-value
  homog.test     <- list(new.common.covariance = newcommon, 
                         new.class.covariances = newcovs, statistic = stat, 
                         dfs = dfs, pValue = pval) 
  
  # test on equal class means in remaining dimension
  new.classmeans2  <- by(as.data.frame(newspace[,(newdim+1):dms]),grouping,colMeans)
  new.classcovs2   <- by(as.data.frame(newspace[,(newdim+1):dms]),grouping,cov)
  if((dms-newdim) > 1) 
    stats <- summary(manova(as.matrix(newspace[,(newdim+1):dms])~grouping), test="Wilks")[[4]][1,]
  if((dms-newdim) == 1) 
    stats <- summary(aov(as.matrix(newspace[,(newdim+1):dms])~grouping))[[1]]
  eqmean.test <- list(new.class.means=new.classmeans2, new.class.covariances=new.classcovs2, stats)
  
  cl <- match.call()
  cl[[1]] <- as.name("hda")
  
  if (crule){
    crule <- naiveBayes(newspace[,1:newdim], grouping)
  }
  
  # compute component accuracies
  comp.accs <- comp.acc(list(grouping = grouping, class.dist = new.classdist))
  if (ncol(comp.accs) == 1) colnames(comp.accs) <- colnames(hda.loadings)[1]
  
  priors <- table(grouping) / length(grouping)  
  objkt <- list(x = x, priors = priors, grouping = grouping, reduced.dimension = newdim, comp.acc = comp.accs, hda.loadings = hda.loadings)
  vlift <- comp.vlifts(objkt)
  rm(objkt)
  
  # computation of componentwise accuracy (w.r.t training priors)
  # [not done earlier as com.acc is used by function comp.vlifts() without total accuracy row]
  tocomp.acc <- priors %*% comp.accs
  rownames(tocomp.acc) <- "total"
  comp.accs <- rbind(comp.accs, tocomp.acc)
  
  result <- list(hda.loadings = hda.loadings, hda.scores = newspace, 
                 grouping = grouping, class.dist = new.classdist, 
                 priors = priors, reduced.dimension = newdim, naivebayes = crule, 
                 comp.acc = comp.accs, vlift = vlift,   
                 reg.lamb = reg.lamb, reg.gamm = reg.gamm, eqmean.test = eqmean.test, 
                 homog.test = homog.test, hda.call = cl, initial.loadings = initial.loadings, 
                 trace.dimensions = trace.newdim)
                 
  class(result) <- "hda"
  return(result) 
}

comp.acc <- function(object){
  #  number of classes and hda discriminant components
  #bject$"class.dist"[[2]]
  ncl     <- length(object$"class.dist"[[1]])
  ncomps  <- length(object$"class.dist"[[1]][[1]])
  priors  <- table(object$"grouping") / length(object$"grouping")
    
  if(ncomps > 1){
    # means in hda space
    cl.means <- object$"class.dist"[[1]][[1]]
    for(i in 2:ncl) cl.means <- rbind(cl.means, object$"class.dist"[[1]][[i]])
    # sdevs in hda space
    cl.sdevs <- sqrt(diag(object$"class.dist"[[2]][[1]]))
    for(i in 2:ncl) cl.sdevs <- rbind(cl.sdevs, sqrt(diag(object$"class.dist"[[2]][[i]])))
  }
  if(ncomps == 1){
    # means in hda space
    cl.means <- matrix(as.vector(object$"class.dist"[[1]]), nrow = ncl, ncol = 1)
    # sdevs in hda space
    cl.sdevs <- matrix(as.vector(object$"class.dist"[[2]]), nrow = ncl, ncol = 1)
  }
  rownames(cl.means) <- names(object$"class.dist"[[1]])
  rownames(cl.sdevs) <- names(object$"class.dist"[[1]])
  
  # compute accurcies of any class in any component
  comp.acc <- matrix(NA, nrow = ncl, ncol = ncomps) 
  rownames(comp.acc) <- names(object$"class.dist"[[1]])
  colnames(comp.acc) <- names(object$"class.dist"[[1]][[1]])
  
  for(i in 1:ncl){ #for all classes
    for(j in 1:ncomps){ #for all components
      # define vector of other classes
      ocls <- ((1:ncl)[-i]) 
      # set vector of decision boundaries for class i (has to be updated for each comparison with any class ocl)
      dbounds <- NULL
      for(ocl in ocls){ 
        m1  <- cl.means[i,j] 
        m2  <- cl.means[ocl,j]
        sd1 <- cl.sdevs[i,j]
        sd2 <- cl.sdevs[ocl,j]
        pr1 <- priors[i]
        pr2 <- priors[ocl]                      
        
        # calculate decision boundaries for equal variances 
        # note: intervals for other class are computed
        if (sd1 == sd2){
          x1 <- (m1 + m2)/2
          if(m1 < m2) dbounds <- rbind(dbounds, c(x1, Inf)) 
          if(m1 > m2) dbounds <- rbind(dbounds, c(-Inf, x1))
        }
        
        # calculate decision boundaries for unequal variances
        if (sd1 != sd2){
          pe <- 2*(m1*sd2^2 - m2*sd1^2)/ (sd1^2-sd2^2)
          # includes prior correction                      
          qu <- (sd2^2*m1^2 - sd1^2*m2^2 - sd1^2 * sd2^2 * 2*(log(sd2)-log(sd1)+log(pr1)-log(pr2))) / (sd2^2-sd1^2)          
          if((pe^2/4 - qu) >= 0){# usual case  
            x1 <- -pe/2 - sqrt(pe^2/4 - qu)
            x2 <- -pe/2 + sqrt(pe^2/4 - qu)
            if (sd1 > sd2) dbounds <- rbind(dbounds, c(x1, x2)) # higher density within boundaries 
            if (sd1 < sd2) dbounds <- rbind(dbounds, c(-Inf, x1), c(x2, Inf)) # higher density outside boundaries
          }          
          if((pe^2/4 - qu) < 0){# exception: for unequal priors an similar m1,m2 and sd1, sd2 one density can dominate the other   
            if (sd1 > sd2) dbounds <- rbind(dbounds, c(Inf, Inf)) # f1(x) > f2(x) for all x: no area with higher likelihood of other class 
            if (sd1 < sd2) dbounds <- rbind(dbounds, c(-Inf, Inf)) # f1(x) < f2(x) for all x: no decision for class one in this component 
          }
        }
      }
      # merge and align univariate classification bounds
      if(nrow(dbounds) > 1){
        sortbounds <- sort(dbounds[,1], index.return = T)
        dbounds <- dbounds[sortbounds$ix,]
      }
      # build union of intervals
      if(nrow(dbounds) > 1){
        k <- 2
        while(k <= nrow(dbounds)){
          if(dbounds[k, 1] <= dbounds[k-1, 2]){# integrate interval into previous one and set to NA
            if(dbounds[k, 2] > dbounds[k-1, 2]){# update upper bound if necessary
              dbounds[k-1, 2] <- dbounds[k, 2]
            }
            if(nrow(dbounds) == 2) dbounds <- matrix(dbounds[-k,], nrow=1) # convert resulting vector into matrix
            if(nrow(dbounds) > 2) dbounds <- dbounds[-k,]            
            k <- k-1
          }
          k <- k+1
        }
      }
      # compute probabilities between class i - intervalls 
      probs <- pnorm(dbounds, m1, sd1)
      probs <- apply(probs, 1, diff)
      err.cli <- sum(probs)
      comp.acc[i,j] <- 1-err.cli
    }
  }
  return(comp.acc)
}



comp.vlifts <- function(obj){
  if (obj$reduced.dimension > 1){
    vlift <- array(NA, c(dim(obj$comp.acc), ncol(obj$x)))
    dimnames(vlift) <- list(levels(obj$grouping), colnames(obj$comp.acc), rownames(obj$hda.loadings))
    if (dim(vlift)[2] == 1) dimnames(vlift)[[2]] <- colnames(obj$hda.loadings)[1]
    for(i in 1:ncol(obj$x)){
      # create loading vector without component i
      hdaloads.i <- obj$hda.loadings
      hdaloads.i[i,] <- 0
      # ...and transformed space
      newspace.woi <- as.data.frame(as.matrix(obj$x) %*% hdaloads.i[,1:ncol(obj$comp.acc)])
      
      # compute means and covariances
      classmeans  <- by(newspace.woi, obj$grouping, colMeans)
      classcovs   <- by(newspace.woi, obj$grouping, cov)
      vlift[,,i] <- comp.acc(list(grouping = obj$grouping, class.dist = list(classmeans, classcovs)))
    }
    # currently vlift contains accuracies of each neutralized variable for any class in each component (array)
    
    # compute overall lift of each component in each component (accuracies weighted by priors)
    vlift.tot <- apply(vlift,2:3,function(z) return(sum(z * as.numeric(obj$priors))))    
    #does not work: for(i in 1:ncol(obj$x)) vlift.tot[,i] <-  (obj$priors %*% obj$comp.acc) / vlift.tot[,i]
    comp.tot.accs <- (apply(obj$comp.acc, 2, function(z) return(sum(z * as.numeric(obj$priors)))))
    for(i in 1:ncol(obj$x)) vlift.tot[,i] <-  comp.tot.accs / vlift.tot[,i]
    
    rownames(vlift.tot) <- colnames(obj$comp.acc)
    colnames(vlift.tot) <- rownames(obj$hda.loadings)
    
    # compute lift of each variable for any class in each component (array)  
    for(i in 1:ncol(obj$x)) vlift[,,i] <-  obj$comp.acc / vlift[,,i]
    
    # create result object 
    vlift <- list("total.lift" = vlift.tot, "classwise.lift" = vlift) 
  }

  if (obj$reduced.dimension == 1){  
    vlift <- matrix(NA, nrow = length(obj$comp.acc), ncol = ncol(obj$x))
    dimnames(vlift) <- list(levels(obj$grouping), rownames(obj$hda.loadings))
    if (dim(vlift)[2] == 1) dimnames(vlift)[[2]] <- colnames(obj$hda.loadings)[1]
    for(i in 1:ncol(obj$x)){
      # create loading vector without component i
      hdaloads.i <- obj$hda.loadings
      hdaloads.i[i,] <- 0
      # ...and transformed space
      newspace.woi <- as.data.frame(as.matrix(obj$x) %*% hdaloads.i[,1])
      
      # compute means and covariances
      classmeans  <- by(newspace.woi, obj$grouping, colMeans)
      classcovs   <- by(newspace.woi, obj$grouping, cov)
      vlift[,i] <- comp.acc(list(grouping = obj$grouping, class.dist = list(classmeans, classcovs)))
    }
    # currently vlift contains accuracies of each neutralized variable for any class in each component
    
    # compute overall lift of each component in each component (accuracies weighted by priors)
    vlift.tot <- apply(vlift,2,function(z) return(sum(z * as.numeric(obj$priors))))    
    #does not work: for(i in 1:ncol(obj$x)) vlift.tot[,i] <-  (obj$priors %*% obj$comp.acc) / vlift.tot[,i]
    comp.tot.accs <- sum(obj$comp.acc * as.numeric(obj$priors))
    vlift.tot <-  comp.tot.accs / vlift.tot    
    names(vlift.tot) <- rownames(obj$hda.loadings)
    
    # compute lift of each variable for any class in each component 
    for(i in 1:ncol(obj$x)) vlift[,i] <-  obj$comp.acc / vlift[,i]
    
    # create result object 
    vlift <- list("total.lift" = vlift.tot, "classwise.lift" = vlift) 
    
  }
  
  return(vlift)
}





#showloadings <- function(object, comps = 1:object$reduced.dimension, loadings = TRUE, ...){
#  if (max(comps) > nrow(object$hda.loadings)) 
#    stop("Component ids have to be <= dimension of object$hda.loadings")
#  vnames <- rownames(object$hda.loadings)
#  
#  #  
#  if(!loadings){
#    lftfrme <- t(object$vlift$classwise.lift[1,,]) 
#    clss <- rep(levels(object$grouping), each = nrow(lftfrme))
#    vnames <- rep(rownames(lftfrme), length(levels(object$grouping)))
#    for(i in 2:length(table(object$grouping))) lftfrme <- rbind(lftfrme, t(object$vlift$classwise.lift[i,,])) 
#    #    # invert lifts
#    #    lftfrme <- 1/lftfrme
#    if (max(comps) > nrow(lftfrme)) 
#      stop("Component ids have to be <= dimension of object$hda.loadings")
#    
#  }
#  
#  d <- length(comps)
#  if(d > 2){
#    op <- par(mfrow=c(d,d))
#    for(i in comps){
#      for(j in comps){
#        if(loadings){
#          plot(object$hda.loadings[,c(j,i)], type="n", ...)
#          for(k in 1:nrow(object$hda.loadings)) 
#            text(object$hda.loadings[k,j], object$hda.loadings[k,i], vnames[k])
#        }
#        if(!loadings){
#          plot(lftfrme[,c(j,i)], type="n", ...)
#          for(k in 1:nrow(lftfrme)) 
#            text(lftfrme[k,j], lftfrme[k,i], vnames[k], col  = as.integer(clss)[k])
#        }
#        
#      }
#    }
#  }
#  if(d==2){
#    op <- par(mfrow=c(1,1))
#    if(loadings){ 
#      plot(object$hda.loadings[,c(comps[1],comps[2])], type="n", ...)
#      for(k in 1:nrow(object$hda.loadings)) 
#        text(object$hda.loadings[k,comps[1]], object$hda.loadings[k,comps[2]], vnames[k]) 
#    }
#    if(!loadings){
#      plot(lftfrme[,c(comps[1],comps[2])], type="n", ...)
#      for(k in 1:nrow(lftfrme)) 
#        text(lftfrme[k,comps[1]], lftfrme[k,comps[2]], vnames[k], col = as.integer(clss)[k])       
#    }
#  }    
#  if(d==1){
#    op <- par(mfrow=c(1,1))
#    if(loadings){       
#      plot(object$hda.loadings[,comps], type="n", xlab="Variable index", ylab=colnames(object$hda.loadings)[comps], ...)
#      for(k in 1:nrow(object$hda.loadings)) 
#        text(k, object$hda.loadings[k,comps], vnames[k])
#    }
#    if(!loadings){
#      plot(lftfrme[,comps], type="n", xlab="Variable index", ylab=colnames(lftfrme)[comps], ...)
#      for(k in 1:nrow(lftfrme)) 
#        text(k, lftfrme[k,comps], vnames[k], col = as.integer(clss)[k])
#      
#    } 
#    
#  }     
#  par(op)
#}
#
