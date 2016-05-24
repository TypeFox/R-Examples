###########################
## convenience functions ##
###########################

## obtain the number/ID for all terminal nodes
terminal_nodeIDs <- function(node) {
  if(node$terminal) return(node$nodeID)
  ll <- terminal_nodeIDs(node$left)
  rr <- terminal_nodeIDs(node$right)
  return(c(ll, rr))
}


#########################
## workhorse functions ##
#########################

### determine which observations go left or right
mob_fit_childweights <- function(node, mf, weights) {

    partvar <- mf@get("part")
    xselect <- partvar[[node$psplit$variableID]]

    ## we need to coerce ordered factors to numeric
    ## this is what party C code does as well!

    if (class(node$psplit) == "orderedSplit") {
        leftweights <- (as.double(xselect) <= node$psplit$splitpoint) * weights
        rightweights <- (as.double(xselect) > node$psplit$splitpoint) * weights
    } else {
        leftweights <- (xselect %in%
            levels(xselect)[as.logical(node$psplit$splitpoint)]) * weights
        rightweights <- (!(xselect %in%
            levels(xselect)[as.logical(node$psplit$splitpoint)])) * weights
    }

    list(left = leftweights, right = rightweights)
}

### setup a new (inner or terminal) node of a tree
mob_fit_setupnode <- function(obj, mf, weights, control) {

    ### control parameters
    alpha <- control$alpha
    bonferroni <- control$bonferroni
    minsplit <- control$minsplit
    trim <- control$trim
    objfun <- control$objfun
    verbose <- control$verbose
    breakties <- control$breakties
    parm <- control$parm

    ### if too few observations: no split = return terminal node
    if (sum(weights) < 2 * minsplit) {
        node <- list(nodeID = NULL, weights = weights,
                     criterion = list(statistic = 0, criterion = 0, maxcriterion = 0),
                     terminal = TRUE, psplit = NULL, ssplits = NULL,
                     prediction = 0, left = NULL, right = NULL,
                     sumweights = as.double(sum(weights)))
        class(node) <- "TerminalModelNode"
        return(node)
    }

    ### variable selection via fluctuation tests
    test <- try(mob_fit_fluctests(obj, mf, minsplit = minsplit, trim = trim,
      breakties = breakties, parm = parm))

    if (!inherits(test, "try-error")) {
        if(bonferroni) {
          pval1 <- pmin(1, sum(!is.na(test$pval)) * test$pval)
          pval2 <- 1 - (1-test$pval)^sum(!is.na(test$pval))
          test$pval <- ifelse(!is.na(test$pval) & (test$pval > 0.01), pval2, pval1)
        }

        best <- test$best
        TERMINAL <- is.na(best) || test$pval[best] > alpha

        if (verbose) {
            cat("\n-------------------------------------------\nFluctuation tests of splitting variables:\n")
            print(rbind(statistic = test$stat, p.value = test$pval))
            cat("\nBest splitting variable: ")
            cat(names(test$stat)[best])
            cat("\nPerform split? ")
            cat(ifelse(TERMINAL, "no", "yes"))
            cat("\n-------------------------------------------\n")    
        }
    } else {
        TERMINAL <- TRUE
        test <- list(stat = NA, pval = NA)
    }

    ### splitting
    na_max <- function(x) {
      if(all(is.na(x))) NA else max(x, na.rm = TRUE)
    }
    if (TERMINAL) {
        node <- list(nodeID = NULL, weights = weights,
	             criterion = list(statistic = test$stat, 
                                      criterion = 1 - test$pval,
				      maxcriterion = na_max(1 - test$pval)),
                     terminal = TRUE, psplit = NULL, ssplits = NULL,
                     prediction = 0, left = NULL, right = NULL, 
                     sumweights = as.double(sum(weights)))
        class(node) <- "TerminalModelNode"
        return(node)
    } else {
        partvar <- mf@get("part")
        xselect <- partvar[[best]]
        thissplit <- mob_fit_splitnode(xselect, obj, mf, weights, minsplit = minsplit, 
                                       objfun = objfun, verbose = verbose)

        ## check if splitting was unsuccessful
        if (identical(FALSE, thissplit)) {
            node <- list(nodeID = NULL, weights = weights,
                         criterion = list(statistic = test$stat, 
                                          criterion = 1 - test$pval, 
                                          maxcriterion = na_max(1 - test$pval)),
                         terminal = TRUE, psplit = NULL, ssplits = NULL,
                         prediction = 0, left = NULL, right = NULL, 
                         sumweights = as.double(sum(weights)))
            class(node) <- "TerminalModelNode"  
            
            ### more confusion than information
	    ### warning("no admissable split found", call. = FALSE)
	    if(verbose)
	      cat(paste("\nNo admissable split found in ", sQuote(names(test$stat)[best]), "\n", sep = ""))	    
	    return(node)
        }

        thissplit$variableID <- best
        thissplit$variableName <- names(partvar)[best]
        node <- list(nodeID = NULL, weights = weights, 
                     criterion = list(statistic = test$stat, 
                                      criterion = 1 - test$pval, 
                                      maxcriterion = na_max(1 - test$pval)),
                     terminal = FALSE,
                     psplit = thissplit, ssplits = NULL, 
                     prediction = 0, left = NULL, right = NULL, 
                     sumweights = as.double(sum(weights)))
        class(node) <- "SplittingNode"
    }
    
    node$variableID <- best
    if (verbose) {
        cat("\nNode properties:\n")
        print(node$psplit, left = TRUE)
        cat(paste("; criterion = ", round(node$criterion$maxcriterion, 3), 
              ", statistic = ", round(max(node$criterion$statistic), 3), "\n",
              collapse = "", sep = ""))
    }
    node
}

### variable selection:
### conduct all M-fluctuation tests of fitted obj 
### with respect to each variable from a set of
### potential partitioning variables in mf
mob_fit_fluctests <- function(obj, mf, minsplit, trim, breakties, parm) {
  ## Cramer-von Mises statistic might be supported in future versions
  CvM <- FALSE
  
  ## set up return values
  partvar <- mf@get("part")
  m <- NCOL(partvar)
  pval <- rep.int(0, m)
  stat <- rep.int(0, m)
  ifac <- rep.int(FALSE, m)

  ## extract estimating functions  
  process <- as.matrix(estfun(obj))
  k <- NCOL(process)
  
  ## extract weights
  ww <- weights(obj)
  if(is.null(ww)) ww <- rep(1, NROW(process))
  n <- sum(ww)
  
  ## drop observations with zero weight
  ww0 <- (ww > 0)
  process <- process[ww0, , drop = FALSE]
  partvar <- partvar[ww0, , drop = FALSE]
  ww <- ww[ww0]
  ## repeat observations with weight > 1
  process <- process/ww
  ww1 <- rep.int(1:length(ww), ww)
  process <- process[ww1, , drop = FALSE]
  stopifnot(NROW(process) == n)

  ## scale process
  process <- process/sqrt(n)
  J12 <- root.matrix(crossprod(process))
  process <- t(chol2inv(chol(J12)) %*% t(process))  

  ## select parameters to test
  if(!is.null(parm)) process <- process[, parm, drop = FALSE]
  k <- NCOL(process)

  ## get critical values for CvM statistic
  if(CvM) {
    if(k > 25) k <- 25 #Z# also issue warning
    critval <- get("sc.meanL2")[as.character(k), ]
  } else {
    from <- if(trim > 1) trim else ceiling(n * trim)
    from <- max(from, minsplit)
    to <- n - from
    lambda <- ((n-from)*to)/(from*(n-to))

    beta <- get("sc.beta.sup")
    logp.supLM <- function(x, k, lambda)
    {
      if(k > 40) {
        ## use Estrella (2003) asymptotic approximation
        logp_estrella2003 <- function(x, k, lambda)
          -lgamma(k/2) + k/2 * log(x/2) - x/2 + log(abs(log(lambda) * (1 - k/x) + 2/x))
        ## FIXME: Estrella only works well for large enough x
	## hence require x > 1.5 * k for Estrella approximation and
	## use an ad hoc interpolation for larger p-values
	p <- ifelse(x <= 1.5 * k, (x/(1.5 * k))^sqrt(k) * logp_estrella2003(1.5 * k, k, lambda), logp_estrella2003(x, k, lambda))
      } else {
        ## use Hansen (1997) approximation
        m <- ncol(beta)-1
        if(lambda<1) tau <- lambda
        else tau <- 1/(1+sqrt(lambda))
        beta <- beta[(((k-1)*25 +1):(k*25)),]
        dummy <- beta[,(1:m)]%*%x^(0:(m-1))
        dummy <- dummy*(dummy>0)
        pp <- pchisq(dummy, beta[,(m+1)], lower.tail = FALSE, log.p = TRUE)
        if(tau==0.5)
          p <- pchisq(x, k, lower.tail = FALSE, log.p = TRUE)
        else if(tau <= 0.01)
          p <- pp[25]
        else if(tau >= 0.49)
          p <- log((exp(log(0.5-tau) + pp[1]) + exp(log(tau-0.49) + pchisq(x,k,lower.tail = FALSE, log.p = TRUE)))*100)
        else
        {
          taua <- (0.51-tau)*50
          tau1 <- floor(taua)
          p <- log(exp(log(tau1 + 1 - taua) + pp[tau1]) + exp(log(taua-tau1) + pp[tau1+1]))
        }
      }
      return(as.vector(p))
    }
  }

  ## compute statistic and p-value for each ordering
  for(i in 1:m) {
    pvi <- partvar[,i]
    pvi <- pvi[ww1]
    if(is.factor(pvi)) {
      proci <- process[ORDER(pvi), , drop = FALSE]
      ifac[i] <- TRUE

      # re-apply factor() added to drop unused levels
      pvi <- factor(pvi[ORDER(pvi)])
      # compute segment weights
      segweights <- as.vector(table(pvi))/n ## tapply(ww, pvi, sum)/n      

      # compute statistic only if at least two levels are left
      if(length(segweights) < 2) {
        stat[i] <- 0
	pval[i] <- NA
      } else {      
        stat[i] <- sum(sapply(1:k, function(j) (tapply(proci[,j], pvi, sum)^2)/segweights))
        pval[i] <- pchisq(stat[i], k*(length(levels(pvi))-1), log.p = TRUE, lower.tail = FALSE)
      }
    } else {
      oi <- if(breakties) {
        mm <- sort(unique(pvi))
	mm <- ifelse(length(mm) > 1, min(diff(mm))/10, 1)
	ORDER(pvi + runif(length(pvi), min = -mm, max = +mm))
      } else {
        ORDER(pvi)
      }    
      proci <- process[oi, , drop = FALSE]
      proci <- apply(proci, 2, cumsum)
      stat[i] <- if(CvM) sum((proci)^2)/n 
        else if(from < to) {
	  xx <- rowSums(proci^2)
	  xx <- xx[from:to]
	  tt <- (from:to)/n
	  max(xx/(tt * (1-tt)))	  
	} else {
	  0
	}
      pval[i] <- if(CvM) log(approx(c(0, critval), c(1, 1-as.numeric(names(critval))), stat[i], rule=2)$y)
        else if(from < to) logp.supLM(stat[i], k, lambda) else NA
    }
  }

  ## select variable with minimal p-value
  best <- which.min(pval)
  if(length(best) < 1) best <- NA
  rval <- list(pval = exp(pval), stat = stat, best = best)
  names(rval$pval) <- names(partvar)
  names(rval$stat) <- names(partvar)
  if (!all(is.na(rval$best)))
      names(rval$best) <- names(partvar)[rval$best]
  return(rval)
}

### split in variable x, either ordered or nominal
mob_fit_splitnode <- function(x, obj, mf, weights, minsplit, objfun, verbose = TRUE) {

    ## process minsplit (to minimal number of observations)
    if (minsplit > 0.5 & minsplit < 1) minsplit <- 1 - minsplit
    if (minsplit < 0.5)
        minsplit <- ceiling(sum(weights) * minsplit)
   
    if (is.numeric(x)) {
    ### for numerical variables
        ux <- sort(unique(x))
        if (length(ux) == 0) stop("cannot find admissible split point in x")
        dev <- vector(mode = "numeric", length = length(ux))

        for (i in 1:length(ux)) {
            xs <- x <= ux[i]
            if (mob_fit_checksplit(xs, weights, minsplit)) {
                dev[i] <- Inf
            } else {
                dev[i] <- mob_fit_getobjfun(obj, mf, weights, xs, objfun = objfun)
            }
        }

        ## maybe none of the possible splits is admissible
        if (all(!is.finite(dev))) return(FALSE)

        split <- list(variableID = NULL, ordered = TRUE, 
                      splitpoint = as.double(ux[which.min(dev)]),
                      splitstatistic = dev, toleft = TRUE)
        class(split) <- "orderedSplit"
    } else {
    ### for categorical variables
        al <- mob_fit_getlevels(x)
        dev <- apply(al, 1, function(w) {
                   xs <- x %in% levels(x)[w]
                   if (mob_fit_checksplit(xs, weights, minsplit)) {
                       return(Inf)
                   } else {
                       mob_fit_getobjfun(obj, mf, weights, xs, objfun = objfun)
                   }
               })

        if (verbose) {
            cat(paste("\nSplitting ", if(is.ordered(x)) "ordered ",
	              "factor variable, objective function: \n", sep = ""))
            print(dev)
        }

        if (all(!is.finite(dev))) return(FALSE)

        ## ordered factors are of storage mode "numeric" in party!
        ## initVariableFrame coerces ordered factors to storage.mode "numeric"
        ## the following is consistent with party
        
        if (is.ordered(x)) {
            split <- list(variableID = NULL, ordered = TRUE,
                          splitpoint = as.double(which.min(dev)),
                          splitstatistic = dev, toleft = TRUE)
            class(split) <- "orderedSplit"
            attr(split$splitpoint, "levels") <- levels(x)
        }  else {
            tab <- as.integer(table(x[weights > 0]) > 0)
            split <- list(variableID = NULL, ordered = FALSE,
                          splitpoint = as.integer(al[which.min(dev),]), 
                          splitstatistic = dev, 
                          toleft = TRUE, table = tab)
            attr(split$splitpoint, "levels") <- levels(x)
            class(split) <- "nominalSplit"
        }
    }
    split
}

### get partitioned objective function for a particular split
mob_fit_getobjfun <- function(obj, mf, weights, left, objfun = deviance) {
  ## mf is the model frame
  ## weights are the observation weights
  ## left is 1 (if left of splitpoint) or 0
  weightsleft <- weights * left
  weightsright <- weights * (1 - left)

  ### fit left / right model 
  fmleft <- reweight(obj, weights = weightsleft)
  fmright <- reweight(obj, weights = weightsright)

  return(objfun(fmleft) + objfun(fmright))
}

### determine all possible splits for a factor, both nominal and ordinal
mob_fit_getlevels <- function(x) {
    nl <- nlevels(x)
    if (inherits(x, "ordered")) {
        indx <- diag(nl)
        indx[lower.tri(indx)] <- 1
        indx <- indx[-nl,]
	rownames(indx) <- levels(x)[-nl]
    } else {
        mi <- 2^(nl - 1) - 1
        indx <- matrix(0, nrow = mi, ncol = nl)
        for (i in 1:mi) { # go though all splits #
            ii <- i
            for (l in 1:nl) {
                indx[i, l] <- ii%%2;
                ii <- ii %/% 2   
            }
        }
        rownames(indx) <- apply(indx, 1, function(z) paste(levels(x)[z > 0], collapse = "+"))
    }
    colnames(indx) <- as.character(levels(x))
    storage.mode(indx) <- "logical"
    indx
}

### check split
mob_fit_checksplit <- function(split, weights, minsplit)
    (sum(split * weights) < minsplit || sum((1 - split) * weights) < minsplit)
