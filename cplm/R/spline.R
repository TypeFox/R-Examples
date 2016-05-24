###########################################################
# Functions used in fitting additive models in cpglmm     #
###########################################################


#######################
# get model frame and factor list 
# for cpglmm with smoothing terms
#######################
frFL <- function (formula, data, family, control = list(),  
    verbose, weights, offset, contrasts, basisGenerators, bySetToZero = T) 
{
    call <- match.call()
    formula <- eval(call$formula)
    tf <- terms.formula(formula, specials = eval(call$basisGenerators, 
        parent.frame(2)))
    f.ind <- unlist(attr(tf, "specials"))
    n.f <- length(f.ind)
    rhs <- safeDeparse(formula[[3]])
    fctterm <- fct <- vector(mode = "list", length = n.f)
    for (i in 1:n.f) 
      fctterm[[i]] <- attr(tf, "variables")[[f.ind[i] + 1]]
    fct <- lapply(fctterm, eval, envir = data, enclos = parent.frame(2))
    for (i in seq_along(fct)) 
      fct[[i]] <- expandBasis(fct[[i]], eval(attr(fct[[i]], "call")$by, data), 
                  eval(attr(fct[[i]], "call")$varying, data), bySetToZero)
    names(fct) <- names(fctterm) <- paste("f.", lapply(fct, 
            function(x) {
                paste(as.character(attr(x, "call")$x), ifelse(!is.null(eval(attr(x, 
                  "call")$varying, data)), paste("X", deparse(attr(x, 
                  "call")$varying), sep = ""), ""), ifelse(eval(attr(x, 
                  "call")$allPen), paste(".", deparse(attr(x, 
                  "call")$by), sep = ""), ""), sep = "")
            }), sep = "")
    rhs <- subFcts(rhs, fctterm, fct, data)
    data <- expandMf(data, fct)
    call[[1]] <- as.name("lmer")
    call$doFit <- FALSE
    call$data <- as.name("data")
    call$formula <- as.formula(paste(formula[[2]], "~", rhs))
    call["basisGenerators"] <- NULL
    m <- eval(call, data)
    #    m$fr$mf <- data
    m <- subAZ(m, fct)
    fctterm <- lapply(fct, function(x) attr(x, "call"))
    return(list(m = m, fct = fct, fctterm = fctterm))
}


### The main event
lmer <-
  function(formula, data, family = NULL, REML = TRUE,
           control = list(), start = NULL, verbose = FALSE, doFit = TRUE,
           subset, weights, na.action, offset, contrasts = NULL,
           model = TRUE, x = TRUE, ...)
    ### Linear Mixed-Effects in R
  {
    mc <- match.call()
    if (!is.null(family)) {             # call glmer
      mc[[1]] <- as.name("glmer")
      return(eval.parent(mc))
    }
    stopifnot(length(formula <- as.formula(formula)) == 3)
    
    fr <- lmerFrames(mc, formula, contrasts) # model frame, X, etc.
    FL <- lmerFactorList(formula, fr, rmInt=FALSE, drop=FALSE) # flist, Zt, dims
    largs <- list(...)
    if (!is.null(method <- largs$method)) {
      warning(paste("Argument", sQuote("method"),
                    "is deprecated.  Use", sQuote("REML"),
                    "instead"))
      REML <- match.arg(method, c("REML", "ML")) == "REML"
      largs <- largs[names(largs) != "method"]
    }
    ### FIXME: issue a warning if the control argument has an msVerbose component    
    FL$dims["LMM"] <- 1L
    FL$dims["mxit"] <- 500L
    FL$dims["mxfn"] <- 1000L
    ans <- list(fr = fr, FL = FL, start = start)
    ans
  }

############################
# function for 2d splines (from Ngo and Wand)
############################
# compute thin plate spline covariance function  
tps.cov <- function(r) {
  r <- as.matrix(r)
  num.row <- nrow(r)
  num.col <- ncol(r)
  r <- as.vector(r)
  nzi <- (1:length(r))[r!=0]
  ans <- rep(0,length(r))
  ans[nzi] <- r[nzi]*r[nzi]*log(abs(r[nzi]))
  if (num.col>1) 
    ans <- matrix(ans,num.row,num.col)
  return(ans)
}

sp2d <- function(x1, x2, k = max(20,min(length(x1)/4,150)), 
                 by = NULL, allPen = FALSE, varying = NULL, 
                 diag = FALSE, knots1 = quantile(x1, probs = 1:k/(k+1)),
                 knots2 = quantile(x1, probs = 1:k/(k+1))) {
  call <- as.list(expand.call(match.call()))
  knots1 <- eval(knots1)
  knots2 <- eval(knots2)
  k <- eval(k)
  # design matrix for fixed effects
  X <- cbind(x1,x2)
  knots <- cbind(knots1,knots2)
  dist.mat <- matrix(0,k,k)
  dist.mat[lower.tri(dist.mat)] <- dist(knots)
  dist.mat <- dist.mat + t(dist.mat)
  Omega <- tps.cov(dist.mat)
  diffs.1 <- outer(x1,knots1,"-")
  diffs.2 <- outer(x2,knots2,"-")
  dists <- sqrt(diffs.1*diffs.1+diffs.2*diffs.2)
  svd.Omega <- svd(Omega)
  sqrt.Omega <- t(svd.Omega$v %*% (t(svd.Omega$u) * sqrt(svd.Omega$d)))
  # design matrix for random effects
  Z <- t(solve(sqrt.Omega,t(tps.cov(dists))))
  res <- list(X=X, Z=Z, knots=knots)
  attr(res, "call") <- as.call(call)  
  return(res)
}


###########################################################
# These are functions copied from the "amer" package      #
###########################################################

getF <- function (object, which, n=100, newdata=NULL, interval = c("NONE", "MCMC", "RW"), addConst = TRUE, varying=1, 
                  level = 0.9, sims = 1000)
{
  
  stopifnot(inherits(object, "cpglmm"), is.null(newdata)||is.data.frame(newdata), n>0, sims>0, level < 1, level > 0)
  
  
  terms <- object@smooths
  interval <- toupper(interval)
  interval <- match.arg(interval)
  
  n <- as.integer(n)
  sims <- as.integer(sims)
  if(missing(which)) which <- seq_along(terms)
  if(is.character(which)) {
    which <- match(which, names(terms))
    if(any(nas <- is.na(which))) 
      warning("entry ", paste(which[nas], collapse=", "), " in 'which' did not match any function names in ", safeDeparse(object) ,".")
    which <- which[!nas]
  }
  if(length(addConst) != length(which)) addConst <- rep(addConst, length=length(which))
    
  
  ci.RW <- function(fhat, base, terms, i, j, object, level, addConst){
    fctV <-
      function(m, indGrp, indPen, indUnpen){
        #Cov(hat.beta, hat.b-b) for bias-adjusted empirical Bayes CIs (s. Ruppert/Wand(2003), Semiparametric Regression, p. 138 f.):
        #use V = cov(hat.fixef, hat.ranef) = sigma.eps^2 (C'C + sigma.eps^2/sigma.b^2 D)^-1; C=[XZ]*sqrt(W), D = blockdiag(0, I_dim(b))
        
        C <- cBind(m@X[,indUnpen, drop=F], t(as.matrix(m@Zt[indPen,])))
        if(length(m@var)) C <- C * sqrt(1/m@var)
        
        V <- crossprod(C)
        # FIXME: this works only for scalar ST and length(indGrp=1)- don't use if allPen=T! 
        if(m@ST[[indGrp]] > 0){
          diag(V)[-(1:length(indUnpen))] <- diag(V)[-(1:length(indUnpen))] + m@ST[[indGrp]]^-2
        } else {
          #FIXME: HACK: V is not invertible for var(ranef)=0, use var(ranef)=10^-9 instead 
          diag(V)[-(1:length(indUnpen))] <- diag(V)[-(1:length(indUnpen))] + 10^9
        }
        return(sigma(m)^2 * solve(V))
      }
    
    z <- qnorm(1-(1-level)/2)
    indUnpen <- if(addConst){
      c(attr(terms[[i]],"indConst")[[j]], attr(terms[[i]],"indUnpen")[[j]])
    } else attr(terms[[i]],"indUnpen")[[j]]	
    cV <- as(chol(fctV(object, attr(terms[[i]],"indGrp")[[j]], 
                       attr(terms[[i]],"indPen")[[j]], indUnpen)), "sparseMatrix")
    C <- as(cBind(base$X[[j]], base$Z[[j]]), "sparseMatrix")
    sd <- apply(C, 1, function(x, cV){
      ctc <- cV %*% as(x, "sparseMatrix")
      return(sqrt(sum(ctc * ctc)))
      #return(sqrt((crossprod(cV%*%as(x, "sparseMatrix")))@x))	
    }, cV=cV)
    
    ci <- cbind(fhat - z*sd, fhat + z*sd)
    colnames(ci) <- c("lo", "hi")
    return(ci)
  }
  
  ans <- vector(mode="list", length=length(which))
  if(interval=="MCMC") attr(ans, "mcmc") <- mcmc 
  names(ans) <- names(terms)[which]
  indWhich <- 1
  
  for(i in which){
    #################################
    # set up / check newdata
    ################################
    if(!is.null(terms[[i]]$by)){
      lvls <- levels(object@frame[, safeDeparse(terms[[i]]$by)])#FIXME: in amerSetup: this will fail if terms[[i]]$x was only in the workspace but not in the supplied data.frame for the original call.
      hasBy <- TRUE
    } else hasBy <-FALSE	
    hasVarying <- !is.null(terms[[i]]$varying)
    
    if(is.null(newdata)){
      grid <- TRUE
      #FIXME: adapt this for 2d/3d-smooths
      #get range of covariates + sequence of values
      lim <- range(object@frame[, safeDeparse(terms[[i]]$x)], na.rm=T) #FIXME: in amerSetup: this will fail if terms[[i]]$x was only in the workspace but not in the supplied data.frame for the original call.
      newX <- seq(lim[1], lim[2], l=n)
      if(hasBy){
        newBy <- factor(rep(lvls, length=n), labels=lvls)
        data <- data.frame(newX, newBy)
        colnames(data) <- c(safeDeparse(terms[[i]]$x), safeDeparse(terms[[i]]$by))
      } else {
        data <- data.frame(newX)
        colnames(data) <- safeDeparse(terms[[i]]$x)
      }
      if(hasVarying){
        #varying covariate is set value of varying
        data <- cbind(data, rep(varying, nrow(data)))
        colnames(data)[NCOL(data)] <- safeDeparse(terms[[i]]$varying)
      }
    } else {
      grid <- FALSE
      data <- newdata
      n <- nrow(newdata)
      vnames <- safeDeparse(terms[[i]]$x)
      if(hasBy) vnames <- c(vnames,safeDeparse(terms[[i]]$by))
      if(hasVarying) vnames <- c(vnames, safeDeparse(terms[[i]]$varying))
      if(any(nas <- is.na(match(vnames, colnames(data))))) 
        stop("variable ", paste(vnames[nas], collapse=", "), "not found in given data.")
      data <- data[, colnames(data) %in% vnames, drop=F]
    }
    
    #################################
    #create basis:
    #################################
    base <- eval(terms[[i]], data)
    base <- expandBasis(base, 
                        by = eval(attr(base, "call")$by, data),  
                        varying = eval(attr(base, "call")$varying, data),
                        bySetToZero = !grid)
    
    # need to modify Z, X for allPen-Fits
    if(eval(terms[[i]]$allPen)){
      nlvl <- length(lvls)
      base0 <- base
      #where are the random effects for the penalized spline functions:
      ##use the first because the spline will have more levels (grps*(p-d)) than the grouping factor (grps)
      indZ <- reinds(object@Gp)[[ attr(terms[[i]],"indGrp") [[1]] [1] ]]  
      #how many penalized spline functions per level of by
      dimOneZ <- length(indZ)/length(lvls)
      useZ <- 1:dimOneZ
      if(grid) fullZ <- expandBasis(eval(terms[[i]], data), 
                                    by = NULL,  
                                    varying = eval(attr(base, "call")$varying, data),
                                    bySetToZero = FALSE)$Z
    }	
    
    nf <- ifelse(hasBy, length(lvls), 1)
    ans[[indWhich]] <- vector(mode="list", length = nf)
    names(ans[[indWhich]]) <-  if(hasBy) paste(safeDeparse(terms[[i]]$by), lvls, sep="") else names(base$X)
    for(j in seq_along(ans[[indWhich]])){
      #################################
      #calculate fits and cis
      #################################	
      
      if(eval(terms[[i]]$allPen)){
        ansInd <- 1
        #need to:
        #-append to Z extra columns for the random effects for base$X	 (+ a random intercept for the by-levels)
        #-set columns in base$Z not relevant for the current by-level to 0
        #-set base$X to zero
        
        base$Z[[1]] <- base0$Z[[1]]
        if(grid){
          #fill up rows having artifical zeroes because of structure of newBy with values
          base$Z[[1]][,useZ] <- fullZ[[1]]
        }	
        base$Z[[1]][,-useZ] <- 0
        
        #step to next block:
        useZ <- useZ + dimOneZ
        lvlInd <- rep(0, nlvl)
        lvlInd[j] <- 1
        
        X <- cBind("(Intercept)"=1, base0$X[[1]])
        if(!grid){
          use <- eval(terms[[i]]$by, data)==lvls[j]
          X[!use,] <- 0
        } 
        if(eval(terms[[i]]$diag)){
          #lmer switches order of terms in X, need to permute X-columns accordingly: 
          uNames <-  unlist(unique(sapply(object@ST[attr(terms[[i]],"indGrp")[[1]][-1]], dimnames)))
          X <- X[,uNames]
        }
        
        base$Z[[1]] <-  as(cBind(base$Z[[1]], kronecker(X, t(lvlInd), FUN = "*")),"sparseMatrix")
        base$X[[1]] <- matrix(0, nrow=n, ncol=0)
      } else {
        ansInd <- j
        #if(!grid && hasBy){
        #        #remove unnecessary rows from design
        #        use <- eval(terms[[i]]$by, data)==lvls[j]
        #        base$X[[ansInd]] <- base$X[[ansInd]][use,,drop=F]
        #        base$Z[[ansInd]] <- base$Z[[ansInd]][use,,drop=F]
        #}
      }	
      
      if(addConst[i]){
        #add columns for constant terms to X, append indUnpen:	
        byColumn <-if(hasBy && paste(safeDeparse(terms[[i]]$by),lvls[j],sep="") %in% names(object@fixef)){
          rep(1, nrow(base$X[[ansInd]]))
        } else numeric(0) 
        base$X[[ansInd]] <- cBind(byColumn, base$X[[ansInd]])	
        if("(Intercept)" %in% names(object@fixef)[attr(terms[[i]],"indConst")[[ansInd]]])
          base$X[[ansInd]] <- cBind(1,base$X[[ansInd]])
        
        indUnpen <- c(attr(terms[[i]],"indConst")[[ansInd]], attr(terms[[i]],"indUnpen")[[ansInd]])
      } else indUnpen <- attr(terms[[i]],"indUnpen")[[ansInd]]	
      
      fhat <- base$X[[ansInd]] %*% 
        object@fixef[indUnpen, drop = F] + 
        base$Z[[ansInd]] %*% 
        object@ranef[attr(terms[[i]],"indPen")[[ansInd]], drop = F] 
      
      if(!eval(terms[[i]]$allPen)){ 
        ci <- switch(interval,
                     "NONE" = matrix(NA, nrow=nrow(base$X[[ansInd]]), ncol=0),
                     "RW" = ci.RW(fhat, base, terms, i, j, object, level, addConst[i]))
      } else {
        #TODO: implement CIs for fits with allPen = T
        if(interval!="NONE" && j==1) warning("CIs for fits with allPen = T not yet implemented.")
        ci <- matrix(NA, nrow=nrow(base$X[[ansInd]]), ncol=0)
      }
      dataJ <- if(grid){
        data[,!(colnames(data)==safeDeparse(terms[[i]]$by)), drop=F]
      } else {
        ## if(!grid && hasBy){
        ##     data[use,] 
        ## } else 
        data  
      }	
      ans[[indWhich]][[j]] <- data.frame(dataJ, fhat= as.matrix(fhat), ci)
    }# end for j
    indWhich <- indWhich + 1
  }#end for i
  return(ans)
}



plotF <- function(object, which, n=100, interval = "RW", addConst = TRUE, trans=I,  
                  level = 0.9, sims = 1000, auto.layout = TRUE, rug = TRUE, legendPos="topright", ...)
{
  # FIXME: add option for centering function estimates at zero/ anchoring at mean of all other covariates?
  # FIXME: valid confints for trans? 	
  terms <- object@smooths
  if(missing(which)) which <- seq_along(terms)
  if(is.character(which)) {
    which <- match(which, names(terms))
    if(any(nas <- is.na(which))) 
      warning("entry ", paste(which[nas], collapse=", "), " in 'which' did not match any function names in ", safeDeparse(object) ,".")
    which <- which[!nas]
  }
  if(length(legendPos) != length(which)) legendPos <- rep(legendPos, length=length(which))
  if(length(addConst) != length(which)) addConst <- rep(addConst, length=length(which))
  
  interval <- toupper(interval)
  
  res <- getF(object = object, which = which, n = n, newdata=NULL, interval = interval, addConst = addConst, level = level, sims = sims)
  allPen <- sapply(object@smooths[which], function(x) eval(x$allPen))
  
  if(auto.layout){
    oldpar <- NULL
    on.exit(par(oldpar))
    nf <- length(res)
    oldpar <- par(mfrow = set.mfrow(Nparms=nf))
  }
  
  
  plot1F <- function(x, interval, legendPos, ...){
    dots <- list(...)
    fhat <- sapply(x, function(x){trans(x$fhat)})
    if(interval) ci <- lapply(x, function(x){trans(cbind(x$lo, x$hi))})
    cov <- sapply(x, function(x){x[,1]})
    if(is.null(dots$ylim)){
      ylim <- range(fhat)
      if(interval) ylim <- range(ylim, ci)
    } else {
      ylim <- dots$ylim
      dots <- dots[names(dots)!="ylim"]
    }	
    if(is.null(dots$xlim)){
      xlim <- range(cov)
    } else {
      xlim <- dots$xlim
      dots <- dots[names(dots)!="xlim"]
    }
    do.call(plot, c(list(x=-2*xlim[1]-10, y=0, ylim = ylim, xlim = xlim, ylab="", xlab=colnames(x[[1]])[1]), dots))
    do.call(matlines, c(list(x = cov, y = fhat, col=1:length(x), lty=1), dots))
    if(interval) for(i in 1:length(ci)) do.call(matlines, c(list(x=cov[,i], y=ci[[i]], col=i, lty=2),dots))
    if(legendPos != "none" && length(x)>1){
      do.call(legend,c(list(x=legendPos, legend=names(x), lty=1, col=1:length(x)),dots))
    }	
  }
  dots <- list(...)
  for(i in seq_along(res)){
    #FIXME: change as CIs for allPen become available
    plot1F(res[[i]], interval = !(interval=="NONE")&&!allPen[i], legendPos = legendPos[i], ...)
    if(is.null(dots$ylab)){
      ylab <- ifelse(addConst[i], paste(names(res)[i], "+ const"), names(res)[i])
      if(any(trans(-2:2)!= (-2:2))) ylab <- paste(safeDeparse(match.call()$trans),"(",ylab,")",sep="")
    } else{
      ylab <-dots$ylab
      dots <- dots[names(dots)!="ylab"]
    }	
    do.call(title, c(list(ylab= ylab), dots))
    if(rug){
      if(length(res[[i]])==1){
        rug(object@frame[,colnames(res[[i]][[1]])[1]], ...)
      } else {
        nlvls <- length(res[[i]]) 
        lvls <- levels(object@frame[, safeDeparse(object@smooths[[i]]$by)])
        for(j in 1:nlvls){
          use <- object@frame[,safeDeparse(object@smooths[[i]]$by)] == lvls[j]
          rug(object@frame[use, colnames(res[[i]][[1]])[1]], col =j, ...)
        }	
      }	
    } 
  }
  invisible(res)
}
