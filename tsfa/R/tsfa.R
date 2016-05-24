#  Notes
#   An  FAmodel has parameters but not factors. It can
#     optionally have stats (about estimation).
#   An  fFAmodel extends an  FAmodel by adding factors.
#   An  TSFmodel extends an  fFAmodel by time.
###################################################

#            Factor Analysis 

###################################################

FAmodel <- function(obj, ...)UseMethod("FAmodel")
FAmodel.FAmodel <- function(obj, ...) obj #extractor

FAmodel.default <- function(obj,  Omega=NULL, Phi=NULL, LB=NULL, LB.std=NULL,
   stats=NULL,  ...)
  {#  obj should be the loadings (hat)B for  x =  B f + e  
   if(!is.matrix(obj))
     stop("FAmodel.default requires a loadings matrix (factor loadings) as the first argument.")
   # do more checking (but only loadings is necessary (for simulation)
   #if(ncol(obj)!= nrow(Omega))
   #    stop("dimensions of obj (loadings) and  Omega do not agree.")
   classed(list(loadings=obj, Omega=Omega, Phi=Phi, LB=LB,  LB.std=LB.std, 
	   stats=stats), "FAmodel") # constructor
  }


factors <- function(x)UseMethod("factors")
# use predict with data to get factors
factors.fFAmodel <- function(x) x$f

nfactors <- function(x) UseMethod("nfactors")
nfactors.FAmodel <- function(x) {ncol(loadings(x))}

factorNames <- function(x) UseMethod("factorNames")
factorNames.FAmodel     <- function(x) dimnames(loadings(x))[[2]]

coef.FAmodel <- function(object, ...) {c(loadings(object), diag(object$Omega))}

explained <- function(object, ...)UseMethod("explained")
explained.FAmodel <- function (object, f=factors(object),
                  names=dimnames(loadings(object))[[1]], ...) {
  # portion of data explained by factors
  r <- t(loadings(object) %*% t(f))
  dimnames(r) <- list(NULL, names)
  r
  }

LedermannBound  <- function(M) {
   if (is.matrix(M)) M <- ncol(M)
   if (1 != length(M)) stop("M must be an integer number of indicator variables.")
      
   ## solve (M^2-M) - (2M+1)k + k^2 = 0 for k
   #r <- polyroot(c(M^2-M, -(2*M+1), 1))
   #if (any(abs(Im(r)) > 10*.Machine$double.eps))
   #         warning("LedermannBound has complex solution")
   #r <- Re(r)
   #r[(0 <=r) & (r <= M)]
   if(M <3) return(0)
   r <- M + 0.5 - (0.5 * sqrt(8*M + 1))
   #fuzz in next is just to assure rounding errror does not give an non-integer
   # when result should be integer. (Especially important if floor might result
   # in the next lower integer.)
   if (1e-10 > abs(r - round(r))) round(r) else r
   }

# FAfitStats could be moved here, but needs to be reorganized to use
#   cov rather than data (as in estFAmodel)

# could add summaryStats 

summary.FAmodel <- function(object, ...)
 {classed(list(k=nfactors(object), M=nrow(loadings(object)),
      Snames=dimnames(loadings(object))[[1]], Fnames=factorNames(object),
      Omega=!is.null(object$Omega), 
      Phi=!is.null(object$Phi), 
      LB=!is.null(object$LB), 
      estConverged=object$stats$estConverged, 
      rotationConverged=object$stats$rotationConverged, 
      orthogonal=object$stats$orthogonal
      ), "summary.TSFmodel")
 }

print.summary.FAmodel <- function (x, ...)
  {cat(x$k, " factors: ",     x$Fnames, "\n")
   cat(x$M, " indicators: ",  x$Snames, "\n")
   cat("Omega ", (if(x$Omega) "is" else "is not"), " specified.\n")
   cat("Phi   ", (if(x$Phi) "is" else "is not"), " specified.\n")
   cat("LB    ", (if(x$LB) "is" else "is not"), " specified.\n")
   if(!is.null(x$estConverged)) cat("loadings estimation ",
           (if(x$estConverged) "converged.\n" else "did not converge.\n"))
   if(!is.null(x$rotationConverged)) cat("rotation ",
           (if(x$rotationConverged) "converged.\n" else "did not converge.\n"))
   if(!is.null(x$orthogonal)) cat("rotation is ",
           (if(x$orthogonal) "orthogonal.\n" else "oblique.\n"))
  }


predict.FAmodel <- function(object, data = NULL, factorNames.=factorNames(object), ...){
      # prediction of factors with data
      if (is.null(data)) stop("data must be supplied.")
      r <- data %*% t(object$LB) #hatf
      dimnames(r) <- list(NULL, factorNames.)
      r
      }

permusign <- function(B, Btarget, Phi=NULL) {

# Selects the permutation and signs of the columns of the factor loadings B
# that resembles the Btarget matrix the most.
# Phi matrix (cov of factors) may need to be reordered for the permutation.
  
  ############################## local functions
  permute <- function(x){
     # with thanks to Bill Venables
     if (length(x) <= 1) as.matrix(x) else{
  	 M <- NULL
  	 for (i in seq(length(x))) M <- rbind(M, cbind(x[i], Recall(x[-i])))
  	 M
  	 }
     }
  signsw <- function(Bprop, Bnew, Btarget){
     # compare distance of Bprop from Btarget and also Bprop with column 
     # signs switched. If the best of these is better than Bnew to Btarget
     #  return the column signs (1 or -1), otherwise return NULL
     signs <- rep(1, ncol(Bprop))
     d1 <- colSums((Btarget - Bprop)^2)
     d2 <- colSums((Btarget + Bprop)^2) # all col signs reversed
     if ( sum(pmin(d1, d2)) < sum((Btarget - Bprop %*% diag(signs))^2))
     	signs <- 2 * (d1 < d2) - 1
     # the fuzz (1e-12) seems to be necessary to avoid rounding error causing
     # T when things should be equal,with the result that random 
     #  permutations occur.
     if (sum((Btarget - Bnew)^2) > 1e-12 +
         sum((Btarget - Bprop %*% diag(signs))^2) ) signs else NULL
     }

  ############################## end local functions

  P    <- permute(seq(ncol(B))) # permutation matrix
  Bnew   <- B
  PhiNew <- Phi
  if( ! is.null(Phi)) {
    for (j in seq(nrow(P))) {
      Bprop <- B[,P[j,]]
      signs <- signsw(Bprop, Bnew, Btarget)
      if(!is.null(signs)){
  	 #cat(j, ":", P[j,], signs)
  	 Bnew	<-  Bprop %*% diag(signs)
  	 PhiNew <-  (Phi[P[j,],P[j,]]) * outer(signs,signs)
  	 }
      }
    }
  list(loadings=Bnew,Phi=PhiNew)
}


estFAmodel <- function(Sigma, p, n.obs=NA,
                est="factanal", 
		estArgs=list(scores="none", control=list(opt=list(maxit=10000))),
		rotation=if(p==1) "none" else "quartimin", rotationArgs=NULL,
		GPFargs=list(Tmat=diag(p), normalize=TRUE, eps=1e-5, maxit=1000),
		BpermuteTarget=NULL,
                factorNames=paste("Factor", seq(p)),
                indicatorNames=NULL) {
  if (1e10 < max(diag(Sigma))/ min(diag(Sigma)) ) warning(
	 "Data variances are very different. Consider rescaling some indicators.")
 
  stds <- sqrt(diag(Sigma))
  if(p < 0)  stop("p (number of factors) must be greater than 0,")
  else if(p == 0) {
      # zero factors
      uniquenesses <- diag(1, nrow(Sigma))
      Omega  <- diag(Sigma)
      if(est != "factanal") 
        stop("Currently Omega  is only correct for factanal estimation.")
      loadings.std <- loadings <- Phi  <- LB <- LB.std <- NULL
      estConverged <- rotationConverged <- orthogonal <- TRUE
      }
  else {
      z <- do.call(est, c(list(covmat = Sigma, 
                 factors=p, n.obs=n.obs, rotation="none"), estArgs))

      estConverged <- z$converged
      
      uniquenesses <- z$uniquenesses
       
      # for debugging compare:  hatOmega - Omega, hatOmega, Omega

      # z$loadings is orth solution
      if (rotation == "none") {
         loadings.std <- z$loadings              
	 Phi  <- NULL
	 rotationConverged <- TRUE
	 orthogonal <- TRUE
 	 }
      else {	 
  	 rotB <- do.call(rotation, c(list(z$loadings), GPFargs, rotationArgs))
    	 loadings.std <-  rotB$loadings
	 rotationConverged <- rotB$convergence
	 orthogonal <- rotB$orthogonal
     
    	 # Make sure columns are ordered properly and have the correct signs.
    	 Phi <- rotB$Phi
    	 if (! is.null(BpermuteTarget)) {
    	    z  <- permusign(diag(stds) %*% loadings.std, BpermuteTarget, Phi=Phi)
    	    loadings.std <- diag(1/stds) %*% z$loadings
    	    Phi  <- z$Phi
    	    }
	 }
      loadings <- diag(stds) %*% loadings.std
      dimnames(loadings) <- list(indicatorNames, factorNames)

      ### Compute Bartlett-predicted factor scores 
      Sigma.std <- if(is.null(Phi)) loadings.std %*% t(loadings.std) + diag(uniquenesses) else
                        loadings.std %*% Phi %*% t(loadings.std) + diag(uniquenesses)
      SinvB <- solve(Sigma.std, loadings.std) 
      LB.std   <- solve(crossprod(loadings.std, SinvB), t(SinvB))

      LB <- LB.std %*% diag(1/stds)
      dimnames(LB) <- list(factorNames, indicatorNames)
      }
  FAmodel(loadings, Omega=diag(stds * uniquenesses * stds), 
      Phi=Phi, LB=LB, LB.std=LB.std,  
      stats=list(estConverged=estConverged, rotationConverged=rotationConverged,
	orthogonal=orthogonal, uniquenesses=uniquenesses, call=match.call()))
  }

###################################################

#        Time Series Factor Analysis 

###################################################


DstandardizedLoadings <- function(x)UseMethod("DstandardizedLoadings")
DstandardizedLoadings.TSFmodel <- function(x){
    r <- diag(1/sqrt(diag(cov(diff(x$data))))) %*% loadings(x)
    dimnames(r) <- dimnames(loadings(x))
    r
    }


# standardizedLoadings <- function(x)UseMethod("standardizedLoadings")
# standardizedLoadings.TSFestModel <- function(x){
#   To transform everything back to
#   the undifferenced scales and then standardize, then somehow a
#   pseudo-Phi and pseudo-Omega (called \Gamma_t and \Psi_t in the TSFA
#   paper) must be computed first.    
#     r <- diag(1/sqrt(diag(cov(x$data)))) %*% loadings(x)
#     dimnames(r) <- dimnames(loadings(x))
#     r
#     }

TSFmodel <- function(obj, ...)UseMethod("TSFmodel")
TSFmodel.TSFmodel <- function(obj, ...) obj #extractor

TSFmodel.default <- function(obj, f=NULL, Omega=NULL, Phi=NULL, LB=NULL,
        positive.data=FALSE, names=NULL, ...)
  {#  arg break.points=NULL, not yet supported
   #  obj should be the loadings (hat)B
   #     x =  B f + e  # NB no mean added to give x
   #   vector processes x, f, and e have times series matrices of data so
   #   calculation is eg t(B %*% t(f))

   if(is.null(f)) stop(" f must be specified.")

   if(!is.matrix(obj))
     stop("TSFmodel.default requires a loadings matrix (factor loadings) as the first argument.")

   if(ncol(obj)!=nseries(f))  stop("dimensions of obj (B) and  f do not agree.")
   if(is.null(names)) names <- paste("Series", seq(nrow(obj)))
   classed(list(loadings=obj, f=f, Omega=Omega, Phi=Phi, LB=LB, 
	   positive.data=positive.data,
           names=names,  #break.points=break.points, 
	   dots=list(...)), c("TSFmodel", "fFAmodel", "FAmodel")) # constructor
  }

TSFmodel.FAmodel <- function(obj, f=NULL, positive.data=FALSE, names=NULL, ...)
  {if(is.null(f)) stop(" f must be specified.")

   if(ncol(loadings(obj))!=nseries(f))  
       stop("dimensions of obj loadings and  f do not agree.")
   if(is.null(names)) names <- paste("Series", seq(nrow(obj$loadings)))
   classed(append(obj, list(f=f, positive.data=positive.data, names=names,...)),
          c("TSFmodel", "fFAmodel", "FAmodel")) # constructor
  }

simulate.TSFmodel <- function(model, f=factors(model), Cov=model$Omega, sd=NULL, noise=NULL,  
	rng=NULL, noise.model=NULL, ...)
   {#  (... further arguments, currently disregarded)
    # tframe and Tobs are taken from factors (f) 
    
    if ( is.null(Cov) & is.null(sd) & is.null(noise) & is.null(noise.model))
      stop("One of Cov, sd, noise, or noise.model, must be specified.")

    p <- nrow(loadings(model))
    noise <- dse::makeTSnoise(Tobs(f), p, 1, noise=noise, rng=rng,
                        Cov=Cov, sd=sd, noise.model=noise.model)

    #use the calculation in explained but discard the class setting
    x <- unclass(explained(model, f=f)) + noise$w  
    if (model$positive.data && any( x < 0 )) {
        warning("negative simulated data values set to zero.")
        x[ x < 0 ] <- 0
        }
    attr(x, "noise") <- noise
    attr(x, "TSFmodel") <- model
    tframed(x, tf=tframe(f), names=seriesNames(model))
}


factors.TSFmodel <- function(x) classed(x$f, c("TSFfactors", class(x$f)))
#factors.TSFestModel <- function(x) {
#   r <- factors(TSFmodel(x))
#   trueM <- attr(x$data, "TSFmodel")
#   if(!is.null(trueM)) attr(r, "true") <- factors(trueM)
#   r
#   }


factors.EstEval  <- function(x)
   {N <- length(x$result)
    r <- vector(mode="list", length=N)
    for (i in 1:N) r[[i]] <- factors(x$result[[i]])
    classed(list(result=r, truth=factors(x$truth)),
            c("factorsEstEval", "EstEval"))
  }

#diff.TSFmodel  <- function (x, ...){
#  x$f <- diff(x$f)
#  x 
#  }

# was diff.TSFestModel
diff.TSFmodel  <- function (x, ...){
  x$f <- diff(x$f)
  x$data <- diff(x$data)
  trueM <- attr(x$data, "TSFmodel")
  if(!is.null(trueM)){
     trueM$f <- diff(factors(trueM))
     attr(x$data, "TSFmodel") <- trueM
     }
  x 
  }

diff.TSFexplained  <- function (x, ...){
  tf <- diff(tframe(x))
  r  <- tframed(diff(unclass(x)), tf=tf)
  d  <- attr(x, "data")
  if(!is.null(d)) attr(r, "data") <- tframed(diff(d), tf=tf)
  classed(r, c("TSFexplained", class(r))) 
  }

diff.TSFfactors  <- function (x, ...){
  tf <- diff(tframe(x))
  r <- tframed(diff(unclass(x)), tf=tf)
  truef <- attr(x, "true")
  if(!is.null(truef)) attr(r, "true") <- tframed(diff(unclass(truef)), tf=tf)
  classed(r, c("TSFfactors", class(r))) 
  }

diff.factorsEstEval  <- function (x, ...){
  N <- length(x$result)
    r <- vector(mode="list", length=N)
    for (i in 1:N) r[[i]] <- diff(x$result[[i]])
    classed(list(result=r, truth=diff(x$truth)),
            c("factorsEstEval", "EstEval"))
  }

nfactors.EstEval <- function(x) {nfactors(x$truth)}
nfactors.TSFfactors <- function(x) {nseries(x)}

seriesNames.TSFmodel <- function(x) {x$names}

factorNames.EstEval     <- function(x) factorNames(x$truth)
factorNames.TSFfactors     <- function(x) seriesNames(x)


tframe.TSFmodel <- function(x) tframe(x$f)

tfplot.TSFmodel <- function(x,..., tf=tfspan(x , ...),   
      start=tfstart(tf), end=tfend(tf), 
      series=seq(nfactors(x)),
      Title="Model factors", 
      lty = 1:5, lwd = 1, pch = NULL, col = 1:6, cex = NULL,
      xlab=NULL, ylab=factorNames(x), xlim = NULL, ylim = NULL, 
      graphs.per.page=5, par=NULL, reset.screen=TRUE) {

   tfplot(factors(x),...,  tf=tf, start=start, end=end, 
	series=series, Title=Title, 
        lty=lty, lwd=lwd, pch=pch, col=col, cex=cex,
        xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
	graphs.per.page=graphs.per.page, 
	par=par, reset.screen=reset.screen)
   }

#tfplot.TSFestModel <- function(x,...)  tfplot(factors(x), ...)

tfplot.TSFfactors <- function(x,..., tf=tfspan(x , ...),  
      start=tfstart(tf), end=tfend(tf), 
      series=seq(nfactors(x)),
      Title="Estimated factors (dashed) and true (solid)", 
      lty = c("dashed", "solid"), lwd = 1, pch = NULL, col = 1:6, cex = NULL,
      xlab=NULL, ylab=factorNames(x), xlim = NULL, ylim = NULL, 
      graphs.per.page=5, par=NULL, reset.screen=TRUE) {

   tfplot(unclass(x), attr(x, "true"),  ..., 
        tf=tf, start=start, end=end, 
	series=series, Title=Title, 
        lty=lty, lwd=lwd, pch=pch, col=col, cex=cex,
        xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
	graphs.per.page=graphs.per.page, 
	par=par, reset.screen=reset.screen)
   }


tfplot.TSFexplained <- function(x,..., tf=tfspan(x, ...), start=tfstart(tf), end=tfend(tf), 
      series=seq(nseries(x)),
      Title="Explained (dashed) and actual data (solid)", 
      lty = c("dashed", "solid"), lwd = 1, pch = NULL, col = 1:6, cex = NULL,
      xlab=NULL, 
      ylab=seriesNames(x), 
      xlim = NULL, ylim = NULL,
      graphs.per.page=5, par=NULL, reset.screen=TRUE) {
    tfplot( unclass(x), attr(x, "data"), 
                Title=Title,
                tf=tf, start=start, end=end, series=series,  
                lty=lty, lwd=lwd, pch=pch, col=col, cex=cex,
                xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
		graphs.per.page=graphs.per.page, 
		par=par, reset.screen=reset.screen)
  }


FAfitStats <- function(object, ...)UseMethod("FAfitStats")

FAfitStats.default <- function(object, diff.=TRUE, 
                N=(nrow(object) - diff.), 
		control=list(lower = 0.0001, opt=list(maxit=1000)), ...) {
    if (!is.matrix(object)) stop("FAfitStats.default expects a matrix object.")
    corDX  <- cor(if (diff.) diff(object) else object)
    #covDX  <- cov(if (diff.) diff(object) else object)
    nvar   <- ncol(object)
    maxfact <- floor(LedermannBound(nvar))
    OmegaTot <- matrix(NA,nvar, maxfact+1)
    
    # Fit statistics could be calculate with either corDX or covDX. Factanal
    # transforms a covariance matrix to a correlation matrix: the loadings and 
    # uniquenesses it returns are for the standardized solution (correlation 
    #  matrix) in whether it is given a cov or cor matrix. This means the
    # loadings an uniquenesses must be converted bac to the unstandardized
    # scale, if the fit staistics are to be calculated with the cov.
    # Commented code below does the calculation using covariances.

    # zero factors
    fitStats   <- FAmodelFitStats(NULL, NULL, diag(corDX), corDX, N)
    OmegaTot[,1] <- diag(corDX)
    nm <- list(names(fitStats), c(0, seq(maxfact), "saturated"))
    
    for (j in 1:maxfact) { # loop over number of factors
     	# estimate parameters
    	FA <- factanal(factors = j, covmat = corDX, n.obs = N,
    			scores = "none", rotation = "none", control=control)
    	#FA <- factanal(factors = j, covmat = covDX, n.obs = N,
    	#		scores = "none", rotation = "none", control=control)
     	OmegaTot[,j+1] <- FA$uniquenesses
        
        # compute fit statistics

 	fitStats <- cbind(fitStats, 
            FAmodelFitStats(FA$loadings, diag(1,j,j), FA$uniquenesses, corDX, N))
	#D <- diag(sqrt(diag(covDX)))
 	#fitStats <- cbind(fitStats, 
	#   FAmodelFitStats(D %*% FA$loadings, diag(1,j,j),
	#                               diag(covDX) * FA$uniquenesses, covDX, N))
        }

    Hey <- apply(control$lower == OmegaTot, 2, any)
    if(any(Hey)) warning("Heywood cases: ", 
                   paste((0:maxfact)[Hey], collapse=" "), " factor model(s)")

    # saturated model
    M <- nvar
    #k <- maxfact
    #nparc <- M * k + M - (k *(k-1))/2 # no. of param's corrected for rotation
    # above is for largest actual model, but saturated model may correspond
    # to a non-integer Ledermann bound
    nparc <- M *(M+1)/2 
    fitStats <- cbind(fitStats, c(
        0,0,NA,0,         # chisq=0, df=0, pval=na, delta=0
        NA,1,1,1,         # RMSEA=na, RNI=1, CFI=1, MCI=1
        1,1,0,            # GFI=1, AGFI=1, AIC=0
        (1 + log(N)) * nparc,  # CAIC
        log(N) * nparc,        # SIC
        (2 * nparc)/N,         # CAK
        2 * nparc/(N-M-1) ))    # CK


    dimnames(fitStats) <- nm
    
    # differences between consecutive models

    seqfitStats <- NULL
    for (j in 1:(maxfact+1))  seqfitStats <- cbind(seqfitStats, 
        c( fitStats["chisq",j] - fitStats["chisq",j+1],
	   fitStats["df",   j] - fitStats["df",   j+1],
	   pchisq(fitStats["chisq",j] - fitStats["chisq",j+1], 
	          fitStats["df",   j] - fitStats["df",   j+1], lower.tail=FALSE)))
    
    nm <- dimnames(fitStats)[[2]]
    dimnames(seqfitStats) <- list(c("chisq", "df", "pval"), 
                                  paste(nm[-length(nm)], "vs", nm[-1]))
        
    list(fitStats=fitStats, seqfitStats=seqfitStats)  #, OmegaTot= OmegaTot)
    }

FAfitStats.TSFmodel <- function(object, diff.=TRUE,
                             N=(nrow(object$data) - diff.), ...) {
    # This uses unstandardized B and Omega (and covDX rather than corDX).
    # This should be the same as standardized for MLE, but not for other
    #  estimation methods. Standardized may be better conditioned numerically,
    #  so might be perferred for ML, but probably does not seem to make much
    #  differene in simple tests. Might consider using both.
    
    X <- if(diff.) diff(object$data) else object$data
    FAmodelFitStats(loadings(object), object$Phi, diag(object$Omega),
                   cov(X), N)
    } 


FAmodelFitStats <- function(B, Phi, omega, S, N) {
  tr <- function(A) {sum(diag(A))} # local function, trace of a matrix

  # Fit statistics of FA model, based on standard likelihood.
  # No consistency checks are made.
  #
  # See, e.g., Wansbeek, T., & Meijer, E. (2000). Measurement error and latent
  #   variables in econometrics. Amsterdam: North-Holland. (W&M below)
  #
  # B     = loadings
  # Phi   = cov. matrix of factors
  # omega = vector of error variances
  # S     = sample covariance matrix, or correlation matrix if B and omega are
  #         standardized. (see the notes in FAfitStats.default)
  # N     = sample size (Many authors prefer sample size - 1.)
  # k     = number of factors (may be zero)

  # numbers of variables and factors

  M <- nrow(S)
  k <- if (is.null(B)) 0 else ncol(B)

  # Saturated model: all elements of cov. matrix are free parameters.

  const <- log(det(S)) + M

  # Null model: independence model (W&M, p. 305).

  Sigma0 <- diag(c(diag(S)))
  chisq0 <- N * (log(det(Sigma0)) + tr(solve(Sigma0) %*% S) - const)
  df0    <- 0.5 * (M^2 - M)
  delta0 <- max(0, chisq0 - df0)

  # Target model, i.e., the one from which B, Phi, and omega are obtained.

  # Model-implied covariance matrix
  if (k == 0)  Sigma <- diag(c(omega))
  else if (is.null(Phi)) Sigma <- B %*% t(B) + diag(c(omega))
  else                   Sigma <- B %*% Phi %*% t(B) + diag(c(omega))


  # Chi-square statistic, its degrees of freedom, and its p-value (W&M, p. 298).
  # Note: the df takes the rotational freedom into account (cf. W&M, p. 169).
  chisq <- N * (log(det(Sigma)) + tr(solve(Sigma) %*% S) - const)
  df    <- 0.5 * ( (M - k)^2  -  (M + k) )
  pval  <- pchisq(chisq, df, lower.tail=FALSE)

  # Estimate of noncentrality parameter (W&M, p. 307).
  delta <- max(0, chisq - df)

  # Comparative fit index (W&M, p. 307).
  if(chisq0 <= df0)       CFI <- 0 # Null model fits very well: extremely unlikely.
  else if(chisq <= df)    CFI <- 1 # Target model fits very well.
  else if(delta0 < delta) CFI <- 0 # Null model fits better than target model: also extremely unlikely.
  else        CFI <- 1 - delta/delta0  # The most common situation: null model fits very badly,
                             # target model fits better, but not perfectly.

  # Root mean square error of approximation (W&M, p. 309).
  RMSEA <- if (df > 0)  sqrt(delta/(N * df)) else  Inf


  # Dozens of other fit indexes possible. See output of LISREK, EQS,
  # Amos, Mplus, and Mx, and the paper by Ogasawara (SEM, ca. 2001).
  # Most are very bad. Here are some possibilities, from W&M (chap. 10),
  # Hu & Bentler (1995), and Bollen (1989, pp. 256-281):

   # NFI  <- 1 - chisq/chisq0 # Normed fit index
   # TLI  <- (chisq0/df0 - chisq/df)/(chisq0/df0 - 1) # Tucker-Lewis Index,
  	       # also called Nonnormed fit index
   # BL86 <- 1 - (chisq/df)/(chisq0/df0) # Bollen (1986)
   # BL89 <- (chisq0 - chisq)/(chisq0 - df) # Bollen (1989)
   RNI  <- 1 - delta/delta0 # Relative noncentrality index
   MCI  <- exp(-0.5 * (chisq - df)/N) # McDonald's centrality index
  
  # # Some indexes by Joreskog & Sorbom (1981, 1986):
  
   T1	<- solve(Sigma) %*% S  # Temporary matrix
   T2	<- T1 - diag(1,M,M)    # Temporary matrix
   GFI  <- 1 - tr(T2 %*% T2)/tr(T1 %*% T1) # Goodness of fit index
   AGFI <- 1 - ((M * (M+1))/(2 * df)) * (1 - GFI) # Adjusted GFI
   #RMR <- sqrt(((sum( (   c(S - Sigma))^2 ) +
   #		 sum( (diag(S - Sigma))^2 ))/(M *(M+1))) #Root mean-square residual

  # # Hoelter's (1983) Critical N; note that the significance level
  # # alpha must be given. E.g.:
  # # alpha <- 0.05
  # # Apparently, this is the definition of Hoelter:
  # CN <- (qnorm(alpha/2, lower.tail=FALSE) + sqrt(2*df-1))^2/(2*chisq/N) + 1
  # # but I've also seen the following, which makes more sense:
  # # CN <- qchisq(alpha, df, lower.tail=FALSE)/(chisq/N) + 1
  # # although both are not very useful.

  # # Some information criteria; accounting for rotational freedom
   nparc <- M * k + M - (k *(k-1))/2 # no. of param's corrected for rotation
   AIC   <- chisq - 2 * df # or chisq + 2 * nparc # Akaike's info. crit.
   CAIC  <- chisq + (1 + log(N)) * nparc # Consistent AIC
   SIC   <- chisq + log(N) * nparc # Schwarz's Bayesian info. crit.
   CAK   <- (chisq + 2 * nparc)/N  # Cudeck & Browne's rescaled AIC
   CK	 <- chisq/N + 2 * nparc/(N-M-1) # Cudeck & Browne's cross-val. index

  r <- c(chisq, df, pval, delta, RMSEA, RNI, CFI,
         MCI, GFI, AGFI, AIC, CAIC, SIC, CAK, CK) #NFI, TLI, BL86, BL89

  names(r) <- c("chisq", "df", "pval", "delta", "RMSEA", "RNI", "CFI",
         "MCI", "GFI", "AGFI", "AIC", "CAIC", "SIC", "CAK", "CK")
  r
}

summary.TSFmodel <- function(object, ...)
 {fitStats <- FAfitStats(object)
  est  <- TSFmodel(object)
   
  barx     <- colMeans(object$data)
  barx.est <- colMeans(explained(est))
  hatk     <- colMeans(factors(est))

  barDx     <-  colMeans(diff(object$data ))
  barDx.est <-  colMeans(diff(explained(est)))
  hatDk     <-  colMeans(diff(factors(est)))

  true <- attr(object$data, "TSFmodel")
  if (is.null(true)) 
      B.true <- hatk.true <- hatDk.true <- NULL 
  else  {
      B.true     <- true$loadings
      hatk.true  <- colMeans(true$f)
      hatDk.true <- colMeans(diff(factors(true)))
      }
  
  classed(list(
      N=Tobs(factors(object)), S=start(factors(object)), E=end(factors(object)),
      Snames=seriesNames(object),Fnames=factorNames(est),
      fitStats=fitStats, B.estimate=est$loadings,   B.true=B.true,
      #stdB.estimate=standardizedLoadings(object), 
      DstdB.estimate=DstandardizedLoadings(object),
      barDx=barDx, barDx.est=barDx.est,   barx=barx,  barx.est=barx.est,
      hatDk=hatDk, hatDk.true=hatDk.true, hatk=hatk,  hatk.true=hatk.true),
	  "summary.TSFmodel")
  }


   #positive.data=object$positive.data
   #cat("positive.data ", (if(x$positive.data) "is" else "is not"), " specified.\n")

print.summary.TSFmodel <- function (x, ...)
  {cat("factors have ", x$N, " observations from:", x$S, " to ", x$E, "\n")
   cat("     Estimated loadings:\n"); print(x$loadings.estimate)
   cat("\n     Standardized (using differenced data covariance):\n")
   print(x$DstdB.estimate)
   #cat("\n     Standardized (using undifferenced data covariance):\n")
   #print(x$stdB.estimate) 
   if (!is.null(x$loadings.true))
     {cat("\n   true loadings:\n"); print(x$loadings.true)
      cat("\n   loadings estimation error:\n"); print(x$loadings.estimate - x$loadings.true)
     }


   z <- rbind(x$barx.est,x$barx,  x$barx.est - x$barx)
   dimnames(z) <- list(c("explained","actual","error"), x$Snames)
   cat("\n                 Mean of data:\n"); print(z)

   z <- rbind(x$barDx.est,x$barDx,  x$barDx.est - x$barDx)
   dimnames(z) <- list(c("explained","actual","error"), x$Snames)
   cat("\n		  Mean of differenced data:\n"); print(z) 

   if (!is.null(x$hatk.true))
     {z <- rbind(x$hatk, x$hatk.true,  x$hatk - x$hatk.true)
      dimnames(z) <- list(c("estimated","true","error"), x$Fnames)
     }
   else
     {z <- x$hatk
      names(z) <- x$Fnames
     }
   cat("\n     Mean of factors:\n"); print(z)
   
   if (!is.null(x$hatDk.true))
     {z <- rbind(x$hatDk, x$hatDk.true,  x$hatDk - x$hatDk.true)
      dimnames(z) <- list(c("estimated","true","error"), x$Fnames)
     }
   else 
     {z <- x$hatDk
      names(z) <- x$Fnames
     }
   cat("\n    Mean of differenced factors:\n"); print(z)

   cat("\n   Fit statistics:\n")
   print(x$fitStats)
   invisible(x)
  }



distribution.factorsEstEval <- function (obj, ..., bandwidth = "nrd0",
        cumulate=TRUE, graphs.per.page = 5, Title=NULL)
  {# if cumulate is true then a distribution is plotted, otherwise,
   # a time series graph of the true and one 1 sd bands
    truth <- obj$truth
    r <- array(NA, c(length(obj$result), dim(truth)))
    otherobj <- list(...)
    obr <- list()
    for (ob in otherobj)
      {if (! testEqual(truth, ob$truth))
                    warning("object true values do not correspond.")
       rx <- r
       for (i in 1:length(ob$result)) rx[i,,] <- ob$result[[i]] - truth
       obr <- append(obr, list(rx))
      }
    for (i in 1:length(obj$result)) r[i,,] <- obj$result[[i]] - truth
    xlab <- "factor "
    old.par <- par(par)
    on.exit(par(old.par)) 
    if (cumulate){
       par(mfcol = c(min(graphs.per.page, ncol(truth)), 1), no.readonly = TRUE)
       for (i in 1:ncol(truth))
         {rd <- density(c(r[,,i]), bw = bandwidth)
          rdy <- rd$y
	  for (rx in obr)
	         rdy <- cbind(rdy, density(c(rx[,,i]), bw = bandwidth)$y)
	  matplot(rd$x, rdy, type = "l",
	     ylab = "density", xlab = paste(xlab, i), main = "")
          if(!is.null(Title) && (i==1) && (is.null(options()$PlotTitles)
                || options()$PlotTitles)) title(main = Title)
	 }
       }
    else
      {rd <- apply(r,c(2,3), FUN="var")^0.5
       tfplot(truth, truth + rd, truth - rd,
                   Title=Title, graphs.per.page = graphs.per.page)
      }
    invisible()
}

checkResiduals.TSFmodel <- function (obj, data=obj$data, diff.=TRUE, ...) {
	res <- if (diff.) diff(explained(obj)) - diff(data)
	           else explained(obj) - data
	seriesNames(res) <- seriesNames(data)
	cat("residual covariance matrix\n")
	cv <- cov(res)
	print(cv)
	cat("\nsum of trace cov: ", sum(diag(cv)), "\n")
	cat("sum of abs (off-diag of cov): ", sum(abs(cv - diag(cv))), "\n")
	checkResiduals(res, ...)
}


summaryStats <- function(object, ...) UseMethod("summaryStats")

summaryStats.TSFmodelEstEval <- function(object, ...) {

  N <- length(object$result)
  if (N <2 ) stop("This requires more than one replication.")

  meanhatf <- sdhatf <- meanhatDf <- sdhatDf <- meanhatPCf <- 
              sdhatPCf <- meanhatB <- sdhatB <- 
	      estConverged <- rotationConverged <- 0

  for (m in object$result) { 
      meanhatf  <- meanhatf   + factors(m)
      sdhatf	<- sdhatf     + factors(m)^2
      meanhatDf <- meanhatDf  + diff(factors(m))
      sdhatDf   <- sdhatDf    + diff(factors(m))^2
      meanhatPCf<- meanhatPCf + percentChange(factors(m))
      sdhatPCf  <- sdhatPCf   + percentChange(factors(m))^2
      meanhatB  <- meanhatB   + TSFmodel(m)$loadings
      sdhatB	<- sdhatB     + TSFmodel(m)$loadings^2
      if (!is.null(TSFmodel(m)$dots$estConverged) && !TSFmodel(m)$dots$estConverged)
          estConverged <- estConverged + 1
      if (!is.null(TSFmodel(m)$dots$rotationConverged) && !TSFmodel(m)$dots$rotationConverged)
          rotationConverged <- rotationConverged + 1
      }
 
  true <- factors(object$truth)
  tf <- tframe(true)
  dtf <- tframe(diff(true))
  
  meanhatf   <- meanhatf   /N
  meanhatDf  <- meanhatDf  /N
  meanhatPCf <- meanhatPCf /N
  meanhatB   <- meanhatB   /N
 
  sdhatf   <- sqrt(sdhatf   /N - meanhatf^2)
  sdhatDf  <- sqrt(sdhatDf  /N - meanhatDf^2)
  sdhatPCf <- sqrt(sdhatPCf /N - meanhatPCf^2)
  sdhatB   <- sqrt(sdhatB   /N - meanhatB^2)
  
  list(true=true, Btrue=TSFmodel(object$truth)$loadings,
  	meanhatf   =  tframed(meanhatf,   tf), 
  	meanhatDf  =  tframed(meanhatDf, dtf), 
  	meanhatPCf =  tframed(meanhatPCf,dtf), 
  	meanhatB   =          meanhatB,   
  	sdhatf     =  tframed(sdhatf,     tf),	
  	sdhatDf    =  tframed(sdhatDf,   dtf),  
  	sdhatPCf   =  tframed(sdhatPCf,  dtf), 
  	sdhatB     =          sdhatB,
        estConverged = estConverged,
        rotationConverged = rotationConverged)
  }


summary.TSFmodelEstEval <- function(object, ...) {
  sm <- summaryStats(object, ...)
  classed(list(
      meanhatf.error  = colMeans(sm$meanhatf  - sm$true),
      meanSDhatf      = colMeans(sm$sdhatf),
      meanhatDf.error = colMeans(sm$meanhatDf - diff(sm$true)),
      meanSDhatDf     = colMeans(sm$sdhatDf),
      meanhatPCf.error= colMeans(sm$meanhatPCf -
                               percentChange(sm$true)),
      meanSDhatPCf    = colMeans(sm$sdhatPCf),
      meanhatB.error  = sm$meanhatB - sm$Btrue, 
      SDhatB          = sm$sdhatB,
      estConverged = sm$estConverged,
      rotationConverged = sm$rotationConverged), "summary.TSFmodelEstEval")
  }

print.summary.TSFmodelEstEval <- function(x, digits = options()$digits, ...) {
  cat("    mean hat{f} error\n") ; print(x$meanhatf.error, digits=digits)
  cat("\n    mean  SD hat{f}\n")   ; print(x$meanSDhatf, digits=digits)
  cat("\n    mean diff hat{f} error\n");  print(x$meanhatDf.error,  digits=digits)
  cat("\n    mean %change hat{f} error\n");  print(x$meanhatPCf.error,  digits=digits)
  cat("\n    mean hat{B} error") ; print(x$meanhatB.error, digits=digits) 
  cat("\n      SD hat{B}")       ; print(x$SDhatB, digits=digits) 
  cat("\n    Estimates NOT converged: ", x$estConverged) 
  cat("\n    Rotations NOT converged: ", x$rotationConverged) 
  cat("\n") 
  invisible(x)
  }


tfplot.TSFmodelEstEval <- function(x, ...,  tf=NULL, start=tfstart(tf), end=tfend(tf), 
		 series=seq(nseries(factors(x))),
		 Title="Monte Carlo Results", 
		 lty = c("solid", "dotdash", "dashed",  "dashed"), lwd = 1, pch = NULL, 
		 col = c("black", "red", "red", "red"), cex = NULL,
		 xlab=NULL, 
		 ylab=seriesNames(factors(x$truth)), 
		 xlim = NULL, ylim = NULL,
		 graphs.per.page=5, par=NULL, reset.screen=TRUE,
		 diff.=FALSE,  percentChange.=FALSE,
                 PCcentered.=FALSE, summary.=TRUE) {

   #if summary. is FALSE then all of the factors are plotted,
   # otherwise the mean and 1 SD bounds are plotted as follows:
   #if diff. is TRUE then the differenced factors are plotted
   #if percentChange. is TRUE then the PC factors are plotted
   #if PCcentered. is TRUE then the PC factors less means are plotted
   # otherwise the undifferenced factors are plotted

   true <- factors(x$truth)
   
   if(!summary.) {
     if(is.null(tf)) tf <- tframe(true)
     tfplot(factors(x), tf = tf, start=start, end=end, series=series,
        truth = true,
        Title = Title,
      ylab = seriesNames(true), remove.mean = FALSE, graphs.per.page = 5,
       par = par, reset.screen = TRUE, ...)
      } 
   else {
      sm <- summaryStats(x)

      if(diff.) { # factor difference
          df <- diff(true)
	  if(is.null(tf)) tf <- tframe(df)
    	  tfplot(df,
    		 sm$meanhatDf, 
    		 sm$meanhatDf  + 1.96 * sm$sdhatDf,  
    		 sm$meanhatDf  - 1.96 * sm$sdhatDf,
		Title=Title,
                tf=tf, start=start, end=end, series=series,  
                lty=lty, lwd=lwd, pch=pch, col=col, cex=cex,
                xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
		graphs.per.page=graphs.per.page, 
		par=par, reset.screen=reset.screen)
    	  }
      else if(percentChange.){ # factor growth rates
          pc <- percentChange(true)
	  if(is.null(tf)) tf <- tframe(pc)
    	  tfplot(pc,
    		 sm$meanhatPCf, 
    		 sm$meanhatPCf  + 1.96 * sm$sdhatPCf,  
    		 sm$meanhatPCf  - 1.96 * sm$sdhatPCf,
		Title=Title,
                tf=tf, start=start, end=end, series=series,  
                lty=lty, lwd=lwd, pch=pch, col=col, cex=cex,
                xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
		graphs.per.page=graphs.per.page, 
		par=par, reset.screen=reset.screen)
    	  }
      else if(PCcentered.){ # factor growth rates: bias (clearer picture?)
    	  growth <- percentChange(true)
          if(is.null(tf)) tf <- tframe(growth)
	  tfplot(sm$meanhatPCf - growth, 
    		sm$meanhatPCf  + 1.96 * sm$sdhatPCf - growth, 
    		sm$meanhatPCf  - 1.96 * sm$sdhatPCf - growth,
		Title=Title,
                tf=tf, start=start, end=end, series=series,  
                lty=lty, lwd=lwd, pch=pch, col=col, cex=cex,
                xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
		graphs.per.page=graphs.per.page, 
		par=par, reset.screen=reset.screen)
          }
      else { # factors
          if(is.null(tf)) tf <- tframe(true)
    	  tfplot(true,
    		 sm$meanhatf, 
    		 sm$meanhatf  + 1.96 * sm$sdhatf, 
    		 sm$meanhatf  - 1.96 * sm$sdhatf,
		Title=Title,
                tf=tf, start=start, end=end, series=series,  
                lty=lty, lwd=lwd, pch=pch, col=col, cex=cex,
                xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
		graphs.per.page=graphs.per.page, 
		par=par, reset.screen=reset.screen)
    	  }
      }
   invisible(x)
   }


predict.TSFmodel <- function(object, data = object$data, factorNames.=factorNames(object), ...){
      # prediction of factors with data
      if (is.null(data)) stop("data must be supplied.")
      tframed(data %*% t(object$LB), tframe(data), names=factorNames.) #hatf
      }


#explained.TSFestModel <- function (object, ...)
# {r <- explained(TSFmodel(object), names=seriesNames(object$data), ...) 
#  attr(r, "data") <- object$data
#  r
# }

#explained.TSFmodel <- function (object, f=factors(object),
#                  names=seriesNames(object), ...) {
#  # portion of data explained by factors
#  classed(tframed(t(loadings(object) %*% t(f)), tf=tframe(f), names=names),
#     "TSFexplained") 
#  }

explained.TSFmodel <- function (object, f=factors(object),
                  names=seriesNames(object), ...) {
  # portion of data explained by factors
  tframed(t(loadings(object) %*% t(f)), tf=tframe(f), names=names)
  }


#######################################

estTSF.ML <- function(y, p, diff.=TRUE, 
                      rotation=if(p==1) "none" else "quartimin", 
		      rotationArgs=NULL, 
		      normalize=TRUE, eps=1e-5, maxit=1000, Tmat=diag(p),
		      BpermuteTarget=NULL,
                      factorNames=paste("Factor", seq(p))) {
   # it would be better to have GPFargs=list(Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit)

      # Estimate parameters using standard (quasi) ML factor analysis
      # (on the correlation matrix and then scaled back).
      # factanal always uses the cor matrix, so standardizing does not affect 
      # the solution. Both standardized and not can be calculated after.
      # With non ML methods this solutions may differ (and working with cov  
      # rather than cor is probabably better.
   
   estTSFmodel(y, p, diff.=diff., 
        est="factanal", 
	estArgs=list(scores="none", control=list(opt=list(maxit=10000))),
        rotation=rotation, 
	rotationArgs=rotationArgs, 
	GPFargs=list(Tmat=diag(p), normalize=normalize, eps=eps, maxit=maxit),
	BpermuteTarget=BpermuteTarget,
        factorNames=factorNames)
   }

estTSFmodel <- function(y, p, diff.=TRUE, 
                est="factanal", 
		estArgs=list(scores="none", control=list(opt=list(maxit=10000))),
                rotation=if(p==1) "none" else "quartimin", 
		rotationArgs=NULL, 
		GPFargs=list(Tmat=diag(p),normalize=TRUE, eps=1e-5, maxit=1000), 
		BpermuteTarget=NULL,
                factorNames=paste("Factor", seq(p))) {
      if (p < 1) stop("number of factors must be a positive integer.")
      indicatorNames <- seriesNames(y)
      zz <- if (diff.) diff(y) else y
      zz <- sweep(zz,2,colMeans(zz), "-")
      Sigma  <- crossprod(zz)/(Tobs(zz) - 1)
      
      z <- estFAmodel(Sigma, p, n.obs=(Tobs(y) - diff.),
                est="factanal", 
		estArgs=estArgs,
                rotation=rotation, rotationArgs=rotationArgs, 
		GPFargs=GPFargs,
		BpermuteTarget=BpermuteTarget,
                factorNames=factorNames,
                indicatorNames=indicatorNames) 

      model <- TSFmodel(z,
                        f=tframed((y %*% diag(1/sqrt(diag(Sigma)))) %*% t(z$LB.std), tframe(y), names=factorNames), #hatf
			positive.data=all(0<y)
		       )
			 
       model$data <- y
       classed(model, c("TSFmodel", "fFAmodel", "FAmodel"))
      }

#######################################
#######################################
