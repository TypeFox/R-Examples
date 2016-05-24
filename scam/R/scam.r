

#############################################################
## the wrapper overall Function to fit scam...             ##
#############################################################

scam <- function(formula, family=gaussian(), data=list(), gamma=1, sp=NULL, 
        weights=NULL, offset=NULL, optimizer="bfgs",
            optim.method=c("Nelder-Mead","fd"),
           scale=0, devtol=1e-8, steptol=1e-8, check.analytical=FALSE, del=1e-4,
              start=NULL, etastart, mustart,keepData=FALSE)
{  ## scale - scale parameter of the exponential deistribution as in gam(mgcv)
   ## devtol - a scalar giving the tolerance at which the relative penalized deviance is considered to be close enougth to 0 to terminate the algorithm
   ## steptol - a scalar giving the tolerance at which the scaled distance between two successive iterates is considered close enough to zero to terminate the algorithm
   ## optimizer - numerical optimization method to use to optimize the smoothing      parameter estimation criterion: "bfgs", "optim", "nlm", "nlm.fd"
   ## optim.method - if optimizer=="optim" then the first argument of optim.method     specifies the method, and the second can be either "fd" for finite-difference approximation of the gradient or "grad" - to use analytical gradient of gcv/ubre
   ## check.analytical - logical whether the analytical gradient of GCV/UBRE should be checked
   ## del - increment for finite differences when checking analytical gradients
   ### --------------------------------------------------------

   ## Setting from mgcv(gam).......
 #  require(mgcv)
    
   G <- gam(formula, family,data, fit=FALSE) 
   n.terms <- length(G$smooth)  ## number of smooth terms in the model
   n <- nrow(G$X)
   intercept <- G$intercept ## TRUE or FALSE
   ## now need to set 'offset' as the above G wouldn't take in 'offset' that is outside of formula..
   gp <- interpret.gam(formula) # interpret the formula 
   cl <- match.call() # call needed in gam object for update to work
   mf <- match.call(expand.dots=FALSE)
   mf$formula <- gp$fake.formula 
   mf$family <- mf$control<-mf$scale<-mf$knots<-mf$sp<-mf$min.sp<-mf$H<-mf$select <-
                 mf$gamma<-mf$method<-mf$fit<-mf$paraPen<-mf$G<-mf$optimizer <- mf$optim.method <- mf$in.out <- mf$...<-NULL
   mf[[1]] <- as.name("model.frame")
   pmf <- mf
   mf <- eval(mf, parent.frame()) 
   G$offset <- as.vector(model.offset(mf)) 
   if (is.null(G$offset)) 
        G$offset <- rep.int(0, n)
   ## done offset 
   
   if (is.null(weights)) 
         weights <- rep.int(1,n)
   fam.name <- G$family[1]
   if (scale == 0) 
      {  if (fam.name == "binomial" || fam.name == "poisson") 
               sig2 <- 1
         else sig2 <- -1
      }
   else { sig2 <- scale }
   if (sig2 > 0) scale.known <- TRUE else scale.known <- FALSE
   ## get penalty matrices and 
   ## vector of identifications for exponentiated model coefficients...
   Q <- penalty_pident(G)
   ## checking sp...
   if (!is.null(sp)) {
         neg <- FALSE
         if (length(sp)!= length(G$off)) {
              warning("Supplied smoothing parameter vector is too short - ignored.")
              sp <- NULL       
         } else if (sum(is.na(sp))) {
               warning("NA's in supplied smoothing parameter vector - ignoring.")
               sp <- NULL
           } else {
                good <- sp < 0
                if (sum(good) > 0) { ## cheking negative values..
                   warning("Supplied smoothing parameter vector has negative values - ignored.")
                   neg <- TRUE
                }
             }                
         if (neg) sp <- NULL
     }
   ## Create new environments with `start' initially empty...
 #  ee <- new.env() ## parent = .GlobalEnv
   env <- new.env() 
   assign("start",rep(0,0),envir=env)
   assign("dbeta.start",rep(0,0),envir=env)
   assign("sp.last",rep(0,0),envir=env)
   
   q.f <- rep(0,n.terms)
   for (i in 1:n.terms) 
            { q.f[i] <- ncol(G$smooth[[i]]$S[[1]]) + 1 }
   G$S <- Q$S
   G$q.f <- q.f
   G$q0 <- G$off[1]-1  ## number of the parameters of the strictly parametric model
   G$p.ident <- Q$p.ident  ## vector of 0's & 1's for the model parameters identification: 
   G$n.terms <- n.terms   ## number of the smooth terms in the SCAM
  # G$intercept <- intercept
   G$weights <- weights
   G$sig2 <- sig2
   G$scale.known <- scale.known

   if (!keepData) rm(data) ## save space

   object <- list() 
   if (is.null(sp)) { 
       ## get initial estimates of the smoothing parameter...
         start <- etastart <- mustart <- NULL
         y <- G$y; family <- G$family
         nobs <- NROW(y)
         eval(family$initialize)
         G$y <- y  ## needed to set factor response values as numeric
         def.sp <- initial.sp.scam (G,Q,q.f=q.f,n.terms=n.terms,family=family,
              intercept=intercept,offset=G$offset, env=env,
              weights=weights, devtol=1e-4, steptol=1e-4)
         rho <- log(def.sp) ## get initial log(sp) ...
         ## minimize GCV/UBRE by optimizer....
         ptm <- proc.time()
         re <- estimate.scam(G=G,optimizer=optimizer,optim.method=optim.method,
               rho=rho, gamma=gamma, env=env,
              check.analytical=check.analytical, del=del, devtol=devtol, steptol=steptol)
         CPU.time <- proc.time()-ptm
         best <- re
         object$gcv.ubre <- re$gcv.ubre
         object$dgcv.ubre <- re$dgcv.ubre
         object$aic <- re$aic
         best$p.ident <- Q$p.ident
         best$S <- Q$S
         object$optimizer <- optimizer
         object$edf1 <- re$edf1
         object$termcode <- re$termcode
         if (optimizer == "bfgs")
            {   object$check.grad <- re$check.grad
                object$dgcv.ubre.check <- re$dgcv.ubre.check
            }
   } else {   ## no GCV minimization if sp is given...
            best <- scam.fit(G=G, sp=sp,gamma=gamma,devtol=devtol, steptol=steptol, env=env)
            object$aic <- best$aic
            object$optimizer <- "NA"           
      }
   ## post-fitting values...
   best$n.smooth <- object$n.smooth <- n.terms
   best$formula <- object$formula <- formula
   best$family <- object$family <- G$family
   best$smooth <- object$smooth <- G$smooth
   best$model <- object$model <- G$mf

   object$R <- best$R
   if (is.null(object$R)){
         rr <- scam.fit(G=G, sp=best$sp,gamma=gamma,devtol=devtol, steptol=steptol, env=env)
         object$R <- rr$R } ## not sure if it's needed?
  
   object$df.residual <- nrow(best$X) - sum(best$edf)

   object$sp <- best$sp
   names(object$sp) <- names(G$sp)
   if (sum(is.na(names(object$sp)))!=0){  ## create names for sp if NA's from G
      for (i in 1:n.terms) names(object$sp)[i] <- object$smooth[[i]]$label
   }
   object$deviance <- best$deviance
   object$residuals <- best$residuals
#   object$X <- best$X

   object$conv <- best$conv # whether or not the inner full Newton method converged
   post <- scam.fit.post(y=G$y,X=G$X,object=best,sig2=sig2,offset = G$offset,
                   intercept=G$intercept, weights=weights,scale.known=scale.known) 

   object$edf <- post$edf
   object$edf1 <- post$edf1
   object$trA <- post$trA
   names(object$edf) <- G$term.names
   names(object$edf1) <- G$term.names

   object$null.deviance <- post$nulldev
   object$var.summary <- G$var.summary 
   object$cmX <- G$cmX ## column means of model matrix --- useful for CIs
   object$model<-G$mf # store the model frame
   
   object$full.sp <- G$full.sp  ## incorrect !!!
   if (!is.null(object$full.sp))   names(object$full.sp) <- names(G$full.sp)

   object$na.action <- attr(G$mf,"na.action") # how to deal with NA's
   object$df.null <- post$df.null
   object$Ve <- post$Ve
   object$Vp <- post$Vb
   object$Ve.t <- post$Ve.t
   object$Vp.t <- post$Vb.t
   object$sig2 <- post$sig2
   object$coefficients <- best$beta
   object$coefficients.t <- best$beta.t
   object$beta <- best$beta
   object$beta.t <- best$beta.t
   object$pterms <- G$pterms
   object$terms <- G$terms
   object$assign <- G$assign
   object$nsdf <- G$nsdf
   object$y <- G$y
 #  object$data <- G$mf
   if (keepData) object$data <- data 
   object$offset <- G$offset
 #  object$scale.known <- scale.known # to be passed in the summary function
   object$scale.estimated <- !scale.known # to be passed in the summary function
   object$prior.weights <-weights # prior weights on observations
   object$weights <- best$w  # final weights used in full Newton iteration
   object$fitted.values <- best$mu
   object$linear.predictors <- best$eta
 #  cl<-match.call() # call needed in gam object for update to work
   object$call <- cl 
   object$p.ident <- Q$p.ident
   object$intercept <- G$intercept
   object$min.edf <- G$min.edf ## Minimum possible degrees of freedom for whole model
   object$gamma <- gamma
   object$iter <- best$iter  # number of iterations of the Full Newton
   if (is.null(sp)) 
         object$CPU.time <- CPU.time
   else 
         object$CPU.time <- NULL

   ## get the optimizer info (smoothing parameter selection).....
  if (is.null(sp))
     {   if (optimizer == "bfgs")
            {   ## get the bfgs info in case of sp selection... 
                object$bfgs.info <- list()
                object$bfgs.info$conv <- re$conv.bfgs  
                object$bfgs.info$iter <- re$iterations 
                object$bfgs.info$grad <- re$dgcv.ubre
            }
         else if (optimizer == "nlm.fd" || optimizer == "nlm")
                 {   object$nlm.info <- list()
                     object$nlm.info$conv <- re$conv 
                     object$nlm.info$iter <- re$iterations 
                     object$nlm.info$grad <- re$dgcv.ubre
                 }
         else if (optimizer=="optim")
                 {   object$optim.info <- list()
                     object$optim.info$conv <- re$conv 
                     object$optim.info$iter <- re$iterations 
                     object$optim.method <- re$optim.method  
                 }
      }
   if (scale.known)
         object$method <- "UBRE"
   else  object$method <- "GCV" 
   if (G$nsdf > 0) 
         term.names <- colnames(G$X)[1:G$nsdf]
   else term.names <- array("", 0)
   if (n.terms) 
        for (i in 1:n.terms) 
            {   k <- 1  
                for (j in G$smooth[[i]]$first.para:G$smooth[[i]]$last.para) 
                   {   term.names[j] <- paste(G$smooth[[i]]$label, ".", 
                                         as.character(k), sep = "")
                       k <- k + 1
                   }
            }
   names(object$coefficients) <- term.names
   names(object$coefficients.t) <- term.names
   ynames <- if (is.matrix(G$y)) 
                 rownames(G$y)
             else names(G$y)
   names(object$residuals) <- ynames
   class(object) <- "scam"   
   object
}  ## end scam



#################################################################
## function to get initial estimates of smoothing parameters...##
#################################################################

initial.sp.scam <- function(G,Q,q.f,n.terms,family,intercept,offset, env= env,
                      weights,devtol=1e-4,steptol=1e-4) 
{  ## function to get initial estimates of smoothing parameters
   ## step 1: set sp=rep(0.5,p) and estimate hessian...
   b <- scam.fit(G=G,sp=rep(0.5,length(G$off)), devtol, steptol, env=env)
   H <- crossprod(b$wX1) - b$E
   ## step 2:...
   n.p <- length(Q$S) ## number of penalty matrices
   def.sp <- array(0,n.p) ## initialize the initial sp values
   j <- 1
   for (i in 1:n.terms)
       {   for (kk in 1:length(G$smooth[[i]]$S))
               {   start <- G$off[j]
                   finish <- start + ncol(G$smooth[[i]]$S[[kk]])-1
                   # matrix norm of the Hessian elements penalized by S[[kk]]...
                   Hi.norm <- sum(H[start:finish,start:finish]*H[start:finish,start:finish]) 
                   Si.norm <- sum(G$smooth[[i]]$S[[kk]]*G$smooth[[i]]$S[[kk]])
                   def.sp[j] <- (Hi.norm/Si.norm)^0.5
                   j <- j+1
               }
       }
   ## Create again new environments with `start' initially empty...
   env <- new.env()
   assign("start",rep(0,0),envir=env)
   assign("dbeta.start",rep(0,0),envir=env)
   assign("sp.last",rep(0,0),envir=env)
   def.sp
}


#########################################################
## function to get list of penalty matrices and        ## 
## vector of parameter identifications .....           ##
#########################################################


penalty_pident <- function(object)
{  ## function to get the list of penalties and vector of model parameters 
   ## identifications from the gam() setting...
   n.terms <- length(object$smooth)  # number of terms in the model
   q <- ncol(object$X)          # total number of parameters
   cons.terms <- rep(0,n.terms) # define whether each term is constrained or not
   for (i in 1:n.terms)
       {   if (!is.null(object$smooth[[i]]$p.ident))
           cons.terms[i] <- 1  
       }
   p.ident <- rep(0,q) # initialize vector of parameter identifications
                      # with `1' - for a parameter to be exponentiated, `0' - otehrwise
   off.terms <- rep(0,n.terms) # starting points for each term
   off <- object$off
   if (n.terms ==length(off))
          off.terms <- off
   else 
      {   off.terms[1] <- off[1]
          k <- 1
          l <- 1
          while (l<length(off))
              {   if (off[l]!=off[l+1])
                     {   off.terms[k+1] <- off[l+1] 
                         k <- k+1; l <- l+1 
                     } 
                  else l <- l+1
              }
     
      }
   for (i in 1:n.terms)
       {   if (cons.terms[i]==1) 
              p.ident[off.terms[i]:(off.terms[i]+ncol(object$smooth[[i]]$S[[1]])-1)] <- 
                 object$smooth[[i]]$p.ident
       }
   ## getting the list of penalty matrices in terms of the full model vector of coefficients...
   S <- list()
   j <- 1
   for(i in 1:n.terms)
        { for (kk in 1:length(object$smooth[[i]]$S))
              {    S[[j]] <- matrix(0,q,q) # initialize penalty matrix
                   S[[j]][off.terms[i]:(off.terms[i]+ncol(object$smooth[[i]]$S[[kk]])-1),
                       off.terms[i]:(off.terms[i]+ncol(object$smooth[[i]]$S[[kk]])-1)] <- object$smooth[[i]]$S[[kk]]
                   j <- j+1       
              }
        }
   object$S <- S 
   object$p.ident <- p.ident
   object
}





#############################################################
## Function to fit SCAM based on Full Newton method        ##     
#############################################################

scam.fit <- function(G,sp, maxit=200,devtol=1e-8, steptol=1e-8,
                gamma=1, start=NULL, etastart=NULL, mustart=NULL, env=env)
   ## G - list of items from gam(...,fit=FALSE) needed to fit a scam
   ## sp- vector of smoothing parameters
   ## maxit - a positive scalar which gives the maximum number of iterations for Newton's method
   ## devtol - a scalar giving the tolerance at which the relative penalized deviance is considered to be close enougth to 0 to terminate the algorithm
   ## steptol - a scalar giving the tolerance at which the scaled distance between two successive iterates is considered close enough to zero to terminate the algorithm 
{ y <- G$y;  X <- G$X;  S <- G$S
  attr(X,"dimnames") <- NULL
  q0 <- G$q0; q.f <- G$q.f
  p.ident <- G$p.ident; n.terms <- G$n.terms
  family <- G$family; intercept <- G$intercept; offset <- G$offset;
  weights <- G$weights;  
  n <- nobs <- NROW(y)
  q <- ncol(X)
  dg <- fix.family.link(family)
  dv <- fix.family.var(family)
  nvars <- ncol(X)
  EMPTY <- nvars == 0
  variance <- family$variance
  linkinv <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv)) 
          stop("'family' argument seems not to be a valid family object")
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  mu.eta <- family$mu.eta
  
 
  if (!is.function(variance) || !is.function(linkinv)) 
        stop("illegal `family' argument")
  valideta <- family$valideta
  if (is.null(valideta)) 
        valideta <- function(eta) TRUE
  validmu <- family$validmu
  if (is.null(validmu)) 
        validmu <- function(mu) TRUE
  if (is.null(mustart)) {
        eval(family$initialize)
  }
  else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
  }
  if (EMPTY) {
      eta <- rep.int(0, nobs) + offset
      if (!valideta(eta)) 
            stop("Invalid linear predictor values in empty model")
      mu <- linkinv(eta)
      if (!validmu(mu)) 
          stop("Invalid fitted means in empty model")
      dev <- sum(dev.resids(y, mu, weights))
      w <- (weights * mu.eta(eta)^2)/variance(mu) ## incorrect for Newton
      residuals <- (y - mu)/mu.eta(eta)
      good <- rep(TRUE, length(residuals))
      boundary <- conv <- TRUE
      coef <- numeric(0)
      iter <- 0
      V <- variance(mu)
      alpha <- dev
      trA <- 0
      GCV <- nobs * alpha/(nobs - gamma * trA)^2
      UBRE <- alpha/nobs - scale + 2 * gamma/n * trA
      scale.est <- alpha/(nobs - trA)
      aic.model <- aic(y, n, mu, weights, dev) +  2 * trA
  } ### end if (EMPTY)
  else {
      eta <- if (!is.null(etastart)) 
          etastart
            else family$linkfun(mustart)
      mu <- as.numeric(linkinv(eta))
      if (!(validmu(mu) && valideta(eta))) 
            stop("Can't find valid starting values: please specify some")
      S.t <- matrix(0,q,q) # define the total sum of the penalties times sp
      n.pen <- length(S) # number of penalties 
      if (length(sp)!=n.pen) stop (paste("length of sp has to be equal to", n.pen))
      for (j in 1:n.pen) S.t <- S.t + sp[j]*S[[j]]
      # get sqrt of total penalty matrix...
      er <- eigen(S.t,symmetric=TRUE);  er$values[er$values<0] <- 0
      rS <- crossprod(sqrt(sqrt(er$values))*t(er$vectors))
      ii <- p.ident==1
      count <- sum(ii)
      iv <- array(0, dim=c(count,1)) # define an index vector for the monotone parameters
      iv <- (1:q)[ii]
      #----------------------------------------
      ## Initialization of parameters start here 
      beta0 <- get("start",envir=env)
      dbeta0 <- get("dbeta.start",envir=env)
      sp.old <- get("sp.last",envir=env)
      if (length(beta0)==0) {
          # list argument to pcls for initializing model coefficients
          M <- list(X=X,p=rep(0.1,q),C=matrix(0,0,0),sp=sp,y=eta-offset,w=y*0+1) 
          M$Ain <- matrix(0,q,q); diag(M$Ain) <- rep(1,q);
          M$bin <- rep(-1e+12,q); M$bin[iv] <- 1e-12
          M$off <- rep(0,n.pen); M$S <- list()
          for (j in 1:n.pen) {M$S[[j]] <- matrix(0,q,q); M$S[[j]] <- S[[j]]}
          beta.t <- pcls(M)      # initialize model coefficients (re-parameterized beta)
          beta <- beta.t         # initialize beta
          beta[iv] <- log(beta.t[iv]) # values of beta of the constrained terms
      }
      else {
          beta <- beta0
          beta.t <- beta               # current beta tilde
          beta.t[iv] <- exp(beta[iv])  # values of re-para beta of the constrained term
      }
      ## Initialization of parameters finishes here 
      #-------------------------------------------
      eta <- as.numeric(X%*%beta.t + offset)  # define initial linear predictor
      mu <- linkinv(eta)  # define initial fitted model
      dev <- sum(dev.resids(y,mu,weights)) # define initial norm/deviance
      pdev <- dev + sum((rS%*%beta)^2) # define initial penalized deviance 
      old.pdev <- pdev       # initialize convergence control for penalized deviance
      pdev.plot <- 0     # define initial pen dev for plotting it 
      E <- matrix(0,q,q)   # define diagonal matrix E- second term of the Hessian
      Cdiag <- rep(1,q);Cdiag[iv] <- beta.t[iv]
      C1diag <- rep(0,q);C1diag[iv] <- beta.t[iv]
      tX1 <- Cdiag*t(X)
      g.deriv <- 1/mu.eta(eta)        # diag(G)
      w1 <- weights/(variance(mu)*g.deriv^2)    # diag(W1)
      Dp.g <- - drop(tX1%*%(w1*g.deriv*(y-mu))) + S.t%*%beta
                     # the gradient vector of the penalized deviance
      Dp.gnorm <- max(abs(Dp.g)) # set convergence on the max value of the Dp.g
      old.beta <- beta
      conv <- FALSE
     # -----------------------------------------------
     # MAIN ITERATIONS START HERE --------------------
     #cat("\nscam.fit iter start")
     for (iter in 1:maxit)  {
         #cat(".")
         good <- weights > 0
         var.val <- variance(mu)
         varmu <- var.val[good]
         if (any(is.na(varmu))) 
              stop("NAs in V(mu)")
         if (any(varmu == 0)) 
              stop("0s in V(mu)")
         mu.eta.val <- mu.eta(eta)
         if (any(is.na(mu.eta.val[good]))) 
              stop("NAs in d(mu)/d(eta)")
         good <- (weights > 0) & (mu.eta.val != 0)
         if (all(!good)) {
              conv <- FALSE
              warning("No observations informative at iteration ", 
                  iter)
              break
         }
         Cdiag[iv] <- C1diag[iv] <- beta.t[iv]
         tX1 <- Cdiag*t(X)
         g.deriv <- 1/mu.eta(eta)  # diag(G)
         w1 <- weights/(variance(mu)*g.deriv^2)    # diag(W1) - Fisher weights
         y.mu <- y - mu
         alpha <- 1+ y.mu*(dv$dvar(mu)/variance(mu)+dg$d2link(mu)/g.deriv) # alpha elements of W
         w <- w1*alpha          # diag(W) - Newton weights
         diag(E) <- drop((C1diag*t(X))%*%(w1*g.deriv*y.mu))
         abs.w <- abs(w)      # absolute values of the diag(W)
         I.minus <- rep(0,nobs)  # define diagonal elements of the matrix I^{-}
         z1 <- g.deriv*y.mu/alpha  # the first term of the pseudodata
         ii <- w < 0;  I.minus[ii] <- 1;z1[ii] <- -z1[ii]
         wX11 <- rbind(sqrt(abs.w)[1:nobs]*t(tX1),rS)
           illcond <- FALSE
           Q <- qr(wX11,LAPACK=TRUE) 
           R <- qr.R(Q)
           rp <- 1:ncol(R)
           rp[Q$pivot] <- rp ## reverse pivot X=Q%*%R[,rp]
           if (Rrank(R)==ncol(R)) { ## no need to truncate, can just use QR
             R.inv <- backsolve(R,diag(ncol(R)))[rp,] ## inverse of unpivoted R
             tR.inv <- t(R.inv)
           } else { ## need SVD step
             R <- R[,rp] ## unpivoted R
             svd.r <- svd(R)
             d.inv <- rep(0,q)  # initial vector of inverse singular values
             good <- svd.r$d >= max(svd.r$d)*.Machine$double.eps^.5
             d.inv[good] <- 1/svd.r$d[good]
             if (sum(!good)>0) illcond <- TRUE
             R <- svd.r$d*t(svd.r$v)
             Q <- qr.qy(Q,rbind(svd.r$u,matrix(0,nobs,q))) ## this is inefficient don't need to extract Q really
             tR.inv <- d.inv*t(svd.r$v)    # inverse of transpose of R
             R.inv <- t(tR.inv)
           }
           
          QtQRER <- tR.inv%*%(diag(E)*R.inv)
          if (sum(I.minus)>0) {
             if (is.qr(Q)) { 
            # QtQRER <- QtQRER + 2*tcrossprod(qr.qty(Q,diag(nrow(wX11))[,(1:nobs)[as.logical(I.minus)]]))
             QtQRER <- QtQRER + 2*crossprod(I.minus*qr.Q(Q)[1:nobs,])  
             } else {
             QtQRER <- QtQRER + 2*crossprod(I.minus*Q[1:nobs,])
             }
         }
         ei <- eigen(QtQRER,symmetric=TRUE)    # eigen-decomposition of QtQRER
         d <- ei$values        # vector of eigenvalues
         ok1 <- sum(d>1) > 0 # checking positive semi-definiteness 
         if (ok1 == TRUE) {# Fisher step in case of not +ve semi-definiteness of penalized loglikelihood
                         # set alpha =1
             eta.t <- drop(t(beta)%*%tX1)     # eta tilde for pseudodata 
             wX11 <- rbind(sqrt(w1)[1:nobs]*t(tX1),rS)   
             z<-g.deriv*y.mu+eta.t      # pseudodata
             wz<-w1^.5*z               # weighted pseudodata
             wz.aug<-c(wz,rep(0,nrow(rS)))   # augmented pseudodata
             Q <- qr(wX11,LAPACK=TRUE) 
             R <- qr.R(Q)
             rp <- 1:ncol(R)
             rp[Q$pivot] <- rp ## reverse pivot X=Q%*%R[,rp]
             if (Rrank(R)==ncol(R)) { ## no need to truncate, can just use QR
               beta <- backsolve(R,qr.qty(Q,wz.aug)[1:q])[rp]
             } else { ## need SVD step
               R <- R[,rp] ## unpivoted R
               s1 <- svd(R)
               d.inv1 <- rep(0,q)
               good1 <- s1$d >= max(s1$d)*.Machine$double.eps^.5
               d.inv1[good1] <- 1/s1$d[good1]
               beta <- s1$v%*%((d.inv1*t(s1$u))%*%qr.qty(Q,wz.aug)[1:q])
             }
         }  ### end of if (ok1) - Fisher step 
         else  {       ##  full Newton step
             Id.inv.r<-1/(1-d)^.5   # vector of inverse values of (1-sv)^.5
             ii <- (1-d) < .Machine$double.eps
             Id.inv.r[ii] <- 0
             eidrop <- t(Id.inv.r*t(ei$vectors))
             wz1<-sqrt(abs.w)*z1  # the first term of the weighted pseudodata
             if (is.qr(Q)) {
               beta <- R.inv%*%(eidrop%*%(t(eidrop)%*%qr.qty(Q,c(wz1,rep(0,nrow(rS))))[1:nrow(eidrop)]))
             } else {
               beta <- R.inv%*%(eidrop%*%(t(eidrop)%*%(t(Q[1:nobs,])%*%wz1)[1:nrow(eidrop)]))
             }
             beta <- old.beta + 
                     drop(beta - R.inv%*%(eidrop%*%(t(eidrop)%*%(tR.inv%*%(S.t%*%old.beta)))))
         }  ###  end of if (!ok1) - Newton step
         delta <- beta-c(old.beta)         # trial step
         step <- 1                      # initial trial step length
         beta <- c(old.beta)+step*delta    # current parameter estimates
         beta.t <- beta                 # current reparameterized beta
         beta.t[iv] <- exp(beta[iv])  # values of re-para beta of the shape constrained term
         eta <- as.numeric(X%*%beta.t + offset)     # linear predictor
         mu <- linkinv(eta)          # fitted values
         dev <- sum(dev.resids(y,mu,weights)) # deviance of the working model
         pdev <- dev + sum((rS%*%beta)^2) # deviance + penalty of the working model
               
         ## `step reduction' approach starts here ---------------------
         ii <- 1 
         div.thresh <- 10*(.1 +abs(old.pdev))*.Machine$double.eps^.5
         while (is.na(pdev) || (pdev-old.pdev) > div.thresh) { # 'step reduction' approach
             if (ii > 200) 
                stop ("step reduction failed")
             ii <- ii+1
             step <- step/2         # decrease step length 
             beta <- c(old.beta)+step*delta   # update current parameter estimates
             beta.t <- beta                # update current re-para beta
             beta.t[iv] <- exp(beta[iv])   
             eta <- as.numeric(X%*%beta.t + offset)    # linear predictor  
             mu <- linkinv(eta)         # fitted values
             dev <- sum(dev.resids(y,mu,weights)) # update deviance of the working model
             pdev <- dev+sum((rS%*%beta)^2) # update pen deviance of the working model
         }
         ## `step reduction' finishes here -----------------------
     
         Dp.g <- -drop(tX1%*%(w1*g.deriv*(y-mu)))+S.t%*%beta # the gradient vector of the penalized deviance
         Dp.gnorm <- max(abs(Dp.g)) 
         pdev.plot[iter] <- pdev      # store penilized deviance of the working model for plotting
          
         ## checking convergence .......
         if (abs(pdev - old.pdev)/(.1 + abs(pdev)) < devtol) {
             if (max(abs(beta - c(old.beta))) > steptol * 
                          max(abs(beta + c(old.beta)))/2) {
                old.beta <- beta
                old.pdev <- pdev
             }
             else {
                conv <- TRUE
                beta <- beta
                break
             }
         }
         else {
             old.pdev <- pdev
             old.beta <- beta
         }
     } ## main iteration procedure is completed here ------------
     ##______________________________________________________________
     ## define matrices at their converged values from the full Newton method------
  
     dev <- sum(dev.resids(y,mu,weights))
     beta.t <- beta               # estimates of re-para beta
     beta.t[iv] <- exp(beta[iv])   
     eta <- as.numeric(X%*%beta.t + offset)      # linear predictor  
     mu <- linkinv(eta)         # fitted values
     Cdiag[iv] <- C1diag[iv] <- beta.t[iv]
     X1 <- t(Cdiag*t(X)) 
     g.deriv <- 1/ mu.eta(eta)        # diag(G)
     w1 <- weights/(variance(mu)*g.deriv^2)    # diag(W1)
     alpha <- 1+(y-mu)*(dv$dvar(mu)/variance(mu)+dg$d2link(mu)/g.deriv) # alpha elements of W
     w <- w1*alpha          # diag(W)
     diag(E) <- drop((C1diag*t(X))%*%(w1*g.deriv*(y-mu))) # diagonal elements of E
     abs.w <- abs(w)      # absolute values of the diag(W)
     I.minus <- rep(0,nobs)  # define diagonal elements of the matrix I^{-}
     I.minus[w<0] <- 1
     wX1 <- sqrt(abs.w)[1:nobs]*X1 ## wX1 actually used later
     wX11 <- rbind(wX1,rS) # augmented model matrix 
      ## Faster version only does SVD when needed (and then only on QR factor)
        illcond <- FALSE
        Q <- qr(wX11,LAPACK=TRUE) 
        R <- qr.R(Q)
        rp <- 1:ncol(R)
        rp[Q$pivot] <- rp ## reverse pivot X=Q%*%R[,rp]
 
        R.out <- R[,rp]  ## unpivoted R, needed for summary function

        if (Rrank(R)==ncol(R)) { ## no need to truncate, can just use QR
           R.inv <- backsolve(R,diag(ncol(R)))[rp,] ## inverse of unpivoted R
           tR.inv <- t(R.inv)
        } else { ## need SVD step
           R <- R[,rp] ## unpivoted R
           svd.r <- svd(R)
           d.inv <- rep(0,q)  # initial vector of inverse singular values
           good <- svd.r$d >= max(svd.r$d)*.Machine$double.eps^.5
           d.inv[good] <- 1/svd.r$d[good]
           if (sum(!good)>0) illcond <- TRUE
           R <- svd.r$d*t(svd.r$v)
           Q <- qr.qy(Q,rbind(svd.r$u,matrix(0,nobs,q))) ## this is inefficient don't need to extract Q really
           tR.inv <- d.inv*t(svd.r$v)    # inverse of transpose of R
           R.inv <- t(tR.inv)
         }
       QtQRER <- tR.inv%*%(diag(E)*R.inv)
         if (sum(I.minus)>0) {
              if (is.qr(Q)) { 
              QtQRER <- QtQRER + 2*crossprod(I.minus*qr.Q(Q)[1:nobs,])  
              } else {
              QtQRER <- QtQRER + 2*crossprod(I.minus*Q[1:nobs,])
              }
         }
      ei <- eigen(QtQRER,symmetric=TRUE)    # eigen-decomposition of QtQRER
      d <- ei$values        # vector of eigenvalues
      ok1 <- sum(d>1)>0
      if (ok1) { ## Fisher step in case of not positive semi-definiteness of penalized loglikelihood
                 ## set alpha =1 ...
         wX1<-sqrt(w1)[1:nobs]*X1
         wX11<-rbind(wX1,rS)
         Q <- qr(wX11,LAPACK=TRUE) 
         R <- qr.R(Q)
         rp <- 1:ncol(R)
         rp[Q$pivot] <- rp ## reverse pivot X=Q%*%R[,rp]

         R.out <- R[,rp]  ## unpivoted R, needed for summary function

         if (Rrank(R)==ncol(R)) { ## no need to truncate, can just use QR
             P <- backsolve(R,diag(ncol(R)))[rp,]
             K <- qr.Q(Q)[1:nobs,]
         } else { ## need SVD step
             R <- R[,rp] ## unpivoted R
             s1 <- svd(R)
             d.inv1 <- rep(0,q)
             good1 <- s1$d >= max(s1$d)*.Machine$double.eps^.5
             d.inv1[good1] <- 1/s1$d[good1]
             P <- t(d.inv1*t(s1$v))
             K <- qr.qy(Q,rbind(s1$u,matrix(0,nobs,q)))[1:nobs,]
          }
     } ## end of if (ok1)
     else  {   ## full Newton step
         Id.inv.r<-1/(1-d)^.5   # vector of inverse values of (1-sv)^1/2
         ii <- (1-d) < .Machine$double.eps
         Id.inv.r[ii] <- 0
         eidrop <- t(Id.inv.r*t(ei$vectors))
         P <- R.inv%*%eidrop  ## ei$vectors%*%diag(Id.inv.r) # define matrix P
         if (is.qr(Q)) {
           K <- qr.qy(Q,rbind(eidrop,matrix(0,nobs,q)))[1:nobs,]
         } else {
           K <- Q[1:nobs,]%*%eidrop ## (ei$vectors%*%diag(Id.inv.r))  # define matrix K 
         }
     }  ## end of if (!ok1)
     # end of calculation of the matrices at their converged values ------
     # -------------------------------------------------------------
     Dp.g <- -t(X1)%*%(w1*g.deriv*(y-mu))+S.t%*%beta # the gradient vector of the penalized deviance
     Dp.gnorm<-max(abs(Dp.g)) 

     # calculating tr(A) ----------------------------------------------
     I.plus <- rep(1,nobs)   # define diagonal elements of the matrix I^{+}
     I.plus[w<0] <- -1
     L <- c(1/alpha)    # define diagonal elements of L=diag(1/alpha)
     ## NOTE PKt is O(np^2) and not needed --- can compute trA as side effect of gradiant
     KtILQ1R <- t(L*I.plus*K)%*%wX1
     edf <- rowSums(P*t(KtILQ1R))
     trA <- sum(edf)

     # ---------------------------------------------------------------------------
     scale.est <- dev/(nobs-trA)  #  scale estimate...
     residuals <- rep.int(NA, nobs)
     residuals <- (y-mu)*g.deriv

     # ---------------------------------------------------------------------------
     ## calculation of the derivatives of beta by the Implicit Function Theorem starts here
     dbeta.rho <- matrix(0,q,n.pen) # define matrix of the parameters derivatives
     for (j in 1:n.pen) {
        dbeta.rho[,j] <- - sp[j]*P%*%(t(P)%*%(S[[j]]%*%beta)) # derivative of beta wrt rho[j]
     }
     # end of calculating the parameters derivatives
     aic.model <- aic(y, n, mu, weights, dev) +  2 * sum(edf)
     assign("start",beta,envir=env)
     assign("dbeta.start",dbeta.rho,envir=env)
     assign("sp.last",sp,envir=env)
  } ### end if (!EMPTY) 

 list(L=L,C1diag=C1diag,E=E,iter=iter, old.beta=old.beta, step=step,gcv=dev*nobs/(nobs-trA)^2,
      sp=sp, mu=mu,X=X, X1=X1,beta=beta,beta.t=beta.t,iv=iv,S=S,S.t=S.t,rS=rS,
      P=P,K=K, KtILQ1R= KtILQ1R,dlink.mu=1/mu.eta(eta),Var=variance(mu), abs.w=abs.w,
      link=family$linkfun(mu),w=as.numeric(w),w1=w1,d2link.mu=dg$d2link(mu),wX1=wX1,I.plus=I.plus,
      dvar.mu=dv$dvar(mu),d2var.mu=dv$d2var(mu),deviance=dev,scale.est=scale.est,
      ok1=ok1,alpha=as.numeric(alpha),d3link.mu=dg$d3link(mu),eta=eta,iter=iter,
      Dp.gnorm=Dp.gnorm, Dp.g=Dp.g,d=d, conv=conv, illcond=illcond,R=R.out, edf=edf,trA=trA,
      residuals=residuals,z=g.deriv*(y-mu)+X1%*%beta,dbeta.rho=dbeta.rho, aic=aic.model)
} ## end of scam.fit




#######################################################################
## function to get null deviance and covariance matrices after fit   ##
#######################################################################


scam.fit.post <- function(y,X,object,sig2,offset,intercept,
                         weights,scale.known)
{  ## Function to compute null deviance and covariance matrices after a scam fit.
   ## covariance matrix should use expected Hessian, so re-computation of factors 
   ## is required.  
   ## object - object from estimate.scam()
   n <- nobs <- NROW(y) # number of observations
   linkinv <- object$family$linkinv
   dev.resids <- object$family$dev.resids
   
   wtdmu <- if (intercept) sum(weights * y)/sum(weights) 
              else linkinv(offset)

   nulldev <- sum(dev.resids(y, wtdmu, weights))
 
   n.ok <- nobs - sum(weights == 0)
   nulldf <- n.ok - as.integer(intercept)
   
   # calculating the approximate covariance matrices 
   # (dealing with the expected Hessian of the log likelihood) ...

   if (!scale.known) sig2 <- object$scale.est
   ## get the inverse of the expected Hessian...
   wX1 <- sqrt(object$w1)[1:n]*object$X1
   wX11 <- rbind(wX1,object$rS)
   q <- ncol(wX1)
   Q <- qr(wX11,LAPACK=TRUE) 
   R <- qr.R(Q)
   rp <- 1:ncol(R)
   rp[Q$pivot] <- rp ## reverse pivot X=Q%*%R[,rp]
   if (Rrank(R)==ncol(R))  ## no need to truncate, can just use QR
      {  P <- backsolve(R,diag(q))[rp,]
         K <- qr.Q(Q)[1:n,]
      } else {  ## need SVD step
                R <- R[,rp] ## unpivoted R
                s1 <- svd(R)
                d.inv1 <- rep(0,q)
                good1 <- s1$d >= max(s1$d)*.Machine$double.eps^.5
                d.inv1[good1] <- 1/s1$d[good1]
                P <- t(d.inv1[good1]*t(s1$v[,good1]))
                K <- qr.qy(Q,rbind(s1$u,matrix(0,n,q))[,good1])[1:n,]         
             }
   Vb <- tcrossprod(P) * sig2 
          ## P%*%t(P)*sig2 # Bayesian posterior covariance matrix for the parameters 
   Ve <- crossprod(K%*%t(P)) *sig2
        #PKt%*%t(PKt)*sig2 # covariance matrix of the parameter estimators 
   ## Delta method to get covariance matrix for the reparametrized parameters...
   df.p <- rep(1,q)
   df.p[object$iv] <- object$beta.t[object$iv]
   Vb.t <- t(df.p*t(df.p*Vb))
   Ve.t <- t(df.p*t(df.p*Ve))

   ## calculating edf and trA...
   KtILQ1R <- t(object$L*object$I.plus*K)%*%wX1
   F <- P%*%(KtILQ1R)
   edf <- diag(F) ## effective degrees of freedom
   edf1 <- 2*edf - rowSums(t(F)*F) ## alternative
   trA <- sum(edf)

   list (nulldev=nulldev, df.null=nulldf,Vb=Vb,Vb.t=Vb.t,Ve=Ve,Ve.t=Ve.t,
        sig2=sig2,edf=edf,edf1=edf1,trA=trA)
}


###############################################################
## loading functions, copied from mgcv() package of Simon Wood
#################################################################

print.scam.version <- function()
{ library(help=scam)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  hello <- paste("This is scam ",version,".",sep="")
  packageStartupMessage(hello)
}


.onAttach <- function(...) { 
  print.scam.version()
}

##.onUnload <- function(libpath) library.dynam.unload("scam", libpath)






