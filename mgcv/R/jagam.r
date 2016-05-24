## (c) Simon Wood 2014. Released under GPL2. 
## jagam code (Just Another Gibbs Additive Model)
## Code offering JAGS/BUGS support for mgcv.
## In particular autogenerates the code and data to fit an mgcv
## style GAM in JAGS, and re-packages the simulation output
## in a form suitable for plotting and prediction. 
## Idea is that the code would be modified to add the sort
## of random effects structure most appropriately handled in JAGS.


write.jagslp <- function(resp,family,file,use.weights,offset=FALSE) {
## write the JAGS code for the linear predictor  
## and response distribution. 
  iltab <- ## table of inverse link functions
    c("eta[i]","exp(eta[i])","ilogit(eta[i])","1/eta[i]","eta[i]^2")
  names(iltab) <- c("identity","log","logit","inverse","sqrt")
  if (!family$link%in%names(iltab)) stop("sorry link not yet handled")
  
  ## code linear predictor and expected response...
  if (family$link=="identity") {
    if (offset) cat("  mu <- X %*% b + offset ## expected response\n",file=file,append=TRUE)
    else cat("  mu <- X %*% b ## expected response\n",file=file,append=TRUE)
  } else {
    if (offset) cat("  eta <- X %*% b + offset ## linear predictor\n",file=file,append=TRUE)
    else cat("  eta <- X %*% b ## linear predictor\n",file=file,append=TRUE)
    cat("  for (i in 1:n) { mu[i] <- ",iltab[family$link],"} ## expected response\n",file=file,append=TRUE)
  }
  ## code the response given mu and any scale parameter prior...
  #scale <- TRUE ## is scale parameter free?
  cat("  for (i in 1:n) { ",file=file,append=TRUE)
  if (family$family=="gaussian") {
    if (use.weights) cat(resp,"[i] ~ dnorm(mu[i],tau*w[i]) } ## response \n",sep="",file=file,append=TRUE)
    else cat(resp,"[i] ~ dnorm(mu[i],tau) } ## response \n",sep="",file=file,append=TRUE)
    cat("  scale <- 1/tau ## convert tau to standard GLM scale\n",file=file,append=TRUE) 
    cat("  tau ~ dgamma(.05,.005) ## precision parameter prior \n",file=file,append=TRUE)
  } else if (family$family=="poisson") {
   # scale <- FALSE
    cat(resp,"[i] ~ dpois(mu[i]) } ## response \n",sep="",file=file,append=TRUE)
    if (use.weights) warning("weights ignored") 
    use.weights <- FALSE
  } else if (family$family=="binomial") {
   # scale <- FALSE
    cat(resp,"[i] ~ dbin(mu[i],w[i]) } ## response \n",sep="",file=file,append=TRUE)
    use.weights <- TRUE
  } else if (family$family=="Gamma") {
    if (use.weights) cat(resp,"[i] ~ dgamma(r*w[i],r*w[i]/mu[i]) } ## response \n",sep="",file=file,append=TRUE)
    else cat(resp,"[i] ~ dgamma(r,r/mu[i]) } ## response \n",sep="",file=file,append=TRUE)
    cat("  r ~ dgamma(.05,.005) ## scale parameter prior \n",file=file,append=TRUE)
    cat("  scale <- 1/r ## convert r to standard GLM scale\n",file=file,append=TRUE) 
  } else stop("family not implemented yet")
  use.weights
} ## write.jagslp

jini <- function(G,lambda) {
## get initial coefficients to initialize JAGS, otherwise
## initialization is hit and miss. 
  y <- G$y; nobs <- length(y); p <- ncol(G$X)
  family <- G$family
  weights <- G$w
  start <- mustart <- etastart <- NULL ## ignore codetools warning - needed for eval
  eval(G$family$initialize)
  w <- as.numeric(G$w * family$mu.eta(family$linkfun(mustart))^2/family$variance(mustart)) 
  w <- sqrt(w)
  z <- c(w*family$linkfun(mustart),rep(0,p)) ## residual is zero, so eta is all there is!
  X <- rbind(w*G$X,matrix(0,p,p)) 
  ## now append square roots of penalties
  uoff <- unique(G$off)
  for (i in 1:length(uoff)) {
    jj <- which(G$off%in%uoff[i])
    S <- G$S[[jj[1]]]*lambda[[jj[1]]]
    m <- length(jj)
    if (m>1) for (j in jj) S <- S +  G$S[[j]]*lambda[j]
    S <- t(mroot(S))
    jj <- nrow(S)
    X[(nobs+1):(nobs+jj),uoff[i]:(uoff[i]+ncol(S)-1)] <- S
    nobs <- nobs + jj
  }
  ## we need some idea of initial coeffs and some idea of 
  ## associated standard error...
  qrx <- qr(X,LAPACK=TRUE)
  rp <- qrx$pivot;rp[rp] <- 1:ncol(X)
  Ri <- backsolve(qr.R(qrx),diag(1,nrow=ncol(X)))[rp,] 
  beta <- qr.coef(qrx,z)
  se <- sqrt(rowSums(Ri^2))*sqrt(sum((z-X%*%beta)^2)/nrow(X))
  list(beta=beta,se=se)
} ## jini

jagam <- function(formula,family=gaussian,data=list(),file,weights=NULL,na.action,
offset=NULL,knots=NULL,sp=NULL,drop.unused.levels=TRUE,control=gam.control(),centred=TRUE,
sp.prior = "gamma",diagonalize=FALSE) {
## rho contains log smoothing params and b the model coefficients, in JAGS
## diagonalize==TRUE actually seems to be faster for high dimensional terms
## in the Gaussian setting (Conjugate updates better than MH), otherwise 
## diagonalize==FALSE faster as block MH is highly advantageous
## WARNING: centred=FALSE is usually a very bad idea!!
  if (is.null(file)) stop("jagam requires a file for the JAGS model specification")
  cat("model {\n",file=file) ## start the model specification
  if (!(sp.prior %in% c("gamma","log.uniform"))) {
    warning("smoothing parameter prior choise not recognised, reset to gamma")
  }
  ## takes GAM formula and data and produces JAGS model and corresponding 
  ## data list...
  if (is.character(family))
            family <- eval(parse(text = family))
  if (is.function(family))
            family <- family()
  if (is.null(family$family))
            stop("family not recognized")
 
  gp <- interpret.gam(formula) # interpret the formula 
  cl <- match.call() # call needed in gam object for update to work
  mf <- match.call(expand.dots=FALSE)
  mf$formula <- gp$fake.formula 
  mf$family <- mf$knots <- mf$sp <- mf$file <- mf$control <- 
  mf$centred <- mf$sp.prior <- mf$diagonalize <- NULL
  mf$drop.unused.levels <- drop.unused.levels
  mf[[1]]<-as.name("model.frame")
  pmf <- mf
 
  pmf$formula <- gp$pf
  pmf <- eval(pmf, parent.frame()) # pmf contains all data for parametric part
  pterms <- attr(pmf,"terms") ## pmf only used for this
  rm(pmf)
 
  mf <- eval(mf, parent.frame()) # the model frame now contains all the data 
  if (nrow(mf)<2) stop("Not enough (non-NA) data to do anything meaningful")
  terms <- attr(mf,"terms")

  ## summarize the *raw* input variables
  ## note can't use get_all_vars here -- buggy with matrices
  vars <- all.vars(gp$fake.formula[-2]) ## drop response here
  inp <- parse(text = paste("list(", paste(vars, collapse = ","),")"))

  ## allow a bit of extra flexibility in what `data' is allowed to be (as model.frame actually does)
  if (!is.list(data)&&!is.data.frame(data)) data <- as.data.frame(data) 

  dl <- eval(inp, data, parent.frame())
  if (!control$keepData) { rm(data)} ## save space
  names(dl) <- vars ## list of all variables needed
  var.summary <- variable.summary(gp$pf,dl,nrow(mf)) ## summarize the input data
  rm(dl)   

  G <- gam.setup(gp,pterms=pterms,
                 data=mf,knots=knots,sp=sp,
                 H=NULL,absorb.cons=centred,sparse.cons=FALSE,select=TRUE,
                 idLinksBases=TRUE,scale.penalty=control$scalePenalty,
                 diagonal.penalty=diagonalize)
  G$model <- mf;G$terms <- terms;G$family <- family;G$call <- cl
  G$var.summary <- var.summary
  ## write JAGS code producing linear predictor and linking linear predictor to 
  ## response....

  use.weights <- if (is.null(weights)) FALSE else TRUE 
  use.weights <- write.jagslp("y",family,file,use.weights,!is.null(G$offset))
  if (is.null(weights)&&use.weights) weights <- rep(1,nrow(G$X))  

  ## start the JAGS data list...

  jags.stuff <- list(y=G$y,n=length(G$y),X=G$X)  
  if (!is.null(G$offset)) jags.stuff$offset <- G$offset
  if (use.weights) jags.stuff$w <- weights

  if (family$family == "binomial") jags.stuff$y <- G$y*weights ## JAGS not expecting observed prob!!
 
  ## get initial values, for use by JAGS, and to guess suitable values for
  ## uninformative priors...

  lambda <- initial.spg(G$X,G$y,G$w,family,G$S,G$off,G$L) ## initial sp values
  jags.ini <- list()
  lam <- if (is.null(G$L)) lambda else G$L%*%lambda
  jin <- jini(G,lam)
  jags.ini$b <- jin$beta
  prior.tau <- signif(0.01/(abs(jin$beta) + jin$se)^2,2)

  ## set the fixed effect priors...
  if (G$nsdf>0) {
    ptau <- min(prior.tau[1:G$nsdf]) 
    cat("  ## Parametric effect priors CHECK tau=1/",signif(1/sqrt(ptau),2),"^2 is appropriate!\n",file=file,append=TRUE,sep="")
    cat("  for (i in 1:",G$nsdf,") { b[i] ~ dnorm(0,",ptau,") }\n",file=file,append=TRUE,sep="")
  }

  ## Work through smooths.
  ## In JAGS terms the penalties should simply define priors.
  ## Any unpenalized term should be given a diffuse prior.  
  ## For diagonalized terms these should be written directly into the code
  ## and there is nothing to pass to JAGS.
  ## For overlapping multi term penalties, a null space penalty needs to
  ## be added and the components of the penalty have to be passed into 
  ## JAGS in the argument list: cbinding the components into one matrix seems sensible.
  ## Smoothing parameters should be in a single vector in the code indexed by 
  ## number.  
  n.sp <- 0 ## count the smoothing parameters....
  for (i in 1:length(G$smooth)) {
    ## Are penalties seperable...
    seperable <- FALSE
    M <- length(G$smooth[[i]]$S)
    p <- G$smooth[[i]]$last.para - G$smooth[[i]]$first.para + 1 ## number of params
    if (M<=1) seperable <- TRUE else {
      overlap <- rowSums(G$smooth[[i]]$S[[1]])
      for (j in 2:M) overlap <- overlap & rowSums(G$smooth[[i]]$S[[j]])
      if (!sum(overlap)) seperable <- TRUE 
    }
    if (seperable) { ## double check that they are diagonal
      if (M>0) for (j in 1:M) {
        if (max(abs(G$smooth[[i]]$S[[j]] - diag(diag(G$smooth[[i]]$S[[j]]),nrow=p)))>0) seperable <- FALSE
      } 
    }
    cat("  ## prior for ",G$smooth[[i]]$label,"... \n",file=file,append=TRUE,sep="")
    if (seperable) {
      b0 <- G$smooth[[i]]$first.para
      if (M==0) {
        cat("  ## Note fixed vague prior, CHECK tau = 1/",signif(1/sqrt(ptau),2),"^2...\n",file=file,append=TRUE,sep="")
        b1 <- G$smooth[[i]]$last.para
        ptau <- min(prior.tau[b0:b1])
        cat("  for (i in ",b0,":",b1,") { b[i] ~ dnorm(0,",ptau,") }\n",file=file,append=TRUE,sep="")
      } else for (j in 1:M) {
        D <- diag(G$smooth[[i]]$S[[j]]) > 0
        b1 <- sum(as.numeric(D)) + b0 - 1
        n.sp <- n.sp + 1
        cat("  for (i in ",b0,":",b1,") { b[i] ~ dnorm(0, lambda[",n.sp,"]) }\n",file=file,append=TRUE,sep="")
        b0 <- b1 + 1
      }
    } else { ## inseperable - requires the penalty matrices to be supplied to JAGS... 
      b0 <- G$smooth[[i]]$first.para; b1 <- G$smooth[[i]]$last.para
      Kname <- paste("K",i,sep="") ## total penalty matrix in JAGS
      Sname <- paste("S",i,sep="") ## components of total penalty in R & JAGS
      cat("  ",Kname," <- ",Sname,"[1:",p,",1:",p,"] * lambda[",n.sp+1,"] ",
          file=file,append=TRUE,sep="")
      if (M>1) { ## code to form total precision matrix...  
        for (j in 2:M) cat(" + ",Sname,"[1:",p,",",(j-1)*p+1,":",j*p,"] * lambda[",n.sp+j,"]",
            file=file,append=TRUE,sep="")
      }
      cat("\n  b[",b0,":",b1,"] ~ dmnorm(zero[",b0,":",b1,"],",Kname,") \n"
           ,file=file,append=TRUE,sep="")
      n.sp <- n.sp + M
      Sc <- G$smooth[[i]]$S[[1]]
      if (M>1) for (j in 2:M) Sc <- cbind(Sc,G$smooth[[i]]$S[[j]])
      jags.stuff[[Sname]] <- Sc
      jags.stuff$zero <- rep(0,ncol(G$X))
    }
  } ## smoothing penalties finished

  ## Write the smoothing parameter prior code, using L if it exists.

  cat("  ## smoothing parameter priors CHECK...\n",file=file,append=TRUE,sep="")
  if (is.null(G$L)) {
    if (sp.prior=="log.uniform") {
      cat("  for (i in 1:",n.sp,") {\n",file=file,append=TRUE,sep="")
      cat("    rho[i] ~ dunif(-12,12)\n",file=file,append=TRUE,sep="") 
      cat("    lambda[i] <- exp(rho[i])\n",file=file,append=TRUE,sep="")
      cat("  }\n",file=file,append=TRUE,sep="")
      jags.ini$rho <- log(lambda)
    } else { ## gamma priors
      cat("  for (i in 1:",n.sp,") {\n",file=file,append=TRUE,sep="")
      cat("    lambda[i] ~ dgamma(.05,.005)\n",file=file,append=TRUE,sep="") 
      cat("    rho[i] <- log(lambda[i])\n",file=file,append=TRUE,sep="")
      cat("  }\n",file=file,append=TRUE,sep="")
      jags.ini$lambda <- lambda
    }
  } else { 
    jags.stuff$L <- G$L
    rho.lo <- FALSE
    if (any(G$lsp0!=0)) {
      jags.stuff$rho.lo <- G$lsp0
      rho.lo <- TRUE
    }
    nr <- ncol(G$L)
    if (sp.prior=="log.uniform") {
      cat("  for (i in 1:",nr,") { rho0[i] ~ dunif(-12,12) }\n",file=file,append=TRUE,sep="")
      if (rho.lo) cat("  rho <- rho.lo + L %*% rho0\n",file=file,append=TRUE,sep="")
      else cat("  rho <- L %*% rho0\n",file=file,append=TRUE,sep="")
      cat("  for (i in 1:",n.sp,") { lambda[i] <- exp(rho[i]) }\n",file=file,append=TRUE,sep="")
      jags.ini$rho0 <- log(lambda)
    } else { ## gamma prior
      cat("  for (i in 1:",nr,") {\n",file=file,append=TRUE,sep="")
      cat("    lambda0[i] ~ dgamma(.05,.005)\n",file=file,append=TRUE,sep="") 
      cat("    rho0[i] <- log(lambda0[i])\n",file=file,append=TRUE,sep="")
      cat("  }\n",file=file,append=TRUE,sep="")
      if (rho.lo) cat("  rho <- rho.lo + L %*% rho0\n",file=file,append=TRUE,sep="")
      else cat("  rho <- L %*% rho0\n",file=file,append=TRUE,sep="")
      cat("  for (i in 1:",n.sp,") { lambda[i] <- exp(rho[i]) }\n",file=file,append=TRUE,sep="")
      jags.ini$lambda0 <- lambda
    }
  } 
  cat("}",file=file,append=TRUE)

  G$formula=formula
  G$rank=ncol(G$X) ## to Gibbs sample we force full rank!
  list(pregam=G,jags.data=jags.stuff,jags.ini=jags.ini)
} ## jagam

sim2jam <- function(sam,pregam,edf.type=2,burnin=0) {
## takes jags simulation output with field, b, containing model coefficients
## and a pregam object from jagam, and attempts to create a fake gam object suitable
## for plotting. This is given a class "jam" since only a limited range of gam 
## methods are appropriate for such models. Ideally...
## vcov, print, plot, predict, model.matrix, ... 
  if (is.null(sam$b)) stop("coefficient simulation data is missing") 
  if (burnin>0) {
    nc <- dim(sam$b)[2] ## chain length
    if (burnin >= nc*.9) {
      warning("burnin too large, reset")
      burnin <- min(nc-1,floor(nc * .9))
    } 
    ind <- (burnin+1):nc
    sam$b <- sam$b[,ind,]
    if (!is.null(sam$mu)) sam$mu <- sam$mu[,ind,]
    if (!is.null(sam$rho)) sam$rho <- sam$rho[,ind,]
    if (!is.null(sam$scale)) sam$scale <- sam$scale[,ind,]
  }
  pregam$Vp <- cov(t(sam$b[,,1]))
  pregam$coefficients <- rowMeans(sam$b[,,1])
  pregam$sig2 <- if (is.null(sam$scale)) 1 else mean(sam$scale)
  n.chain <- dim(sam$b)[3]
  if (n.chain>1) { 
    for (i in 2:n.chain) {
      pregam$Vp <- pregam$Vp +  cov(t(sam$b[,,i]))
      pregam$coefficients <-  pregam$coefficients + rowMeans(sam$b[,,i])
    }
    pregam$Vp <- pregam$Vp/n.chain
    pregam$coefficients <-  pregam$coefficients/n.chain
  }
  ## NOTE: 3 edf versions... 
  ##       0. diag((X'X+S)^{-1}X'X)
  ##       1. diag((X'WX+S)^-1X'WX)
  ##       2. diag(VbX'WX)/scale Vb by simulation. mu used for W may also be by sim.
  if (edf.type<2&&is.null(sam$rho)) {
    edf.type <- 2
    warning("rho missing from simulation data edf.type reset to 2")
  }
  if (edf.type > 0) { ## use X'WX not X'X
    if (is.null(sam$mu)) {
      eta <- pregam$X %*% pregam$coefficients
      mu <- pregam$family$linkinv(eta)
    } else { 
      mu <- rowMeans(sam$mu)
      eta <- pregam$family$linkfun(mu)
    }
    w <- as.numeric(pregam$w * pregam$family$mu.eta(eta)^2/pregam$family$variance(mu))
    XWX <- t(pregam$X) %*% (w*pregam$X) 
  } else XWX <- t(pregam$X) %*% (pregam$X) 
  if (edf.type < 2) { ## tr((X'WX + S)^{-1}X'WX
    rho <- rowMeans(sam$rho);lambda <- exp(rho)
    XWXS <- XWX
    for (i in 1:length(lambda)) {
      ind <- pregam$off[i]:(pregam$off[i]+ncol(pregam$S[[i]])-1)
      XWXS[ind,ind] <-  XWXS[ind,ind] + pregam$S[[i]] * lambda[i]
    } 
    pregam$edf <-  diag(solve(XWXS,XWX))
  } else pregam$edf <- rowSums(pregam$Vp*t(XWX))/pregam$sig2 ## tr(Vb%*%XWX)/scale
  class(pregam) <- "jam"
  pregam
} ## sim2jam

## method functions. Simple wrappers for gam methods
## idea is to limit options to those generally computable...

print.jam <- function(x,...) print.gam(x,...)
vcov.jam <- function(object,...) vcov.gam(object,...)
plot.jam <- function(x,rug=TRUE,se=TRUE,pages=0,select=NULL,scale=-1,
              n=100,n2=40,pers=FALSE,theta=30,phi=30,jit=FALSE,xlab=NULL,
              ylab=NULL,main=NULL,ylim=NULL,xlim=NULL,too.far=0.1,
              shade=FALSE,shade.col="gray80",
              shift=0,trans=I,seWithMean=FALSE,
              scheme=0,...) {
  ## residuals, unconditional, by.resids and all.terms not supported...
  arg.names <- names(list(...))
  if (length(arg.names)>0) {
    if ("residuals"%in% arg.names) stop("residuals argument not supported")
    if ("unconditional"%in% arg.names) stop("unconditional argument not meaningful here")
    if ("by.resids"%in% arg.names) stop("by.resids argument not supported")
    if ("all.terms"%in% arg.names) stop("all.terms argument not supported")
  }
  plot.gam(x,residuals=FALSE,rug=rug,se=se,pages=pages,select=select,scale=scale,
              n=n,n2=n2,pers=pers,theta=theta,phi=phi,jit=jit,xlab=xlab,
              ylab=ylab,main=main,ylim=ylim,xlim=xlim,too.far=too.far,
              all.terms=FALSE,shade=shade,shade.col=shade.col,
              shift=shift,trans=trans,seWithMean=seWithMean,
              unconditional=FALSE,by.resids=FALSE,
              scheme=scheme,...)
} ## plot.jam

predict.jam <- function(object,newdata,type="link",se.fit=FALSE,terms=NULL,
             block.size=NULL,newdata.guaranteed=FALSE,na.action=na.pass,...) {
  class(object) <- "gam" ## cheat!
  arg.names <- names(list(...))
  if (length(arg.names)>0) {
    if ("unconditional"%in% arg.names) warning("unconditional argument not meaningful here")
  }
  predict.gam(object,newdata,type=type,se.fit=se.fit,terms=terms,
             block.size=block.size,newdata.guaranteed=newdata.guaranteed,
             na.action=na.action,unconditional=FALSE,...)
} ## predict.jam

