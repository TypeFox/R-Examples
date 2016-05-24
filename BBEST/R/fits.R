#####################################################################################
#
#
#   FUNCTIONS TO PERFORM BACKGROUND FIT VIA DIFFERENTIAL EVOLUTION ALGORITHM
#
#


#####################################################################################
# set.control(CR, F, NP, itermax, parallelType)
#
# returns:         wrapper to set.DEoptim.control  
# arguments
#   CR:            crossover probability from interval [0,1]
#   F:             differential weighting factor from interval [0,2]
#   NP:            number of population members
#   itermax:       the maximum iteration (population generation) allowed
#   parallelType:  defines the type of parallelization to employ. 
set.control <- function(CR=.85, F=.7, NP=300, itermax=2000, parallelType=1){
  control <- list()
  control$CR <- CR
  control$F <- F
  control$NP <- NP	
  control$itermax <- itermax
  control$parallelType <- parallelType
	  
  return(control)  
}
	
	

#####################################################################################
# set.DEoptim.control(control, knots.n)
#
# returns:         DEoptim control parameters  
# arguments
#   control:       return value of set.control
#   knots.n:       knots number
set.DEoptim.control <- function(control=list(), knots.n){
  if(is.null(control$CR))
    CR <- 0.85  
  else
    CR <- control$CR
  if(is.null(control$F))
    F <- 0.7  
  else
    F <- control$F
  if(is.null(control$NP))
    NP <- round(2.5*knots.n*10)  
  else    
    NP <- control$NP
  if(is.null(control$itermax))
    itermax <- 2000 
  else  
  	itermax <- control$itermax
  if(is.null(control$parallelType))
    parallelType <- 1 
  else  
	parallelType <- control$parallelType
	
  DE.control <- DEoptim.control(CR=CR, F=F, NP=NP, itermax=itermax, parallelType=parallelType, 
#	packages = list("PerformanceAnalytics"), 
	parVar=list("basisMatrix", "basisSpline", 
		  "logPosterior", "logPriorF", "DMatrix", "logLikelihoodBkg", "logLikelihoodSignal",
		  "get.bkg", "logProbabilityBkgR", "logLikelihoodGrGauss", "logLikelihoodGrCorr", 
		  "logPriorBkgRSmooth", "logPriorBkgRGP", "Dx", "invert.order", "set.SB", 
		  "logPosteriorAnalyt", "bkg.analyt", "logLikelihoodBkgAnalyt", "logLikelihoodSignalAnalyt"))
					  
  DE.control				  
}


#####################################################################################
# do.fit(data, bounds.lower, bounds.upper, scale, knots.x, knots.n, stdev, control, p.bkg, save.to)
#
# returns:          Performs evolutionary global background optimization via DEoptim  
# arguments
#   data:           object of type data. Contatins experimental data and fit parameters 
#   bounds.lower,   lower and upeer boundaries for background 
#   bounds.upper:    
#   scale:          scale factor  
#   knots.x:        the spline knot x-values. 
#   knots.n:        number of knots. 
#   stdev:          Whether to calculate uncertainty in background estimation 
#   control:        he return value of set.control
#   p.bkg:          the probability that a single pixel contains "only" background 
#   save.to:        the name of the file where the results will be saved 
do.fit <- function(data, bounds.lower, bounds.upper, scale=c(1,1), knots.x=NA, 
           knots.n=NA, analytical=FALSE, stdev=TRUE, control=list(), p.bkg=.5, save.to=""){

# 1. Prepare data to fit...	
  cat("Preparing data... \n")
  if( (scale[1]==1) && (scale[2]==1) ) scale <- NA
    cov.r <- ADP <- knots.y <- pars <- NA
	alpha <- 1
	if(is.null(data$Gr))
	  Gr <- NA
    else{
	  Gr <- data$Gr
	  data$Gr <- NA
    }
	if(is.null(data$SB))
	  data$SB <- rep(0, length(data$x))
	else if (is.na(data$SB[1]))
	  data$SB <- rep(0, length(data$x))

# 2. prepare fit params...
  if(is.na(knots.x[1]) && !is.na(knots.n))  
    knots.x <- seq(min(data$x), max(data$x), length=knots.n)
	if(is.na(knots.n))
	  knots.n <- length(knots.x)
    
	if(analytical==TRUE){
    iter.0 <- as.relistable(list(pars=rep(0, 6)))
    bounds.lower <- c(pars=c(-abs(bounds.lower[1]*10), -5, -5, -20, 0, 0))
    bounds.upper <- c(pars=c(abs(bounds.upper[1]*10), 10, 5, 2000, 
	                      max(data$x), max(data$x)*2))
    knots.n <- 6
  }
	else{
    iter.0 <- as.relistable(list(knots.y=rep(0, knots.n)))
    bounds.lower <- c(knots.y=rep(bounds.lower[1], knots.n))
    bounds.upper <- c(knots.y=rep(bounds.upper[1], knots.n))
  }
	
	if(!is.na(scale[1])){
	  iter.0$alpha <- 0
    bounds.lower <- c(bounds.lower, alpha=scale[1])
    bounds.upper <- c(bounds.upper, alpha=scale[2])	  
	}
	if(!is.null(data$fitADP)){
	  if(data$fitADP$oneADP==TRUE)
	    nn <- 1
	  else
	    nn <- length(data$fitADP$n.atoms)
	    iter.0$ADP <- rep(0, nn)
      bounds.lower <- c(bounds.lower, ADP=rep(data$fitADP$ADP.lim[1], nn))
      bounds.upper <- c(bounds.upper, ADP=rep(data$fitADP$ADP.lim[2], nn))	  
	}			
  DE.control <- set.DEoptim.control(control, knots.n)
 
# 3. Starting fit! 
#
# initial guess
#  cc <- matrix(nrow=DE.control$NP, ncol=length(knots.x))
#  Phi <- basisMatrix(data$x, knots.x)
#  for(ii in 1:DE.control$NP){
#    spar=(0.9-0.5)/DE.control$NP*ii+0.5
#    lowpass.spline <- smooth.spline(data$x,data$y-data$SB, spar = spar) ## Control spar for amount of smoothing
#    cc[ii,] <- c(solve(t(Phi)%*%Phi) %*% t(Phi) %*% (predict(lowpass.spline, data$x)$y) )
#  }
#  DE.control$initialpop=cc
  cat("Starting DifEv algorithm... \n")
  doPbkgIter <- FALSE
  if(p.bkg==-1){
    doPbkgIter <- TRUE
    p.bkg <- 0.02
  }    
  DEoptim.fit <- DEoptim(get.posterior, lower = bounds.lower, upper = bounds.upper, 
                         control= DE.control, skel=iter.0, data=data, 
                         knots.x=knots.x, Gr=Gr, p.bkg=p.bkg)
  if(doPbkgIter){
    cat("\n\n P.bkg iteration... \n\n")
    if(analytical==TRUE){
      pars <- DEoptim.fit$optim$bestmem[1:6]
      bkg <- bkg.analyt(pars=pars, x=data$x)
    } 
    else{	
      knots.y <- DEoptim.fit$optim$bestmem[1:knots.n]
      bkg <- get.bkg(x=data$x, knots.x=knots.x, knots.y=knots.y)  
    }
    dev.norm <- (data$y-bkg-data$SB)/data$sigma
    p.bkg <- 0.5*50^((1-dev.norm^2)/8)
    
    DEoptim.fit <- DEoptim(get.posterior, lower = bounds.lower, upper = bounds.upper, 
                         control= DE.control, skel=iter.0, data=data, 
                         knots.x=knots.x, Gr=Gr, p.bkg=p.bkg)       
  }
  
	if(analytical==TRUE){
    pars <- DEoptim.fit$optim$bestmem[1:6]
    bkg <- bkg.analyt(pars=pars, x=data$x)
  } 
  else{	
    knots.y <- DEoptim.fit$optim$bestmem[1:knots.n]
    bkg <- get.bkg(x=data$x, knots.x=knots.x, knots.y=knots.y)  
	}
	if(!is.null(data$fitADP)){
    ADP <- tail(DEoptim.fit$optim$bestmem, nn)
	  data <- set.SB(data, n.atoms=data$fitADP$n.atoms, 
	                   scatter.length=data$fitADP$scatter.length, ADP=ADP)
  }
	data$Gr <- Gr
	if(!is.na(scale[1])){ # recalculating signal and params...
	  alpha <- DEoptim.fit$optim$bestmem[knots.n+1]
	  data$y <- (data$y-bkg-data$SB)*alpha + data$SB + bkg
	  data$lambda <- data$lambda*alpha
	  data$sigma <- data$sigma*alpha
	  if(!is.na(Gr[1]))
        data$Gr$sigma.r <- Gr$sigma.r * alpha
	}  
	
  cat("Background estimation complete! \n")
  
  if(stdev==TRUE && analytical==FALSE){
    cat("Calculating uncertainty in background... \n")
    stdev <- get.hess(data, knots.x, knots.y, Gr=data$Gr, r=seq(0, 2, 0.01), p.bkg=p.bkg)
  }
  
  curves <- list(y=data$y, bkg=bkg, SB=data$SB)
  knots <- list(x=knots.x, y=knots.y)
  fit.details <- list(lambda=data$lambda, sigma=data$sigma, knots.n=knots.n, 
                      control=control, Gr=data$Gr, n.atoms=data$fitADP$n.atoms, 
                      scatter.length=data$fitADP$scatter.length, id=data$id,
                      bounds.lower=bounds.lower['knots.y1'],  
                      bounds.upper=bounds.upper['knots.y1'])
  fit.results <- list(x=data$x, curves=curves, uncrt=stdev, knots=knots, 
                    pars=pars, scale=alpha, ADP=ADP, fit.details=fit.details)
  
  if(save.to!=""){
      cat("Saving results to file ", save.to, "\n")
      save(fit.results, file=save.to)
  }
  cat("...done! \n")
  return(fit.results)
}


#####################################################################################
# do.fit(data, bounds.lower, bounds.upper, knots.n.left, knots.n.right, x.boundary, control, save.to)
#
# returns:          Performs evolutionary global background optimization via DEoptim 
#                   (wrapper to do.fit) for several banks. 
# arguments
#   data:           object of type data. Contatins experimental data and fit parameters 
#   bounds.lower,   lower and upeer boundaries for background 
#   bounds.upper:    
#   knots.n.left,   specify knot positions. knots.n.left and knots.n.right knots are created on the 
#   knots.n.right,  left and on the right of x.boundary point, respectively       
#   x.boundary:            
#   control:        he return value of set.control
#   save.to:        the name of the file where the results will be saved 
do.fit.banks <- function(data, bounds.lower, bounds.upper, knots.n.left=NA,
       knots.n.right=NA, x.boundary=NA, analytical=FALSE, control, save.to=""){
  N <- length(data)
  fit.res <- list()
  knots.x <- NA
  for(i in 1:N){
    cat("\n\n ===================================\n\n")
    
    cat("Fitting bank # ",i,"\n\n")
    x <- data[[i]]$x
    if(analytical==FALSE){
        dx <- ((max(x)-x.boundary)/(knots.n.right-1) +  x.boundary/(knots.n.left-1) )/2
        knots.x <- seq(0, x.boundary, length=knots.n.left)
        knots.x <- c(knots.x, seq(x.boundary+dx, max(x), length=knots.n.right))
    }
    fit.res[[i]] <- do.fit(data[[i]], bounds.lower, bounds.upper, knots.x=knots.x,
        analytical=analytical, stdev=FALSE, control=control, save.to="")
    fit.res[[i]]$fit.details$id <- data[[i]]$id     
  }
  
  if(save.to!=""){
    cat("Saving results to file ", save.to, "\n")
    save(fit.res, file=save.to)
  }
  fit.res
}



#####################################################################################
# do.fit(data, bounds.lower, bounds.upper, knots.n.left, knots.n.right, x.boundary, control, save.to)
#  
do.iter <- function(fit.results, local=TRUE, eps=1e-4, n.iter=10000, save.to=""){

  dat <- list(x=fit.results$x, y=fit.results$curves$y, SB=fit.results$curves$SB, 
              lambda=fit.results$fit.details$lambda, sigma=fit.results$fit.details$sigma)
  knots.x <- fit.results$knots$x
  knots.y <- fit.results$knots$y
  
  
  cat("Adjusting baseline...", save.to, "\n")
  if(local)
    cc <- grad.descent(data=dat, knots.x, knots.y, Gr=NA, p.bkg=0.5, eps=eps, N=n.iter)
  else{
    control <- fit.results$fit.details$control
    bounds.lower <- fit.results$fit.details$bounds.lower
    bounds.upper <- fit.results$fit.details$bounds.upper
    cat("\n\n Starting DifEv to find bkg with no low-r constraints... \n\n") 
    ff <- do.fit(dat, bounds.lower, bounds.upper, knots.x=knots.x, 
                 stdev=FALSE, control=control)
    cc <- ff$knots$y                
  }
 
  if(any(is.na(cc))) stop("perhaps adjust parameters")
  
  l <- fit.results$curves$bkg
  l.no.r <- get.bkg(dat$x, knots.x, cc) 

  r <- seq(0, 1.0, 0.005)
  gr <- sineFT(f.Q=dat$y-l-1, Q=dat$x, r=r)
  gr.no.r <- sineFT(f.Q=dat$y-l.no.r-1, Q=dat$x, r=r)

  d <- gr.no.r- gr   # < 0!!! Difference between fits with and without Gr info

  l.corr <- 0
  dr <- r[2]-r[1]
  for(i in 1:length(dat$x)){
    l.corr[i] <- sum(d*sin(dat$x[i]*r)*dr/dat$x[i])
  }

  dat$SB <- dat$SB - l.corr
 
  cat("\n\n Estimating background for the corrected data... \n\n") 
  if(local)
    cc2 <- grad.descent(data=dat, knots.x, cc, Gr=fit.results$fit.details$Gr, p.bkg=0.5, eps=eps, N=n.iter) 
  else{
    Gr <- fit.results$fit.details$Gr
    dat <- set.Gr(dat, r1=Gr$r1, r2=Gr$r2, rho.0=Gr$rho.0, 
                   type1=Gr$type1, type2=Gr$type2)
    ff2 <- do.fit(dat, bounds.lower, bounds.upper, knots.x=knots.x, 
                  stdev=FALSE, control=control)
    cc2 <- ff2$knots$y                
  }
  
  if(any(is.na(cc2))) stop("perhaps adjust parameters")
  
  l2 <- get.bkg(dat$x, knots.x, cc2) 
  
  fit.results$curves$bkg <- l2
  fit.results$knots$y <- cc2
  fit.results$curves$corr <- l.corr
  
  if(length(fit.results$uncrt)>1){
    cat("Calculating uncertainty in background... \n")
    stdev <- get.hess(dat, knots.x, knots.y, Gr=fit.results$fit.details$Gr, r=seq(0, 2, 0.01), p.bkg=.5)
  }
  if(save.to!=""){
    cat("Saving results to file ", save.to, "\n")
    save(fit.results, file=save.to)
  }
  return(fit.results)
}

# MATMUL
#matmul <- function(A, B){
#  if(!is.loaded("matmul"))
#    dyn.load("./matmul.dll")
	
#  matrix(.C("matmul_matmul_R",
#            heightA=as.integer(nrow(A)), 
#			widthA=as.integer(ncol(A)), 
#			widthB=as.integer(ncol(B)), 
#			A=as.double(c(A)), 
#			rstrideA=as.integer(ncol(A)), 
#			B=as.double(c(B)), 
#			rstrideB=as.integer(ncol(B)),
#			C=as.double(rep(0, nrow(A)*ncol(B))), 
#			rstrideC=as.integer(ncol(B)))$C,
#		  nrow=nrow(A), ncol=ncol(B))
#}