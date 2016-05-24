#' Make a reference function in bayou
#'
#' This function generates a reference function from a mcmc chain for use in marginal likelihood
#' estimation.
#' 
#' @param chain An mcmc chain produced by \code{bayou.mcmc()} and loaded with \code{load.bayou()}
#' @param prior The prior function used to generate the mcmc chain
#' @param burnin The proportion of the mcmc chain to be discarded when generating the reference function
#' @param plot Logical indicating whether or not a plot should be created
#' 
#' @details Distributions are fit to each mcmc chain and the best-fitting distribution is chosen as
#' the reference distribution for that parameter using the method of Fan et al. (2011). For positive 
#' continuous parameters \code{alpha, sigma^2, halflife, Vy, w2, Ne}, Log-normal, exponential, gamma and weibull
#' distributions are fit. For continuous distributions \code{theta}, Normal, Cauchy and Logistic distributions 
#' are fit. For discrete distributions, \code{k}, negative binomial, poisson and geometric distributions are fit.
#' Best-fitting distributions are determined by AIC. 
#' 
#' @export
#' @return Returns a reference function of class "refFn" that takes a parameter list and returns the log density 
#' given the reference distribution. If \code{plot=TRUE}, a plot is produced showing the density of variable parameters
#' and the fitted distribution from the reference function (in red).
make.refFn <- function(chain, prior, burnin=0.3, plot=TRUE){
  model <- attributes(prior)$model
  contdists <- c("norm", "cauchy", "logis")
  poscontdists <- c("lnorm", "exp", "gamma", "weibull")
  discdists <- c("nbinom", "pois", "geom")
  bounddists <- c("beta")
  parorder <- switch(model,"QG"=c("dh2","dP","dw2","dNe","dsb","dk","dtheta","dloc"), "OU"=c("dalpha","dsig2","dsb","dk","dtheta","dloc"),"OUrepar"=c("dhalflife","dVy","dsb","dk","dtheta","dloc"))
  postburn <- round(burnin*length(chain[[1]]),0):length(chain[[1]])
  dists <- attributes(prior)$distributions
  dists <- dists[parorder]
  varDists <- dists[which(dists!="fixed")]
  refFx <- list()
  refNames <- list()
  dists <- list()
  parameters <- list()
  x <- NULL
  for(i in 1:length(varDists)){
    parname <- gsub('^[a-zA-Z]',"",names(varDists)[i])
    xx <- unlist(chain[[parname]][postburn])
    if(parname %in% c("P", "w2", "Ne", "alpha", "sig2", "halflife", "Vy")){
      tmpFits <- lapply(poscontdists, function(x) suppressWarnings(try(fitdist(xx, x), silent=TRUE)))
      tmpFits <- tmpFits[sapply(tmpFits, function(x) !(class(x)=="try-error"))]
      aic <- sapply(tmpFits, function(x) x$aic)
      fit <- tmpFits[[which(aic==min(aic,na.rm=TRUE))]]
      fitPars <- as.list(fit$estimate)
      fitPars$log <- TRUE
      fitName <- fit$distname
      fitfx <- get(paste("d",fitName, sep=""))
      refFx[[i]] <- .set.defaults(fitfx, defaults=fitPars)
    }
    if(parname %in% c("h2")){
      tmpFits <- lapply(bounddists, function(x) suppressWarnings(try(fitdist(xx, x), silent=TRUE)))
      tmpFits <- tmpFits[sapply(tmpFits, function(x) !(class(x)=="try-error"))]
      aic <- sapply(tmpFits, function(x) x$aic)
      fit <- tmpFits[[which(aic==min(aic,na.rm=TRUE))]]
      fitPars <- as.list(fit$estimate)
      fitPars$log <- TRUE
      fitName <- fit$distname
      fitfx <- get(paste("d",fitName, sep=""))
      refFx[[i]] <- .set.defaults(fitfx, defaults=fitPars)
    }
    if(parname %in% c("theta")){
      tmpFits <- lapply(contdists, function(x) suppressWarnings(try(fitdist(xx, x), silent=TRUE)))
      tmpFits <- tmpFits[sapply(tmpFits, function(x) !(class(x)=="try-error"))]
      aic <- sapply(tmpFits, function(x) x$aic)
      fit <- tmpFits[[which(aic==min(aic,na.rm=TRUE))]]
      fitPars <- as.list(fit$estimate)
      fitPars$log <- TRUE
      fitName <- fit$distname
      fitfx <- get(paste("d",fitName, sep=""))
      refFx[[i]] <- .set.defaults(fitfx, defaults=fitPars)
    }
    if(parname %in% c("k")){
      tmpFits <- lapply(discdists, function(x) fitdist(xx, x))
      aic <- sapply(tmpFits, function(x) x$aic)
      fit <- tmpFits[[which(aic==min(aic,na.rm=TRUE))]]
      fitPars <- as.list(fit$estimate)
      fitPars$log <- TRUE
      fitName <- fit$distname
      if(fitName=="nbinom") fitPars$prob=fitPars$size/(fitPars$size+fitPars$mu)
      fitfx <- get(paste("d",fitName, sep=""))
      refFx[[i]] <- .set.defaults(fitfx, defaults=fitPars)
    }
    if(parname %in% c("loc", "sb")){
      fitName <- paste("d", parname,sep="")
      refFx[[i]] <- attributes(prior)$functions[[fitName]]
    }
    dists[[i]] <- fitName
    parameters[[i]] <- fitPars
  }
  names(refFx) <- gsub('^[a-zA-Z]',"", names(varDists))
  par.names <- names(refFx)
  names(dists) <- par.names
  names(parameters) <- par.names
  if(plot){
    pars2plot <- par.names[(par.names %in% c("h2", "P", "w2", "Ne", "alpha", "halflife","sig2","Vy", "k", "theta"))]
    par(mfrow=c(ceiling(length(pars2plot)/2),2))
    for(i in 1:length(pars2plot)){
      plot(density(unlist(chain[[pars2plot[i]]][postburn])), main=pars2plot[i])
      if(pars2plot[i]=="k"){
        points(seq(ceiling(par('usr')[1]),floor(par('usr')[2]),1), refFx[[pars2plot[i]]](seq(ceiling(par('usr')[1]),floor(par('usr')[2]),1),log=FALSE),pch=21,bg="red")
      } else {x <- NULL; curve(refFx[[pars2plot[i]]](x,log=FALSE), add=TRUE, col="red")}
    }
  }
  refFUN <- function(pars,cache){
    if(any(!(par.names %in% names(pars)))) stop(paste("Missing parameters: ", paste(par.names[!(par.names %in% names(pars))],collapse=" ")))
    pars.o <- pars[match(par.names,names(pars))]
    pars.o <- pars.o[!is.na(names(pars.o))]
    densities <- sapply(1:length(pars.o),function(x) refFx[[x]](pars.o[[x]]))
    names(densities) <- par.names
    lnprior <- sum(unlist(densities,F,F))
    return(lnprior)
  }
  attributes(refFUN) <- list("model"=model,"parnames"=par.names,"distributions"=dists,"parameters"=parameters,"functions"=refFx)
  class(refFUN) <- c("refFn","function")
  return(refFUN)
}


#' Makes a power posterior function in bayou
#' 
#' This function generates a power posterior function for estimation of marginal likelihood using the stepping stone method
#' 
#' @param k The step in the sequence being estimated
#' @param Bk The sequence of steps to be taken from the reference function to the posterior
#' @param priorFn The prior function to be used in marginal likelihood estimation
#' @param refFn The reference function generated using \code{make.refFn()} from a preexisting mcmc chain
#' 
#' @details For use in stepping stone estimation of the marginal likelihood using the method of Fan et al. (2011).
#' @export
#' @return A function of class "powerposteriorFn" that returns a list of four values: \code{result} (the log density of the power posterior), 
#' \code{lik} (the log likelihood), \code{prior} (the log prior), \code{ref} the log reference density. 
make.powerposteriorFn <- function(k, Bk, priorFn, refFn){
  model <- attributes(priorFn)$model
  if(model != attributes(refFn)$model) stop("Error: prior and reference function are not of same type")
  powerposteriorFn <- function(k, Bk, pars, cache, dat, SE, model){
    lik <- bayou.lik(pars, cache, dat, model=model)$loglik
    prior <- priorFn(pars, cache)
    ref <- refFn(pars, cache)
    coeff <- c(Bk[k],Bk[k],(1-Bk[k]))
    result <- c(lik, prior, ref)
    result[coeff==0] <- 0
    result <- result*coeff
    result <- sum(result)
    return(list(result=result, lik=lik, prior=prior, ref=ref))
  }
  class(powerposteriorFn) <- c("powerposteriorFn", "function")
  return(powerposteriorFn)
}

.steppingstone.mcmc <- function(k, Bk, powerposteriorFn, tree, dat, SE=0, prior, ngen=10000, samp=10, chunk=100, control=NULL, tuning=NULL, new.dir=TRUE, plot.freq=500, outname="bayou", ticker.freq=1000, tuning.int=c(0.1,0.2,0.3), startpar=NULL, moves=NULL, control.weights=NULL){
  model <- attributes(prior)$model
  fixed <- gsub('^[a-zA-Z]',"",names(attributes(prior)$distributions)[which(attributes(prior)$distributions=="fixed")])
  if("loc" %in% fixed){
    fixed <- c(fixed,"slide")
  }
  if(is.null(moves)){
    moves <- switch(model,"QG"=list(h2=".multiplierProposal",P=".multiplierProposal",w2=".multiplierProposal",Ne=".multiplierProposal",k=".splitmerge",theta=".adjustTheta",slide=".slide"),
                    "OU"=list(alpha=".multiplierProposal",sig2=".multiplierProposal",k=".splitmerge",theta=".adjustTheta",slide=".slide"),
                    "OUrepar"=list(halflife=".multiplierProposal",Vy=".multiplierProposal",k=".splitmerge",theta=".adjustTheta",slide=".slide")) #,"OUcpp"=list(alpha=".multiplierProposal",sig2=".multiplierProposal",sig2jump=".multiplierProposal",k=".splitmergeSimmap",theta=".adjustTheta",slide=".slideCPP"), "QGcpp"=list(h2=".multiplierProposal",P=".multiplierProposal",w2=".multiplierProposal",Ne=".multiplierProposal",sig2jump=".multiplierProposal",k=".splitmergeSimmap",theta=".adjustTheta",slide=".slideCPP"),"OUreparcpp"=list(halflife=".multiplierProposal",Vy=".multiplierProposal",sig2jump=".multiplierProposal",k=".splitmergeSimmap",theta=".adjustTheta",slide=".slideCPP"))
    moves <- moves[which(!(names(moves) %in% fixed))]
  }
  
  cache <- .prepare.ou.univariate(tree,dat, SE=SE)
  dat <- cache$dat
  if(is.null(startpar)){
    if(any(fixed %in% c("h2", "P", "w2", "Ne", "halflife", "Vy"))){
      stop(paste("Parameters '", paste(fixed[fixed %in% c("h2", "P", "w2", "Ne", "halflife", "Vy")], collapse=" "), "' are set to be fixed but no starting values are supplied. 
                 Please specify starting parameter values",sep=""))
    }
    startpar <- priorSim(prior,cache$phy,model,nsim=1,plot=FALSE, exclude.branches=NULL)$pars[[1]]
    if(length(fixed)>0){
      assumed <- sapply(fixed, function(x) switch(x, "slide"="", "sb"="sb=numeric(0)", "k"= "k=0", "alpha"="alpha=0", "sig2"="sig2=0", "loc"="0.5*edge.length"))
      print(paste("Warning: Fixed parameters '", paste(fixed,collapse=", "), "' not specified, assuming values: ", paste(assumed,collapse=", "),sep="" ))
    }
    } 
  if(length(fixed)==0 & is.null(control.weights)){
    ct <- .buildControl(startpar, prior)
  } else {
    if(is.null(control.weights)){
      control.weights <- switch(model,"OU"=list("alpha"=4,"sig2"=2,"theta"=4,"slide"=2,"k"=10),"QG"=list("h2"=5,"P"=2,"w2"=5,"Ne"=5,"theta"=5,"slide"=3,"k"=20),"OUrepar"=list("halflife"=5,"Vy"=3,"theta"=5,"slide"=3,"k"=20))
      #"OUcpp"=list("alpha"=3,"sig2"=3,"sig2jump"=3,"theta"=3,"slide"=5,"k"=10),"QGcpp"=list("h2"=1,"P"=1,"w2"=2,"Ne"=2,"sig2jump"=3,"theta"=3,"slide"=5,"k"=10),"OUreparcpp"=list("halflife"=3,"Vy"=3,"sig2jump"=3,"theta"=3,"slide"=5,"k"=10)
      control.weights[fixed[fixed %in% names(control.weights)]] <- 0
    } else {control.weights <- control.weights}
    ct <- .buildControl(startpar, prior, move.weights=control.weights)
  }
  
  if(is.null(tuning)){
    D <- switch(model, "OU"=list(alpha=1, sig2= 1, k = 4,theta=2,slide=1), "QG"=list(h2=1, P=1, w2=1, Ne=1, k = 4, theta=2, slide=1), "OUrepar"=list(halflife=1, Vy=1, k=4, theta=2, slide=1))#,"OUcpp"=list(alpha=1, sig2= 1,sig2jump=2, k = 4,theta=2,slide=1),"QGcpp"=list(h2=1,P=1,w2=1,Ne=1,sig2jump=2,k=4,theta=2,slide=1),"OUreparcpp"=list(halflife=1,Vy=1,sig2jump=2,k=4,theta=2,slide=1))
  } else {D <- tuning}
  
  if(is.logical(new.dir)){
    if(new.dir){
      dir.name <- paste(sample(LETTERS,10),collapse="")
      dir <- paste(tempdir(),"/",dir.name,"/",sep="")
      dir.create(dir)
    } else {
      dir <- paste(getwd(),"/",sep="")
      dir.name <- sapply(strsplit(dir,'/'), function(x) x[length(x)])
      }
    } else {
      dir.name <- paste(sample(LETTERS,10),collapse="")
      dir <- paste(new.dir,"/",dir.name,"/",sep="")
      dir.create(dir)
    }
  
  #mapsb <<- file(paste(dir, outname,".", k, ".sb",sep=""),open="w")
  #mapsloc <<- file(paste(dir, outname,".", k, ".loc",sep=""),open="w")
  #mapst2 <<- file(paste(dir, outname,".", k, ".t2",sep=""),open="w")
  #pars.output <<- file(paste(dir, outname,".", k, ".pars",sep=""),open="w")
  files <- list(mapsb=file(paste(dir, outname,".",k,".sb",sep=""),open="a"), 
                mapsloc=file(paste(dir, outname,".",k,".loc",sep=""),open="a"),
                mapst2=file(paste(dir, outname,".",k,".t2",sep=""),open="a"),
                pars.output=file(paste(dir, outname,".",k,".pars",sep=""),open="a"))
  
  
  oldpar <- startpar
  store <- list("out"=list(), "sb"=list(), "loc"=list(), "t2"=list())
  
  
  lik.fn <- bayou.lik
  oll  <- lik.fn(oldpar, cache, dat, model=model)$loglik
  pr1 <- prior(oldpar,cache)
  parorder <- switch(model,"QG"=c("h2","P","w2","Ne","k","ntheta","theta"), "OU"=c("alpha","sig2","k","ntheta","theta"),"OUrepar"=c("halflife","Vy","k","ntheta","theta"),"OUcpp"=c("alpha","sig2","sig2jump","k","ntheta","theta"))#,"QGcpp"=c("h2","P","w2","Ne","sig2jump","k","ntheta","theta"),"OUreparcpp"=c("halflife","Vy","sig2jump","k","ntheta","theta"))
  
  accept.type <- NULL
  accept <- NULL
  if(!is.null(plot.freq)){
    tr <- .toSimmap(.pars2map(oldpar, cache),cache)
    tcols <- makeTransparent(rainbow(oldpar$ntheta),alpha=100)
    names(tcols)<- 1:oldpar$ntheta
    phenogram(tr,dat,colors=tcols,ftype="off", spread.labels=FALSE)
    plot.dim <- list(par('usr')[1:2],par('usr')[3:4])
  }
  #tuning.int <- round(tuning.int*ngen,0)

  pB.old <- powerposteriorFn(k, Bk, oldpar, cache, dat, SE=0, model=model)
  Ref <- NULL
  
  for (i in 1:ngen){
    ct <- .updateControl(ct, oldpar, fixed)
    u <- runif(1)
    prop <- .proposalFn(u,ct,D,moves,cache,oldpar)
    new.pars <- prop$pars
    #new.cache <- prop$prop$cache
    accept.type <- c(accept.type,paste(prop$decision,prop$move,sep="."))
    pr2 <- prior(new.pars,cache)
    hr <- prop$hr
    pB.new <- powerposteriorFn(k, Bk, new.pars, cache, dat, SE=SE, model=model)
    bhr <- ifelse(!is.finite(hr), hr, Bk[k]*hr)
    if (runif(1) < exp(pB.new$result-pB.old$result+bhr)){
      oldpar <- new.pars
      oldpar <- new.pars
      pr1 <- pB.new$prior
      oll <- pB.new$lik
      pB.old <- pB.new
      accept <- c(accept,1)
      if(i %% samp ==0){
        Ref <- c(Ref,pB.new$ref)
      }
    } else {
      accept <- c(accept,0)
      if(i %% samp ==0){
        Ref <- c(Ref,pB.new$ref)
      }
    }
    store <- .store.bayou(i, oldpar, oll, pr1, store, samp, chunk, parorder, files)
    if(!is.null(plot.freq)){
      if(i %% plot.freq==0){
        tr <- .toSimmap(.pars2map(oldpar, cache),cache)
        tcols <- makeTransparent(rainbow(oldpar$ntheta),alpha=100)
        names(tcols)<- 1:oldpar$ntheta
        plot(plot.dim[[1]],plot.dim[[2]],type="n",xlab="time",ylab="phenotype")
        mtext(paste("gens = ",i," lnL = ",round(oll,2)),3)
        #regime.plot probably doesn't work for simmaps
        #try(regime.plot(oldpar,tr$tree,tcols,type="density",model=model),silent=TRUE)
        phenogram(tr,dat,colors=tcols,ftype="off",add=TRUE, spread.labels=FALSE)
      }
    }
    #if(i %in% tuning.int){
    #  D <- tune.D(D,accept,accept.type)$D
    #}
    if(i%%ticker.freq==0){
      alpha <- switch(model,"QG"=QG.alpha(oldpar),"OU"=oldpar$alpha,"OUrepar"=OU.repar(oldpar)$alpha,"OUcpp"=oldpar$alpha,"QGcpp"=QG.alpha(oldpar),"OUrepar"=OU.repar(oldpar)$alpha)
      sig2 <- switch(model,"QG"=QG.sig2(oldpar),"OU"=oldpar$sig2,"OUrepar"=OU.repar(oldpar)$sig2,"OUcpp"=oldpar$sig2,"QGcpp"=QG.sig2(oldpar),"OUrepar"=OU.repar(oldpar)$sig2)
      tick <- c(round(Bk[k],2), i,oll,pr1,pB.old$ref,log(2)/alpha,sig2/(2*alpha),oldpar$k,tapply(accept,accept.type,mean))
      tick[-1] <- round(tick[-1],2)
      names(tick)[1:8] <- c('Bk_k', 'gen','lnL','prior','ref','half.life','Vy','K')
      if(i==ticker.freq){
        cat(c(names(tick),'\n'),sep='\t\t\t')
      }
      cat(c(tick,'\n'),sep='\t\t\t')
    }
  }
  lapply(files, close)
  out <- list('model'=model, 'dir.name'=dir.name,'dir'=dir, 'outname'=paste(outname,".",k,sep=""), 'accept'=accept,'accept.type'=accept.type, 'ref'=Ref)
  class(out) <- c("ssbayouFit", "list")
  return(out)
}


#' Stepping stone estimation of the marginal likelihood for a bayou model
#' 
#' Estimates the marginal likelihood of a bayou model by generating mcmc chains for power posteriors for a series of steps from 0 to 1, 
#' progressing from a reference distribution to the posterior distribution.
#' 
#' @param Bk A vector sequence from 0 to 1 that gives the exponents of the power posterior distribution to take from the reference distribution 
#' to the reference distribution. (See details)
#' @param chain A mcmc chain used to generate the reference distribution (see \code{make.refFn} for details).
#' @param tree A phylogenetic tree of class "phylo"
#' @param dat A named vector of continuous trait data
#' @param SE A vector giving the standard error of the trait data. If a single number is given, standard errors are assumed to be constant across the phylogeny.
#' @param prior The prior function used to generate the mcmc chain.
#' @param startpar The starting parameter values to be used. If any parameters are set as "fixed", this should be specified. If \code{NULL}, then the parameters are 
#' drawn from the prior distribution.
#' @param burnin The initial proportion of the provided mcmc chain to be discarded when generating the reference function
#' @param ngen The number of mcmc generations to be run for each step of the stepping stone alogrithm. 
#' @param powerposteriorFn The power posterior function to be used. If \code{NULL}, this is generated from the provided mcmc chain.
#' @param parallel A logical indicating whether or not the chains should be run in parallel. 
#' @param ... Other parameters passed to the mcmc algorithm, see \code{bayou.mcmc()}.
#' 
#' @details This function estimates the marginal likelihood of a bayou model by using stepping stone estimation from a reference distribution to the posterior 
#' distribution using the method of Fan et al. (2011). The vector \code{Bk} provides a sequence from 0 to 1. The length of this sequence determines the number
#' of mcmc chains that will be run, and the values are used as the exponents of the power posterior function, stepping from purely the reference distribution (k=0)
#' to purely the posterior distribution (k=1). These chains can be run in parallel if \code{parallel} is set to \code{TRUE}. The number of cores available is determined
#' by a call to \code{detectCores}, and can be set by .... Note that when run in parallel, progress within each
#' of the individual mcmc chains will not be reported, and if \code{ngen} is high, it may take a considerable amount of time to run. Furthermore, if many samples 
#' are saved from each mcmc run, and a number of steps along \code{Bk} is large, the returned object may require a substantial amount of memory. 
#' 
#' @export
#' @return A list of class "ssMCMC" that provides the log marginal likelihood \code{lnr}, a list of the individual normalizing constants estimated at each step \code{lnrk},
#' a list of the mcmc chains used for importance sampling to estimating the marginal likelihood at each step \code{chains}, and mcmc fit data from each of the runs \code{fits}.
#' Note that this object may become quite large if a number of chains are run for many generations. To reduce the number of samples taken, increase the parameter \code{samp} (default = 10)
#' which sets the frequency at which samples are saved in the mcmc chain. 
steppingstone <- function(Bk, chain, tree, dat, SE=0, prior, startpar=NULL, 
                          burnin=0.3, ngen=10000, powerposteriorFn=NULL, 
                          parallel=FALSE, ...){
    model <- attributes(prior)$model
    if(is.null(powerposteriorFn)){
      cat("Making power posterior function from provided mcmc chain...\n")
      ref <- suppressWarnings(make.refFn(chain, prior, plot=TRUE))
      ppost <- make.powerposteriorFn(1, Bk=seq(0,1,length.out=10), prior, ref)
    } else {ppost <- powerposteriorFn}
    cat("Running mcmc chains...\n")
    if(parallel==TRUE){
      k <- NULL; i <- NULL
      ssfits <- foreach(k = 1:length(Bk)) %dopar% {
        .steppingstone.mcmc(k=k, Bk=Bk, tree=tree, dat=dat, SE=SE, prior=prior, powerposteriorFn=ppost, startpar=startpar, plot.freq=NULL, ngen=ngen,  ...)
      }
    } else {
      k <- NULL; i <- NULL; ssfits <- list()
      for(k in 1:length(Bk)){
        ssfits[[k]] <- .steppingstone.mcmc(k=k, Bk=Bk, tree=tree, dat=dat, SE=SE, prior=prior, powerposteriorFn=ppost, startpar=startpar, plot.freq=NULL, ngen=ngen,  ...)
      }
    }
  cat("Loading mcmc chains...\n")
  if(parallel){
    Kchains <- foreach(i = 1:length(Bk)) %dopar% {
      load.bayou(ssfits[[i]], save.Rdata=FALSE, cleanup=FALSE)
    }
  } else {
    Kchains <- list()
    for (i in 1:length(Bk)){
      Kchains[[i]] <- load.bayou(ssfits[[i]], save.Rdata=FALSE, cleanup=FALSE)
    }
  }
  Kchains <- lapply(1:length(Kchains), function(x){Kchains[[x]]$ref <- ssfits[[x]]$ref; Kchains[[x]]})
  postburn <- round(burnin*length(Kchains[[1]][[1]]),0):length(Kchains[[1]][[1]])
  lnr <- .computelnr(Kchains, ssfits, Bk, postburn)
  out <- list(lnr= lnr$lnr, lnrk = lnr$lnrk, Bk=Bk, chains=Kchains, fits=ssfits)
  class(out) <- c("ssMCMC", "list")
  return(out)
}

#' S3 method for printing ssMCMC objects
#' 
#' @param x An ssMCMC object
#' @param ... Optional arguments passed to print
#' 
#' @export
#' @method print ssMCMC
print.ssMCMC <- function(x, ...){
  cat("Stepping stone estimation of marginal likelihood\n")
  cat("Marginal Likelihood:\n")
  print(x$lnr, ...)
  cat(paste("A total of ", length(x$Bk), " power posteriors were run along the sequence: ",paste(round(x$Bk,5), collapse="\t\t"), "\n", sep=""))
  cat("lnr_k", round(unlist(x$lnrk),2))  
}
#' S3 method for plotting ssMCMC objects
#' 
#' @param x An 'ssMCMC' object
#' @param ... Additional arguments passed to \code{plot}
#' 
#' @details Produces 4 plots. The first 3 plot the prior, reference function and likelihood. Different colors
#' indicate different power posteriors for each. These chains should appear to be well mixed. The final plot
#' shows the sum of the marginal likelihood across each of the steps in the stepping stone algorithm. 
#' 
#' @export
#' @method plot ssMCMC
plot.ssMCMC <- function(x, ...){
  par(mfrow=c(2,2))
  if(is.null(attributes(x)$burnin)){
    start <- 1
  } else {
    start <- round(attributes(x)$burnin*length(x$chains[[1]][[1]]),0)
  }
  postburn <- start:length(x$chains[[1]][[1]])
  lnL <- lapply(x$chains, function(x) x$lnL[postburn])
  rangelnL <- c(min(unlist(lnL))-2, max(unlist(lnL))+2)
  plot(0,0,type="n", xlim=c(0,length(unlist(lnL))), ylim=rangelnL,xaxt="n",xlab="",ylab="lnL", main="lnL",...)
  xindex <- lapply(1:length(lnL), function(x) (x-1)*length(lnL[[1]]) + 1:length(lnL[[1]]))
  sapply(1:length(lnL), function(x) lines(xindex[[x]], lnL[[x]], col=x))
  abline(v=seq(0,length(unlist(lnL)), length.out=length(lnL)+1),lty=2)
  
  pr <- lapply(x$chains, function(x) x$prior[postburn])
  rangepr <- c(min(unlist(pr))-2, max(unlist(pr))+2)
  plot(0,0,type="n", xlim=c(0,length(unlist(pr))), ylim=rangepr,xaxt="n",xlab="",ylab="Ln prior", main="ln prior",...)
  xindex <- lapply(1:length(pr), function(x) (x-1)*length(pr[[1]]) + 1:length(pr[[1]]))
  sapply(1:length(pr), function(x) lines(xindex[[x]], pr[[x]], col=x))
  abline(v=seq(0,length(unlist(pr)), length.out=length(pr)+1),lty=2)
  
  ref <- lapply(x$chains, function(x) x$ref[postburn])
  rangeref <- c(min(unlist(ref))-2, max(unlist(ref))+2)
  plot(0,0,type="n", xlim=c(0,length(unlist(ref))), ylim=rangeref,xaxt="n",xlab="",ylab="Ln ref", main="ln ref",...)
  xindex <- lapply(1:length(ref), function(x) (x-1)*length(ref[[1]]) + 1:length(ref[[1]]))
  sapply(1:length(ref), function(x) lines(xindex[[x]], ref[[x]], col=x))
  abline(v=seq(0,length(unlist(ref)), length.out=length(ref)+1),lty=2)
  
  plot(x$Bk, c(0, cumsum(x$lnrk)), ylab="ln r", xlab="power posterior",pch=21, bg=1:length(ref),cex=1.5, ...)
  lines(x$Bk, c(0,cumsum(x$lnrk)))
}

.pull.rsample <- function(samp, chain, fit, refFn){
  #pars.list <- lapply(samp,function(y) pull.pars(y,chain,model=model))
  #emap.list <- lapply(samp,function(y) read.emap(chain$branch.shift[[y]],chain$location[[y]],chain$t2[[y]],cache$phy)$emap)
  L <- chain$lnL[samp]+chain$prior[samp]-fit$ref[samp]
  Lmax <- max(L)
  Lfactored <- L-Lmax
  return(list(Lmax=Lmax,Lfactored=Lfactored))
}

#' Compute marginal likelihood
#' 
#' \code{computelnr} computes the marginal likelihood of a set of chains estimated via stepping stone
#' sampling and produced by the function \code{steppingstone}
.computelnr <- function(Kchains,ssfits,Bk,samp){
  lnr <- list()
  for(i in 1:(length(Bk)-1)){
    Lk <- .pull.rsample(samp, Kchains[[i]],ssfits[[i]])
    lnr[[i]] <- (Bk[i+1]-Bk[i])*Lk$Lmax+log(1/length(Lk$Lfactored)*sum(exp(Lk$Lfactored)^(Bk[i+1]-Bk[i])))
  }
  return(list("lnr"=sum(unlist(lnr)),"lnrk"=lnr))
} 
