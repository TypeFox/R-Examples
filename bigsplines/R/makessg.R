makessg <-
  function(formula,family,data,type=NULL,nknots=NULL,rparm=NA,
           lambdas=NULL,skip.iter=TRUE,se.lp=FALSE,rseed=1234,
           gcvopts=NULL,knotcheck=TRUE,gammas=NULL,weights=NULL,
           gcvtype=c("acv","gacv","gacv.old")){
    ###### Makes Generalized Smoothing Spline Anova models
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: March 10, 2015
    
    ### get initial info 
    mf <- match.call()
    m <- match(c("formula","data"),names(mf),0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())          # mf contains response and predictor info
    mt <- attr(mf, "terms")                 # mt contains model info and terms 
    et <- attr(mt,"factors")                # et is effects table
    mfdim <- dim(et)                        # dim of effects table
    ndpts <- dim(mf)[1]                     # total number of data points
    nxvar <- mfdim[1]-1L                    # number of predictors
    nterm <- mfdim[2]                       # number of model terms
    xynames <- row.names(et)
    xnames <- xynames[2:(nxvar+1L)]
    tnames <- colnames(et)
    if(any(colSums(et>0L)>3L)){stop("Four-way (and higher-order) interactions are not supported.")}
    
    ### check gcvtype
    gcvtype <- gcvtype[1]
    if(!any(gcvtype==c("acv","gacv","gacv.old"))){stop("Input 'gcvtype' must be one of three specific options.")}
    
    ### check response
    yvar <- as.matrix(model.response(mf, "numeric")+0.0)   # response variable
    if(ncol(yvar)>1){stop("Response must be unidimensional (vector).")}
    if(nrow(yvar)!=ndpts){stop("Response vector must have same length as predictors.")}
    ysm <- sum(yvar)
    nunewr <- ndpts
    fweights <- rep.int(1L,ndpts)
    yty <- yvar^2
    
    ### check weights
    if(is.null(weights[1])==FALSE){
      weights <- as.numeric(weights)
      if(length(weights)!=ndpts){stop("Weights vector must have same length as response.")}
      if(any(weights<=0)){stop("Weights must be positive.")}
      if(is.na(rparm[1])==FALSE & min(weights)<=0.005){stop("Minimum weight is too small (less than 0.005).")}
      wtchk <- TRUE
    } else {
      weights <- 1
      wtchk <- FALSE
    }
    
    ### check family
    if(length(family)==1L){
      if(class(family[[1]])!="character"){stop("Input 'family' should be character class.")}
      family <- family[[1]]
      dispersion <- NULL
    } else if(length(family)==2L){
      dispersion <- as.numeric(family[[2]])
      if(class(family[[1]])!="character"){stop("Input 'family' should be character class.")}
      family <- family[[1]]
    } else{stop("Input 'family' should be 1 or 2 element vector or list.")}
    if(family=="gaussian"){stop("Use 'bigssa' function for gaussian data.")}
    if(is.na(match(family,c("binomial","Gamma","inverse.gaussian","nb","negbin","poisson")))){stop("Invalid input for 'family'. See documentation.")}
    if(family=="nb"){family="negbin"}
    if(family=="binomial"){
      yrng <- range(yvar)
      if(yrng[1]<0){stop("Response should be between 0 and 1 for family='binomial'.")}
      if(yrng[2]>1){stop("Response should be between 0 and 1 for family='binomial'.")}
      weights <- as.integer(weights)
      if(any(weights<1)){stop("Weights should be total number of trials for family='binomial'.")}
      yvar <- matrix(round(yvar*weights))
    }
    if(family==c("poisson")){
      yvar <- matrix(round(yvar))
      ymin <- min(yvar)
      if(ymin[1]<0){stop("Response must be nonnegative integer for family='poisson'.")}
    }
    if(family==c("negbin")){
      yvar <- matrix(round(yvar))
      ymin <- min(yvar)
      if(ymin[1]<0){stop("Response must be nonnegative integer for family='negbin'.")}
    }
    if(family=="Gamma"){
      ymin <- min(yvar)
      if(ymin[1]<=0){stop("Response must be positive real for family='Gamma'.")}
    }
    if(family=="inverse.gaussian"){
      ymin <- min(yvar)
      if(ymin[1]<=0){stop("Response must be positive real for family='inverse.gaussian'.")}
    }
    
    ### check gcvopts
    if(is.null(gcvopts[1])==FALSE){
      gmatch <- match(c("maxit","gcvtol","alpha","inmaxit","intol","insub"),names(gcvopts))
      if(any(is.na(gmatch))){stop("Invalid input for 'gcvopts'.")}
      gcvopts$maxit <- as.integer(gcvopts$maxit[1])
      gcvopts$gcvtol <- gcvopts$gcvtol[1]
      gcvopts$alpha <- gcvopts$alpha[1]
      gcvopts$inmaxit <- as.integer(gcvopts$inmaxit[1])
      gcvopts$intol <- gcvopts$intol[1]
      gcvopts$insub <- as.integer(gcvopts$insub[1])
      ugcv <- unlist(gcvopts)
      if(any(ugcv<=0)){stop("Invalid input for 'gcvopts'.")}
      if(gcvopts$maxit<1L){stop("Invalid input for gcvopts$maxit")}
      if(gcvopts$inmaxit<1L){stop("Invalid input for gcvopts$inmaxit")}
      if(gcvopts$insub<1L){stop("Invalid input for gcvopts$insub")}
    }
    
    ### check formula
    testform <- strsplit(as.character(formula),":")
    for(k in 1:length(testform)){
      if(length(testform[[k]])>1L){stop("Formula syntax 'x1:x2' is not supported (i.e., main effects must accompany interactions). \n  Use 'x1*x2' to include interaction with main effects.")}
    }
    
    ### get constant (for log-likelihood evaluation)
    slogy <- NULL
    if(family=="binomial"){
      slogy <- sum(log(choose(weights,yvar)))
    } else if(family=="poisson"){
      slogy <- sum(apply(yvar,1,function(x) {ifelse(x>0,sum(log(1:x)),0)}))
    } else if(family=="Gamma"){
      slogy <- sum(log(yvar))
    } else if(family=="inverse.gaussian"){
      slogy <- c(sum(1/yvar),sum(log(yvar^3)))
    } else if(family=="negbin"){
      slogy <- sum(lgamma(yvar+1))
    }
    
    ### check order of things for nonparametric predictors
    if(any(type=="prm")){stop("Parametric predictors are not supported. Use 'bigssp' function.")}
    if(nxvar>1){
      snidx <- match(xnames,names(type))
      if(any(is.na(snidx))){stop('Variable names in "type" must match variable names in "formula".')}
      type <- type[snidx]
      if(!is.na(rparm[1])){
        rnidx <- match(xnames,names(rparm))
        if(any(is.na(rnidx))){stop('Variable names in "rparm" must match variable names in "formula".')}
        rparm <- rparm[rnidx]
      }
      if(!is.null(gammas[1])){
        gnidx <- match(xnames,names(gammas))
        if(any(is.na(gnidx))){stop('Variable names in "gammas" must match variable names in "formula".')}
        gammas <- unlist(gammas[gnidx])
        if(any(gammas<=0)){stop('Input "gammas" must be nonnegative smoothing parameters.')}
      }
    } else{
      type <- type[[1]]
      if(is.null(type)){type <- "cub"}
      rparm <- rparm[[1]]
    }
    
    ### collect predictors
    xvars <- vector("list",nxvar)
    for(jj in 1:nxvar){xvars[[jj]] <- mf[,jj+1L]}
    rm(mf)
    names(xvars) <- xnames
    
    ### check predictor types
    flvls <- xrng <- vector("list",nxvar)
    xdim <- vector("integer",nxvar)
    for(k in 1:nxvar){
      if(any(type[[k]]==c("cub","cub0","per","tps"))){
        xvars[[k]] <- as.matrix(xvars[[k]]+0.0)
        xdim[k] <- ncol(xvars[[k]])
        if(type[[k]]=="tps"){
          if(xdim[k]>3){stop(paste("Too many predictors for",xnames[k]))}
        } else{
          if(xdim[k]>1){stop(paste("Too many predictors for",xnames[k]))}
        }
        if(!is.na(rparm[1])){xrng[[k]] <- apply(xvars[[k]],2,range)}
        flvls[[k]] <- NA
      } else if (type[[k]]=="ord"){
        xvars[[k]] <- factor(xvars[[k]],ordered=TRUE)
        flvls[[k]] <- levels(xvars[[k]])
        xvars[[k]] <- matrix(as.integer(xvars[[k]]))
        xrng[[k]] <- matrix(c(1,length(flvls[[k]])),2,1)
        xdim[k] <- 1L
      } else if (type[[k]]=="nom"){
        xvars[[k]] <- as.factor(xvars[[k]])
        flvls[[k]] <- levels(xvars[[k]])
        xvars[[k]] <- matrix(as.integer(xvars[[k]]))
        xrng[[k]] <- matrix(c(1,length(flvls[[k]])),2,1)
        xdim[k] <- 1L
      } else{
        stop('You must set type to either "cub", "cub0", "nom", "per", or "tps" for each effect in formula.')
      }
    } # end for(k in 1:nxvar)
    
    ### check knots and random seed
    if(is.null(nknots)){
      if(nxvar==1L){
        if(xdim[1]==1L){nknots <- 30} else if(xdim[1]==2L){nknots <- 100} else {nknots <- 200}
      } else if(nxvar==2L) {
        if(any(xdim>1)){nknots <- 200} else {nknots <- 100}
      } else {nknots <- 300} # end if(nxvar==1L)
    }
    if(!is.null(rseed) && length(nknots)==1L){set.seed(rseed)}
    theknots <- vector("list",nxvar)
    if(length(nknots)>1){
      kidx <- as.integer(nknots)
      nknots <- length(kidx)
      if(any(kidx<1L) || any(kidx>nunewr)){stop("Input 'nknots' out of data index range.")}
    } else {kidx <- NA}
    
    ### check for rounding
    if(is.na(rparm[1])){
      for(k in 1:nxvar){
        if(any(type[[k]]==c("cub","cub0","per"))){
          xrng[[k]] <- range(xvars[[k]])
          xvars[[k]] <- (xvars[[k]]-xrng[[k]][1])/(xrng[[k]][2]-xrng[[k]][1])
        }
        if(length(kidx)>1){for(k in 1:nxvar){theknots[[k]]=as.matrix(xvars[[k]][kidx,])}}
      } # end for(k in 1:nxvar)
      xorig <- yorig <- worig <- NA
    } else{
      # check rparms
      rpall <- NULL
      for(k in 1:nxvar){rpall <- c(rpall,rparm[[k]])}
      for(k in 1:length(rpall)){
        rplog <- log(c(rpall[k],rpall[k]/2,rpall[k]/5),base=10)
        rpchk <- rep(FALSE,3)
        for(jj in 1:3){rpchk[jj] <- (rplog[jj]==as.integer(rplog[jj]))}
        if(any(rpchk)==FALSE){stop("Must set input 'rparm' such that rparm=a*(10^-b) with a in {1,2,5} and b>=1 (integer).")}
      }
      # save original variables
      xorig <- xvars
      yorig <- yvar
      worig <- weights
      # get rounded point indices
      gvec <- matrix(1,ndpts,1)
      kconst <- 1
      for(k in 1:nxvar){
        if(type[[k]]=="nom"){
          gvec <- gvec + kconst*(xvars[[k]]-1L)
          kconst <- kconst*xrng[[k]][2]
        } else if(type[[k]]=="ord"){
          gvec <- gvec + kconst*(xvars[[k]]-1L)
          kconst <- kconst*xrng[[k]][2]
        } else if(type[[k]]=="tps"){
          if(length(rparm[[k]])!=xdim[k]){
            rprep <- rep(rparm[[k]][1],xdim[k])
            if(nxvar==1L){rparm <- list(rprep)} else {rparm[[k]] <- rprep}
          }
          rxrng <- xrng[[k]]
          for(j in 1:xdim[k]){
            rxrng[,j] <- round(rxrng[,j]/rparm[[k]][j])*rparm[[k]][j]
            gvec <- gvec + kconst*round((xvars[[k]][,j]-rxrng[1,j])/rparm[[k]][j])
            kconst <- kconst*round(1+(rxrng[2,j]-rxrng[1,j])/rparm[[k]][j])
            xvars[[k]][,j] <- as.matrix(round(xvars[[k]][,j]/rparm[[k]][j]))*rparm[[k]][j]
          }
        } else{
          gidx <- round(((xvars[[k]]-xrng[[k]][1])/(xrng[[k]][2]-xrng[[k]][1]))/rparm[[k]])
          gvec <- gvec + kconst*gidx
          kconst <- kconst*round(1+1/rparm[[k]])
          xvars[[k]] <- as.matrix(gidx*rparm[[k]]);  rm(gidx)
        }
      } # end for(k in 1:nxvar)
      # account for weights
      if(wtchk){
        rmw <- round(range(weights)/.01)*.01
        gvec <- gvec + kconst*round((weights-rmw[1])/.01)
        kconst <- kconst*round(1+(rmw[2]-rmw[1])/.01)
        weights <- round(weights/.01)*.01
      }
      
      # get unique points, frequencies, and sums
      gvec <- as.factor(gvec)
      glindx <- split(cbind(1:ndpts,yvar),gvec)
      if(!is.na(kidx[1])){for(k in 1:nxvar){theknots[[k]] <- as.matrix(xvars[[k]][kidx,])}}
      fs <- matrix(unlist(lapply(glindx,unifqsumg)),ncol=4,byrow=TRUE)
      for(k in 1:nxvar){xvars[[k]] <- as.matrix(xvars[[k]][fs[,1],])}
      yvar <- fs[,2]
      fweights <- fs[,3]
      yty <- fs[,4]
      nunewr <- length(yvar)
      if(wtchk){
        weights <- weights[fs[,1]]
        if(family=="binomial"){yty <- yty/(weights^2)}
      }
    } # end if(is.na(rparm[1]))
    
    ### get knot indices
    if(length(kidx)==1L){
      nknots <- as.integer(min(c(nunewr,nknots)))
      if(nknots<1){stop("Input 'nknots' must be positive integer.")}
      kidx <- sample(nunewr,nknots,prob=(fweights/ndpts))
      for(k in 1:nxvar){theknots[[k]] <- as.matrix(xvars[[k]][kidx,])}
    } 
    
    ### check knots
    if(knotcheck){
      matknots <- NULL
      for(k in 1:nxvar){matknots <- cbind(matknots,theknots[[k]])}
      matknots <- unique(matknots)
      if(nrow(matknots)<nknots){
        if(nxvar>1){
          csdim <- c(0,cumsum(xdim))
          for(k in 1:nxvar){theknots[[k]] <- as.matrix(matknots[,(csdim[k]+1):(csdim[k]+xdim[k])])}
        } else {theknots[[1]] <- as.matrix(matknots)}
        nknots <- nrow(matknots)
      }
    }
    
    ### make marginal reproducing kernel matrices
    rks <- makerkm(xvars,type,theknots,xrng)
    
    ### collect output
    ylist <- list(yvar)
    names(ylist) <- xynames[1]
    ssgmk <- list(et=et,n=c(nunewr,ndpts),xvars=c(xvars,ylist),nxvar=nxvar,
                  yty=yty,fweights=fweights,type=type,xdim=xdim,theknots=theknots,
                  nknots=nknots,lambdas=lambdas,rks=rks[1:4],gcvopts=gcvopts,
                  xorig=xorig,yorig=yorig,se.lp=se.lp,skip.iter=skip.iter,
                  ysm=ysm,rparm=rparm,xrng=xrng,flvls=flvls,tpsinfo=rks$tpsinfo,
                  formula=formula,gammas=gammas,weights=weights,worig=worig,family=family,
                  slogy=slogy,dispersion=dispersion,gcvtype=gcvtype)
    class(ssgmk) <- "makessg"
    return(ssgmk)
    
  }