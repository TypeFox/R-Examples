##
## MCMC sampling for the AR models
## time.data format: col-1: year, col-2: day
##
spAR.Gibbs<-function(formula, data=parent.frame(), time.data, coords, 
           priors=NULL, initials=NULL, nItr, nBurn=0, report=1, 
           tol.dist=2, distance.method="geodetic:km", cov.fnc="exponential",
           scale.transform="NONE", spatial.decay, X.out=TRUE, Y.out=TRUE)
{
    start.time<-proc.time()[3]
  # 
  #
    if(nBurn >= nItr){
         stop(paste("\n Error: iterations < nBurn\n Here, nBurn = ",nBurn," and iterations = ",nItr,"."))
    }
  #  
  #  
    if (missing(formula)) {
         stop("\n Error: formula must be specified \n")
    }
  #
    if (class(formula) != "formula") {
         stop("\n Error: equation must be in formula-class \n ...")
    }
   #
         XY <- Formula.matrix(formula, data)
         Y <- XY[[1]]
         X <- as.matrix(XY[[2]])
         x.names <- XY[[3]]
         Xsp <- XY[[4]]
         x.names.sp <- XY[[5]]
         Xtp <- XY[[6]]
         x.names.tp <- XY[[7]]
         if((!is.null(x.names.sp)) | (!is.null(x.names.tp))){
           stop("\n## \n# Error: spatially and/or temporally varying approach is not available for the AR model\n ##\n")
         }
   #
    if (missing(coords)) {
         stop("\n Error: need to specify the coords \n")
    }
    if ( !is.matrix(coords) ) {
         stop("\n Error: coords must be a (n x 2) matrix of xy-coordinate locations \n")
    }
    if ( (!is.numeric(coords[,1])) | (!is.numeric(coords[,2]))) {
         stop("\n Error: coords columns should be numeric \n")
    }
   #
     method <- distance.method
     spT.check.sites.inside(coords, method, tol=tol.dist)
   #
   if(method=="geodetic:km"){
     coords.D <- as.matrix(spT.geodist(Lon=coords[,1],Lat=coords[,2], KM=TRUE))
   }
   else if(method=="geodetic:mile"){
     coords.D <- as.matrix(spT.geodist(Lon=coords[,1],Lat=coords[,2], KM=FALSE))
   }
   else {
     coords.D <- as.matrix(dist(coords, method, diag = TRUE, upper = TRUE))
   }
   #
   # check time.data
   if(is.null(time.data)){
     #time.data<-c(1,0,length(Y)/length(coords[,1]))
     time.data<-list(1,length(Y)/length(coords[,1]))
   }
   else{
     time.data<-time.data
   }
   #
   #
         n <- length(coords[,1])            # number of sites
         r <- time.data[[1]]                  # number of years
         T <- time.data[[2]]                  # number of days
         # check for T
         if(r > 1){ 
            if(length(T) != r){         
              T<-rep(T,r) 
            }
         }
         #  
         # checking unequal T
         if(length(T) > 1){
           rT <- sum(T)
         }
         else{
           rT <- r*T
         }
         N <- n*rT

    #
    if (N != length(Y)) {
         stop(" Error: Years, Months, and Days are misspecified,\n i.e., total number of observations in the data set should be equal to N\n  : N = n * r * T \n   where, N = total number of observations in the data,\n          n = total number of sites,\n          r = total number of years,\n          T = total number of days.\n# Check the function spT.time.\n#\n")
    }
    # 
         p <- length(x.names)          # number of covariates
    #
         priors<-priors.checking.ar(priors)
    #
         shape_e<-N/2+priors$prior_a
         shape_eta<-N/2+priors$prior_a
         shape_0<-n/2+priors$prior_a
    #

         zm <- matrix(Y,rT,n)
         #cat("\n-- Initial Imputation using Amelia --\n")
         #zm <- imputation.z(zm)
         zm <- apply(zm,1,median,na.rm=TRUE)
         zm <- rep(zm,n)
         zm <- cbind(Y,c(zm))
         zm[is.na(zm[,1]),1] <- zm[is.na(zm[,1]),2]
         zm[is.na(zm[,1]),1] <- median(zm[,2],na.rm=TRUE)
    #
         flag <- matrix(NA,n*rT,2)
         flag[,1] <- c(Y)
         flag[!is.na(flag[,1]),2] <- 0
         flag[is.na(flag[,1]),2] <- 1
    #
    if(cov.fnc=="exponential"){
         cov <- 1
    }
    else if(cov.fnc=="gaussian"){
         cov <- 2
    }
    else if(cov.fnc=="spherical"){
         cov <- 3
    }
    else if(cov.fnc=="matern"){
         cov <- 4
    }
    else{
         stop("\n Error: cov.fnc is not correctly specified \n")
    }
    #
    if(scale.transform=="NONE"){
         zm[,1] <- zm[,1]
    }
    else if(scale.transform=="SQRT"){
         zm[,1] <- sqrt(zm[,1])
    }
    else if(scale.transform=="LOG"){
         zm[,1] <- log(zm[,1])
    }
    else{
         stop("\n Error: scale.transform is not correctly specified \n")
    }
    #
      initials<-initials.checking.ar(initials,zm[,1],X,Xsp,n,r,T,coords.D)
    #
         o <- zm[,1]
    #
    #
    if(spatial.decay$type=="FIXED"){
         spdecay <- 1
		 if(is.null(spatial.decay$value)){
		 spatial.decay$value <- (3/max(c(coords.D))) 
		 }
         init.phi <- spatial.decay$value 
         tuning <- 0; phis<-0; phik<-0; 
		 phi_a <- 0; phi_b <- 0; 
    }
    else if(spatial.decay$type=="DISCRETE"){
         spdecay <- 2
         init.phi <- initials$phi 
         tuning <-0; 
         phis<-spatial.decay$value; 
         phik<-spatial.decay$segments;
		 phi_a<-0; phi_b<-0
    }
    else if(spatial.decay$type=="MH"){
         spdecay <- 3
         init.phi <- initials$phi 
         tuning <- spatial.decay$tuning
         phis<-0; phik<-0; 
		 phi_a<-spatial.decay$val[1]
		 phi_b<-spatial.decay$val[2]
    }
    else{
         stop("\n Error: spatial.decay is not correctly specified \n")
    }
    # 
      if(length(initials$mu_l) != length(initials$sig_l0)){
         stop("Error: check the parameters with year labels")
      }  
      for(i in 1:length(initials$mu_l)){
      if(is.na(initials$mu_l[i])){
         stop("Error: mu_l must be specified correctly")
      }
      if(is.na(initials$sig_l0[i])){
         stop("Error: sig2l must be specified correctly")
      }
      }
    #
      if (length(initials$mu_l) != r){
         stop("Error: need to specify correct number of years for mu.")
      }
      if (length(initials$beta) != p){
         stop("Error: need to specify correct number of parameters for beta.")
      }
    #    
    if(length(x.names.sp)==0){
    # non-spatial beta
      out<-.C('GIBBS_ar',as.double(flag[,2]),as.integer(nItr), 
           as.integer(nBurn), as.integer(n),as.integer(T),as.integer(r),
           as.integer(rT),as.integer(p),as.integer(N),as.integer(report),
           as.integer(cov),as.integer(spdecay),as.double(shape_e),
           as.double(shape_eta),as.double(shape_0),
		   as.double(phi_a),as.double(phi_b),
		   as.double(priors$prior_a),as.double(priors$prior_b),
		   as.double(priors$prior_sig),as.double(init.phi), 
           as.double(tuning),as.double(phis),as.integer(phik),
           as.double(coords.D),as.integer(1),as.double(initials$sig2eps),
           as.double(initials$sig2eta),as.double(initials$sig_l0),as.double(initials$mu_l),
           as.double(initials$rho),as.double(initials$beta),
           as.double(X),as.double(zm[,1]),as.double(o), 
           phip=double(nItr),accept=double(1),nup=double(nItr),sig_ep=double(nItr),sig_etap=double(nItr), 
           rhop=double(nItr),betap=matrix(double(nItr*p),p,nItr),mu_lp=matrix(double(r*nItr),r,nItr),
           sig_l0p=matrix(double(r*nItr),r,nItr),op=matrix(double(nItr*N),N,nItr),wp=matrix(double(nItr*N),N,nItr),
           fit=matrix(double(2*N),N,2),gof=double(1),penalty=double(1))[36:49]
    }
    else{
         stop("\n#\n## Error: \n#")
    }
    #
      accept <- round(out$accept/nItr*100,2)
    #
      output <- NULL
    #
     if(X.out==TRUE){
      output$X <- X
     }
     if(Y.out==TRUE){
      output$Y <- Y
     }
    #
      output$accept <- accept
      output$call<-formula
    #
           output$phip <- as.matrix(out$phip[(nBurn+1):nItr])
           if(cov==4){
           output$nup <- as.matrix(out$nup[(nBurn+1):nItr])
           }
           output$sig2ep <- as.matrix(out$sig_ep[(nBurn+1):nItr])
           output$sig2etap <- as.matrix(out$sig_etap[(nBurn+1):nItr])
           output$sig2lp <- matrix(out$sig_l0p[1:r,(nBurn+1):nItr],r,length((nBurn+1):nItr))
           output$rhop <- as.matrix(out$rhop[(nBurn+1):nItr])
           output$betap <- matrix(out$betap[1:p,(nBurn+1):nItr],p,length((nBurn+1):nItr))
           if(length(x.names.sp) != 0){          
           output$betasp <- matrix(out$betasp[1:(n*q),(nBurn+1):nItr],n*q,length((nBurn+1):nItr))
           }
           output$mu_lp <- matrix(out$mu_lp[1:r,(nBurn+1):nItr],r,length((nBurn+1):nItr))
           output$op <- out$op[1:N,(nBurn+1):nItr]
           output$wp <- out$wp[1:N,(nBurn+1):nItr]
           dimnames(out$fit)[[2]] <- c("Mean","SD")
           output$fitted <- out$fit
           output$tol.dist<-tol.dist
           output$distance.method<-method
           output$cov.fnc<-cov.fnc
           output$scale.transform<-scale.transform
           output$sampling.sp.decay<-spatial.decay
           output$covariate.names<-x.names
           if(length(x.names.sp) != 0){          
           output$sp.covariate.names<-c(x.names.sp)
           }
           output$Distance.matrix <- coords.D
           output$coords <- coords
           output$n <- n
           output$r <- r
           output$T <- T
           output$p <- p
           output$initials <- initials	
           output$priors <- priors	
           output$gof <- round(out$gof,2)
           output$penalty <- round(out$penalty,2)
           tmp <- matrix(c(output$gof,output$penalty,output$gof+output$penalty),1,3)
           dimnames(tmp)[[2]]<-c("Goodness.of.fit","Penalty","PMCC")
           dimnames(tmp)[[1]]<-c("values:")
           output$PMCC <- tmp
           tmp <- NULL
           output$gof <- NULL
           output$penalty <- NULL
    #
          output$iterations <- nItr	
          output$nBurn <- nBurn	
    #
     rm(out)
    #
      cat("##","\n")
      cat("# nBurn = ",nBurn,", Iterations = ",nItr,".", "\n")
      cat("# Overall Acceptance Rate (phi) = ",output$accept,"%", "\n")
      cat("##","\n")
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   output$computation.time<-comp.time
    #
     #class(output) <- "spAR"
     output
    #
    #
}
##
## Prediction of Z_lt for the AR models using MCMC samples
##
 spAR.prediction<-function(nBurn, pred.data, pred.coords, posteriors,
                  tol.dist, Summary, fitted.X=NULL)
{
       start.time<-proc.time()[3] 
       if(is.null(fitted.X)==FALSE){
          if(!is.matrix(fitted.X)) {
           stop("Error: fitted.X must be a MATRIX of order (N x p),\n  p = number of parameters, N = n*r*T.")
          }
        }
      #
        if (missing(posteriors)) {
           stop("Error: need to provide the posterior MCMC samples.")
        }
      #
        if (is.null(posteriors$X)) {
           stop("Error: need to provide the fitted covariate values.")
        }
        if (!is.null(posteriors$X)) {
           fitted.X <- posteriors$X
        }
      #
           #tol.dist<-posteriors$tol.dist
           method<-posteriors$distance.method
           scale.transform<-posteriors$scale.transform
           coords<-posteriors$coords
      #  	
        if (missing(coords)) {
           stop("Error: need to specify the coords.")
        }
        if (!is.matrix(coords)) {
           stop("Error: coords must be a (n x 2) matrix of xy-coordinate locations.")
        }
        if ( (!is.numeric(coords[,1])) | (!is.numeric(coords[,2]))) {
           stop("\n Error: coords columns should be numeric \n")
        }
        if (missing(pred.coords)) {
           stop("Error: need to specify the prediction coords.")
        }
        if (!is.matrix(pred.coords)) {
           stop("Error: prediction coords must be a (n x 2) matrix of xy-coordinate locations.")
        }
        if ( (!is.numeric(pred.coords[,1])) | (!is.numeric(pred.coords[,2]))) {
           stop("\n Error: prediction coords columns should be numeric \n")
        }
      #
           coords.all <- rbind(coords,pred.coords)
           spT.check.locations(coords, pred.coords, method, tol=tol.dist)
           tn.fitsites <- length(coords[, 1])
           nfit.sites <- 1:tn.fitsites
           tn.predsites <- length(coords.all[, 1]) - tn.fitsites
           npred.sites <- (tn.fitsites + 1):(length(coords.all[, 1]))
      #
      if(posteriors$cov.fnc=="exponential"){
         cov <- 1
      }
      else if(posteriors$cov.fnc=="gaussian"){
         cov <- 2
      }
      else if(posteriors$cov.fnc=="spherical"){
         cov <- 3
      }
      else if(posteriors$cov.fnc=="matern"){
         cov <- 4
      }
      else{
         stop("\n Error: cov.fnc is not correctly specified \n")
      }
      #
      #
      if(method=="geodetic:km") {
          coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=TRUE))
      }
      else if(method=="geodetic:mile"){
           coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=FALSE))
      }
      else{
       coords.D <- as.matrix(dist(coords.all, method, diag = T, upper = T))
      }  
      #
           d12 <- coords.D[nfit.sites, npred.sites]
           d <- coords.D[1:tn.fitsites, 1:tn.fitsites]
           dns <- coords.D[(tn.fitsites+1):length(coords.all[,1]), (tn.fitsites+1):length(coords.all[,1])]
      #
           nItr <- posteriors$iterations
           if((nBurn+posteriors$nBurn) >= nItr){
              stop(paste("\n Error: burn-in >= iterations\n Here, nBurn = ",nBurn+posteriors$nBurn," and iterations = ",nItr,"."))
           }
           nItr <-(nItr-posteriors$nBurn) 
           itt <-(nItr-nBurn)
       # 
           n <- tn.fitsites
           r <- posteriors$r
           T <- posteriors$T
           p <- posteriors$p
           if(length(T)>1){
              rT <- sum(T)
           }
           else{
              rT <- r*T
           }
           N <- n*rT
      #
           nsite <- tn.predsites
           predN <- nsite*rT
      #
    #
    # adding the formula in to the prediction dataset
      if(!is.data.frame(pred.data) & !is.null(pred.data)){
        stop("#\n# Error: pred.data should be in data format\n#")
      }
      call.f<-posteriors$call  
      call.f<-as.formula(paste("tmp~",paste(call.f,sep="")[[3]]))
      if(is.data.frame(pred.data)){
      if((nsite*rT)!=dim(pred.data)[[1]]){
        print("#\n # Check the pred.data \n#\n")
      }
      pred.data$tmp<-rep(1,nsite*rT)
      }
      if(is.null(pred.data)){
        pred.data<-data.frame(tmp=rep(1,nsite*rT))
      }
      pred.x<-Formula.matrix(call.f,data=pred.data)[[2]]
    #
    #
           phip<-posteriors$phip[(nBurn+1):nItr,]
           if(cov==4){
           nup<-posteriors$nup[(nBurn+1):nItr,]
           }
           else{
           nup<-0
           } 
           sig_ep<-posteriors$sig2ep[(nBurn+1):nItr,]
           sig_etap<-posteriors$sig2etap[(nBurn+1):nItr,]
           rhop<-posteriors$rhop[(nBurn+1):nItr,]
           betap<-posteriors$betap[,(nBurn+1):nItr]
           mu_lp<-posteriors$mu_lp[,(nBurn+1):nItr]	
           sig_l0p<-posteriors$sig2lp[,(nBurn+1):nItr]	
           op<-posteriors$op[,(nBurn+1):nItr] 	
      #
       out<-matrix(.C('z_pr_its_ar',as.integer(cov),as.integer(itt),as.integer(nsite), 
            as.integer(n),as.integer(r),as.integer(rT),as.integer(T),
            as.integer(p),as.integer(N),as.double(d), 
            as.double(d12),as.double(phip),as.double(nup),as.double(sig_ep),as.double(sig_etap),
            as.double(sig_l0p),as.double(rhop),as.double(betap),
            as.double(mu_lp),as.double(fitted.X),
            as.double(pred.x),as.double(op),as.integer(1),
            out=double(itt*predN))$out,predN,itt)
      #
       output <- NULL
          #
          if(scale.transform=="NONE"){ 
           output$pred.samples <- out          
          }
          if(scale.transform=="SQRT"){ 
           output$pred.samples <- out^2
          }
          if(scale.transform=="LOG"){ 
           output$pred.samples <- exp(out)
          }
          # 
       out<-NULL
       output$pred.coords <- pred.coords
       output$distance.method <- posteriors$distance.method
       output$Distance.matrix.pred <- dns
       output$cov.fnc <- posteriors$cov.fnc
       output$predN <- predN
      #
        if(Summary == TRUE){
         if(itt < 40){
          cat("##", "\n")
          cat("# Summary statistics are not given, because\n# the number of predicted samples (",itt,") are too small.\n# nBurn = ",nBurn+posteriors$nBurn,". Iterations = ",posteriors$iterations,".", "\n")
          cat("##", "\n")
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   output$computation.time<-comp.time
    #
          #class(output) <- "spAR"
          output
         }
         else { 
          cat("##", "\n")
          cat("# Predicted samples and summary statistics are given.\n# nBurn = ",nBurn+posteriors$nBurn,". Iterations = ",posteriors$iterations,".", "\n")
          cat("##", "\n")
          #
          szp<-spT.Summary.Stat(output$pred.samples[,])
          # 
          output$Mean <- matrix(szp$Mean,rT, nsite)
          output$Median <- matrix(szp$Median,rT, nsite)
          output$SD <- matrix(szp$SD,rT, nsite)
          output$Low <- matrix(szp[,4],rT, nsite)
          output$Up <- matrix(szp[,5],rT, nsite)
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   output$computation.time<-comp.time
    #
          #class(output) <- "spAR"
          output
         }
      }
      else {
          cat("##", "\n")
          cat("# Predicted samples are given.\n# nBurn = ",nBurn+posteriors$nBurn,", Iterations = ",posteriors$iterations,".", "\n")	
          cat("##", "\n")
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   output$computation.time<-comp.time
    #
          #class(output) <- "spAR"
          output
      }
}
##
## K-step forecasts for the AR models (1st function)
##
 spAR.forecast <- function(nBurn, K, fore.data, 
        fore.coords, posteriors, pred.samples.ar=NULL,
        tol.dist, Summary=TRUE)
{
      start.time<-proc.time()[3]
      #
      if (missing(posteriors)) {
           stop("\n# Error: need to provide the posterior MCMC samples.")
      }
      #
           nItr <- posteriors$iterations
           if((nBurn+posteriors$nBurn) >= nItr){
              stop(paste("\n Error: burn-in >= iterations\n Here, nBurn = ",nBurn+posteriors$nBurn," and iterations = ",nItr,"."))
           }
           nItr <-(nItr-posteriors$nBurn) 
           itt <-(nItr-nBurn)
     #
     #
       if(is.null(fore.coords)){
              stop("\n# Error: need to provide fore.coords ")
       }
       if(missing(fore.coords)) {
         stop("\n# Error: need to specify the fore.coords")
       }
       if ( !is.matrix(fore.coords) ) {
         stop("\n# Error: fore.coords must be a matrix of xy-coordinate locations")
       }
       if ( (!is.numeric(fore.coords[,1])) | (!is.numeric(fore.coords[,2]))) {
           stop("\n Error: fore.coords columns should be numeric \n")
       }
      #
      #
           n <- posteriors$n
           r <- posteriors$r
           T <- posteriors$T
           if(length(T)>1){
             rT <- sum(T)
           } 
           else{
             rT <- r*T
           }
           N <- n*rT
           p <- posteriors$p
           coords<-posteriors$coords
      #
      if(posteriors$cov.fnc=="exponential"){
         cov <- 1
      }
      else if(posteriors$cov.fnc=="gaussian"){
         cov <- 2
      }
      else if(posteriors$cov.fnc=="spherical"){
         cov <- 3
      }
      else if(posteriors$cov.fnc=="matern"){
         cov <- 4
      }
      else{
         stop("\n# Error: cov.fnc is not correctly specified \n")
      }
      #
      #
      if(is.null(pred.samples.ar)[1]) {
           nsite <- dim(fore.coords)[[1]]
           predN <- nsite*rT
           coords.D <- posteriors$Distance.matrix
           coords.f.D<-as.matrix(spT.geodist(Lon=c(fore.coords[,1],coords[,1]),
                                 Lat=c(fore.coords[,2],coords[,2]),KM=TRUE))
           coords.f.D<-coords.f.D[1:nsite,(nsite+1):(nsite+n)] 
      }
      #
      if(!is.null(pred.samples.ar)[1]){
        if (!is.matrix(pred.samples.ar$pred.samples[,])) {
           stop("\n# Error: MCMC samples of predictions must be a (samples x iterations) matrix.")
        }
        #
        if(itt != dim(pred.samples.ar$pred.samples)[[2]]){
           stop("\n# Error: number of nBurn should be same as the number of nBurn in the predictions.")
        }  	
           predN <- length(pred.samples.ar$pred.samples[,])/itt
           nsite <- predN/(rT)
           if(nsite!=dim(fore.coords)[[1]]){
              stop("#\n# Error: predAR should be equal to NULL for insite temporal forecast. \n  Or correctly define the forecast coords.\n#")
           }
           coords.D <- posteriors$Distance.matrix
           coords.f.D<-as.matrix(spT.geodist(Lon=c(fore.coords[,1],coords[,1]),
                                 Lat=c(fore.coords[,2],coords[,2]),KM=TRUE))
           coords.f.D<-coords.f.D[1:nsite,(nsite+1):(nsite+n)] 
      }
      #
      #
    #
    # adding the formula in to the forecast dataset
      if(!is.data.frame(fore.data) & !is.null(fore.data)){
        stop("#\n# Error: fore.data should be in data format\n#")
      }
      call.f<-posteriors$call  
      call.f<-as.formula(paste("tmp~",paste(call.f,sep="")[[3]]))
      if(is.data.frame(fore.data)){
        if((nsite*r*K)!=dim(fore.data)[[1]]){
          stop("\n# Check the newdata and/or newcoords for forecasts and/or spT.time function#\n")
          #stop("\n# Check the fore.data and/or spT.time function#\n")
        }
        fore.data$tmp<-rep(1,nsite*r*K)
      }
      if(is.null(fore.data)){
        fore.data<-data.frame(tmp=rep(1,nsite*r*K))
      }
      fore.x<-Formula.matrix(call.f,data=fore.data)[[2]]
    #
    #
           phip<-posteriors$phip[(nBurn+1):nItr]
           if(cov==4){
           nup<-posteriors$nup[(nBurn+1):nItr,]
           }
           else{
           nup<-0
           } 
           sig_ep<-posteriors$sig2ep[(nBurn+1):nItr]
           sig_etap<-posteriors$sig2etap[(nBurn+1):nItr]
           rhop<-posteriors$rhop[(nBurn+1):nItr]
           betap<-matrix(posteriors$betap,p,nItr)
           betap<-betap[,(nBurn+1):nItr]
      #
        w<-apply(posteriors$wp,1,median)
        w<-matrix(w,rT,n)
        w<-apply(w,2,median)
        #w<-apply(posteriors$wp[,(nBurn+1):nItr],2,function(x)apply(matrix(x,rT,n),2,median))
      #
      if(is.null(pred.samples.ar)[1]){
      out<-matrix(.C('zlt_fore_ar_its_anysite',as.integer(cov),
           as.integer(itt),as.integer(K),as.integer(nsite),as.integer(n),
           as.integer(r),as.integer(p),as.integer(rT),as.integer(T),
           as.integer(r*K),as.integer(nsite*r*K),as.double(coords.D),as.double(t(coords.f.D)),
           as.double(phip),as.double(nup),as.double(sig_ep),as.double(sig_etap),
           as.double(rhop),as.double(fore.x),as.double(betap),
           as.double(posteriors$op[,]),as.double(w),
           as.integer(1),foreZ=double(nsite*r*K*itt))$foreZ,nsite*r*K,itt)
      }
      #
      if(!is.null(pred.samples.ar)[1]){
      #
          if(posteriors$scale.transform=="NONE"){ 
pred.samples<-pred.samples.ar$pred.samples[,]
          }
          if(posteriors$scale.transform=="SQRT"){ 
pred.samples<-sqrt(pred.samples.ar$pred.samples[,])          
          }
          if(posteriors$scale.transform=="LOG"){ 
pred.samples<-log(pred.samples.ar$pred.samples[,])
          }
      # 
      out<-matrix(.C('zlt_fore_ar_its_anysite',as.integer(cov),
           as.integer(itt),as.integer(K),as.integer(nsite),as.integer(n),
           as.integer(r),as.integer(p),as.integer(rT),as.integer(T),
           as.integer(r*K),as.integer(nsite*r*K),as.double(coords.D),as.double(t(coords.f.D)),
           as.double(phip),as.double(nup),as.double(sig_ep),as.double(sig_etap),
           as.double(rhop),as.double(fore.x),as.double(betap),
           as.double(pred.samples),as.double(w),
           as.integer(1),foreZ=double(nsite*r*K*itt))$foreZ,nsite*r*K,itt)
      }
      #
      output <- NULL
          #
          if(posteriors$scale.transform=="NONE"){ 
      output$fore.samples <- out
          }
          if(posteriors$scale.transform=="SQRT"){ 
      output$fore.samples <- out^2
          }
          if(posteriors$scale.transform=="LOG"){ 
      output$fore.samples <- exp(out)
          }
          # 
      out<-NULL
      output$fore.coords <- fore.coords
      output$distance.method<-posteriors$distance.method  
      output$cov.fnc<-posteriors$cov.fnc  

      output$obsData<-matrix(posteriors$Y,rT,n)  
      output$fittedData<-matrix(posteriors$fitted[,1],rT,n) 
      if(posteriors$scale.transform=="SQRT"){output$fittedData<-output$fittedData^2}
      else if(posteriors$scale.transform=="LOG"){output$fittedData<-exp(output$fittedData)}
      else {output$fittedData<-output$fittedData}
      output$residuals<-matrix(c(output$obsData)-c(output$fittedData),rT,n)

      #
      if(Summary == TRUE){
         if(itt < 40){
          cat("##", "\n")
          cat("# Summary statistics are not given for the prediction sites,\n#   because the number of forecast samples (",itt,") are too small.\n# nBurn = ",nBurn,". Iterations = ",nItr,".", "\n")
          cat("##", "\n")
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   output$computation.time<-comp.time
    #
          #class(output) <- "spAR"
          output
         }
         else { 
          cat("##", "\n")
          cat("# Forecast samples and summary statistics are given.\n# nBurn = ",nBurn,". Iterations = ",nItr,".", "\n")
          cat("##", "\n")
          #
          szp<-spT.Summary.Stat(output$fore.samples[,])
          # 
          output$Mean <- matrix(szp$Mean,r*K, nsite)
          output$Median <- matrix(szp$Median,r*K, nsite)
          output$SD <- matrix(szp$SD,r*K, nsite)
          output$Low <- matrix(szp[,4],r*K, nsite)
          output$Up <- matrix(szp[,5],r*K, nsite)
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   output$computation.time<-comp.time
    #
          #class(output) <- "spAR"
          output
         }
      }
      else{
          cat("##", "\n")
          cat("# Forecast samples are given.\n# nBurn = ",nBurn,", Iterations = ",itt,".", "\n")
          cat("##", "\n")
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   output$computation.time<-comp.time
    #
          #class(output) <- "spAR"
          output
      }
}
##
## MCMC fit and predictions for the AR models
##
spAR.MCMC.Pred<-function(formula, data=parent.frame(), time.data, 
            coords, pred.coords, priors, initials,
            pred.data, nItr, nBurn=0, report=1, tol.dist=2,
            distance.method="geodetic:km", cov.fnc="exponential",
            scale.transform="NONE", spatial.decay, annual.aggregation="NONE")
{
    start.time<-proc.time()[3]
  #
    if(nBurn >= nItr){
         stop(paste("\n Error: iterations < nBurn\n Here, nBurn = ",nBurn," and iterations = ",nItr,"."))
    }
  #  
    if (missing(formula)) {
         stop("\n Error: formula must be specified \n")
    }
  #
    if (class(formula) != "formula") {
         stop("\n Error: equation must be in formula-class \n ...")
    }
  #
         XY <- Formula.matrix(formula, data)
         Y <- XY[[1]]
         X <- as.matrix(XY[[2]])
         x.names <- XY[[3]]
         Xsp <- XY[[4]]
         x.names.sp <- XY[[5]]

  #
    if (missing(coords)) {
         stop("Error: need to specify the coords")
    }
    if ( !is.matrix(coords) ) {
         stop("\n Error: coords must be a (n x 2) matrix of xy-coordinate locations \n")
    }
    if ( (!is.numeric(coords[,1])) | (!is.numeric(coords[,2]))) {
         stop("\n Error: coords columns should be numeric \n")
    }
    #
      method <- distance.method
      spT.check.sites.inside(coords, method, tol=tol.dist)
   #
    if(method=="geodetic:km"){
      coords.D <- as.matrix(spT.geodist(Lon=coords[,1],Lat=coords[,2], KM=TRUE))
    }
    else if(method=="geodetic:mile"){
      coords.D <- as.matrix(spT.geodist(Lon=coords[,1],Lat=coords[,2], KM=FALSE))
    }
    else {
       coords.D <- as.matrix(dist(coords, method, diag = T, upper = T))
    }  
   #
   # check time.data
   if(is.null(time.data)){
     #time.data<-c(1,0,length(Y)/length(coords[,1]))
     time.data<-list(1,length(Y)/length(coords[,1]))
   }
   else{
     time.data<-time.data
   }
   #
         n <- length(coords[,1])            # number of sites
         r <- time.data[[1]]                  # number of years
         T <- time.data[[2]]                  # number of days
         # check for T
         if(r > 1){ 
            if(length(T) != r){         
              T<-rep(T,r) 
            }
         }
         #  
         # checking unequal T
         if(length(T) > 1){
           rT <- sum(T)
         }
         else{
           rT <- r*T
         }
         N <- n*rT
   #
    if (N != length(Y)) {
         stop(" Error: Years, Months, and Days are misspecified,\n i.e., total number of observations in the data set should be equal to N\n  : N = n * r * T \n   where, N = total number of observations in the data,\n          n = total number of sites,\n          r = total number of years,\n          T = total number of days.\n# Check the function spT.time.\n#\n")
    }
    #
         p <- length(x.names)          # number of covariates
    #
         priors<-priors.checking.ar(priors)
    #
         shape_e<-N/2+priors$prior_a
         shape_eta<-N/2+priors$prior_a
         shape_0<-n/2+priors$prior_a
    #
         zm <- matrix(Y,rT,n)
         zm <- apply(zm,1,median,na.rm=TRUE)
         zm <- rep(zm,n)
         zm <- cbind(Y,zm)
         zm[is.na(zm[,1]),1]=zm[is.na(zm[,1]),2]
         zm <- zm[,1]
    #
         flag <- matrix(NA,n*rT,2)
         flag[,1] <- c(Y)
         flag[!is.na(flag[,1]),2] <- 0
         flag[is.na(flag[,1]),2] <- 1
         flag <- flag[,2]
    #
    #
    if(cov.fnc=="exponential"){
         cov <- 1
    }
    else if(cov.fnc=="gaussian"){
         cov <- 2
    }
    else if(cov.fnc=="spherical"){
         cov <- 3
    }
    else if(cov.fnc=="matern"){
         cov <- 4
    }
    else{
         stop("\n Error: cov.fnc is not correctly specified \n")
    }
    #
    if (scale.transform == "NONE") {
         zm <- zm
         trans <- 0
         scale.transform = "NONE"
    }
    else if (scale.transform == "SQRT") {
         zm <- sqrt(zm)
         trans <- 1
         scale.transform = "SQRT"
    }
    else if (scale.transform == "LOG") {
         zm <- log(zm)
         trans <- 2
         scale.transform = "LOG"
    }
    else{
         stop("\n Error: scale.transform is not correctly specified \n")
    }
    #
      initials<-initials.checking.ar(initials,zm,X,Xsp,n,r,T,coords.D)
    #
         o <- zm
    #
    #
    if(spatial.decay$type=="FIXED"){
         spdecay <- 1
		 if(is.null(spatial.decay$value)){
		 spatial.decay$value <- (3/max(c(coords.D))) 
		 }
         init.phi <- spatial.decay$value 
         tuning <- 0; phis<-0; phik<-0; 
		 phi_a <- 0; phi_b <- 0; 
    }
    else if(spatial.decay$type=="DISCRETE"){
         spdecay <- 2
         init.phi <- initials$phi 
         tuning <-0; 
         phis<-spatial.decay$value; 
         phik<-spatial.decay$segments;
		 phi_a<-0; phi_b<-0
    }
    else if(spatial.decay$type=="MH"){
         spdecay <- 3
         init.phi <- initials$phi 
         tuning <- spatial.decay$tuning
         phis<-0; phik<-0; 
		 phi_a<-spatial.decay$val[1]
		 phi_b<-spatial.decay$val[2]
    }
    else{
         stop("\n Error: spatial.decay is not correctly specified \n")
    }
    # 
    #
      if(length(initials$mu_l) != length(initials$sig_l0)){
         stop("Error: check the parameters with year labels")
      }  
      for(i in 1:length(initials$mu_l)){
      if(is.na(initials$mu_l[i])){
         stop("Error: mu_l must be specified correctly")
      }
      if(is.na(initials$sig_l0[i])){
         stop("Error: sig2l must be specified correctly")
      }
      }
    #
      if (length(initials$mu_l) != r){
         stop("Error: need to specify correct number of years for mu.")
      }
      if (length(initials$beta) != p){
         stop("Error: need to specify correct number of parameters for beta.")
      }
    #
      #if (missing(nBurn)) {
      #   stop("Error: nBurn must be specified")
      #}
    #
    # prediction part
    #
        if (missing(pred.coords)) {
           stop("Error: need to specify the prediction coords.")
        }
        if (!is.matrix(pred.coords)) {
           stop("Error: prediction coords must be a (n x 2) matrix of xy-coordinate locations.")
        }
        if ( (!is.numeric(pred.coords[,1])) | (!is.numeric(pred.coords[,2]))) {
           stop("\n Error: prediction coords columns should be numeric \n")
        }
      #
           coords.all <- rbind(coords,pred.coords)
           spT.check.locations(coords, pred.coords, method, tol=tol.dist)
           tn.fitsites <- length(coords[, 1])
           nfit.sites <- 1:tn.fitsites
           tn.predsites <- length(coords.all[, 1]) - tn.fitsites
           npred.sites <- (tn.fitsites + 1):(length(coords.all[, 1]))
     #                                                                         
      if(method=="geodetic:km"){
           coords.D.all <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=TRUE))
      }
      else if(method=="geodetic:mile"){
           coords.D.all <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=FALSE))
      }
      else {
           coords.D.all<- as.matrix(dist(coords.all, method, diag = T, upper = T))
      }
      #
           d12 <- coords.D.all[nfit.sites, npred.sites]
           dns <- coords.D.all[(tn.fitsites+1):length(coords.all[,1]), (tn.fitsites+1):length(coords.all[,1])]
      #
           nsite <- tn.predsites
           predN <- nsite*rT
      #
    #
    # adding the formula in to the prediction dataset
      if(!is.data.frame(pred.data) & !is.null(pred.data)){
        stop("#\n# Error: pred.data should be in data format\n#")
      }
      call.f<-formula  
      call.f<-as.formula(paste("tmp~",paste(call.f,sep="")[[3]]))
      if(is.data.frame(pred.data)){
      if((nsite*rT)!=dim(pred.data)[[1]]){
        print("#\n # Check the pred.data \n#\n")
      }
      pred.data$tmp<-rep(1,nsite*rT)
      }
      if(is.null(pred.data)){
        pred.data<-data.frame(tmp=rep(1,nsite*rT))
      }
      pred.x<-Formula.matrix(call.f,data=pred.data)[[2]]
    #
    #
       if(annual.aggregation=="NONE"){
         aggtype<-0
       }
       else if(annual.aggregation=="ave"){
         aggtype<-1
       }
       else if(annual.aggregation=="an4th"){
         aggtype<-2
       }
       else if(annual.aggregation=="w126"){
         aggtype<-3
         if(T != 214){
          stop("\n# Error: day/time should have 214 values for annual.aggregation = w126.\n")        
         }
       }
       else{
        stop("Error: correctly define annual.aggregation.")
       }         
    #
      #
      ########
      #
        out <- NULL
        out <- .C("GIBBS_sumpred_txt_ar", as.integer(aggtype),as.double(flag), as.integer(nItr),
            as.integer(nBurn), as.integer(n), as.integer(T), as.integer(r), 
            as.integer(rT), as.integer(p), as.integer(N), as.integer(report),
            as.integer(cov), as.integer(spdecay), as.double(shape_e), as.double(shape_eta), as.double(shape_0), 
			as.double(phi_a), as.double(phi_b),
			as.double(priors$prior_a), as.double(priors$prior_b), 
            as.double(priors$prior_sig), as.double(init.phi), as.double(tuning),
            as.double(phis), as.integer(phik), as.double(coords.D), 
            as.integer(1), as.double(initials$sig2eps), as.double(initials$sig2eta),
            as.double(initials$sig_l0), as.double(initials$mu_l), as.double(initials$rho), 
            as.double(initials$beta), as.double(X), as.double(zm), as.double(o),
            as.integer(nsite), as.integer(predN), as.double(d12), as.double(pred.x), 
            as.integer(trans), accept=double(1), gof=double(1),penalty=double(1))[42:44] 
     #
      out$accept <- round(out$accept/nItr*100,2)
      out$call<-formula
     #
      cat("##","\n")
      cat("# MCMC output and predictions are given in text format ... \n")
      cat("# nBurn = ",nBurn,". Iterations = ",nItr,".", "\n")
      cat("# Acceptance rate: (phi) = ",out$accept,"%", "\n")
      cat("##","\n")
     #
     #
           tmp<-read.table('OutAR_Values_Parameter.txt',sep='',header=FALSE)
           tmp<-tmp[(nBurn+1):nItr,]
           tmp<-spT.Summary.Stat(t(tmp))
           if(cov==4){
           tmp<-tmp[1:length(c(x.names,"rho","sig2eps","sig2eta","phi","nu")),]
           row.names(tmp)<- c(x.names,"rho","sig2eps","sig2eta","phi","nu")
           }
           else{
           tmp<-tmp[1:length(c(x.names,"rho","sig2eps","sig2eta","phi")),]
           row.names(tmp)<- c(x.names,"rho","sig2eps","sig2eta","phi")
           } 
           #row.names(tmp)<- c(x.names,"rho","sig2eps","sig2eta", paste("sig2",1:r),paste("mu",1:r),"phi") # when to get summary of sigl and mul
           out$parameters<-round(tmp, 4)
           tmp<-NULL
#
           tmp<-read.table('OutAR_Stats_FittedValue.txt', sep=',',header=FALSE)
           names(tmp)<- c("Mean","SD")
           out$fitted<-round(tmp, 4)
           tmp<-NULL
#
           tmp<-read.table('OutAR_Stats_PredValue.txt', sep=',',header=FALSE)
           names(tmp)<- c("Mean","SD")
           out$prediction<-round(tmp, 4)
           tmp<-NULL
#
       if(annual.aggregation=="ave"){
           tmp<-read.table('OutAR_Annual_Average_Prediction.txt', sep='',header=FALSE)
           tmp<-spT.Summary.Stat(t(tmp))
           out$an.agr.pred.average<-round(tmp, 4)
           tmp<-NULL
       }   
       if(annual.aggregation=="an4th"){
           tmp<-read.table('OutAR_Annual_4th_Highest_Prediction.txt', sep='',header=FALSE)
           tmp<-spT.Summary.Stat(t(tmp))
           out$an.agr.pred.an4th<-round(tmp, 4)
           tmp<-NULL
       }   
       if(annual.aggregation=="w126"){
           tmp<-read.table('OutAR_Annual_w126_Prediction.txt', sep='',header=FALSE)
           tmp<-spT.Summary.Stat(t(tmp))
           out$an.agr.pred.w126<-round(tmp, 4)
           tmp<-NULL
       }   
#
           out$tol.dist <- tol.dist	
           out$distance.method <- distance.method
           out$cov.fnc <- cov.fnc
           out$scale.transform <- scale.transform
           out$sampling.sp.decay<-spatial.decay
           out$covariate.names<-x.names
           out$gof<-round(out$gof,2)
           out$penalty<-round(out$penalty,2) 
           tmp <- matrix(c(out$gof,out$penalty,out$gof+out$penalty),1,3)
           dimnames(tmp)[[2]] <- c("Goodness.of.fit","Penalty","PMCC")
           dimnames(tmp)[[1]] <- c("values:")
           out$PMCC <- tmp
           tmp <- NULL
           out$gof<-NULL
           out$penalty<-NULL 
           out$n <- n
           out$pred.n<-nsite
           out$r <- r
           out$T <- T
           out$Y <- Y
           out$initials <- initials	
           out$priors <- priors 
           out$iterations <- nItr 
           out$nBurn <- nBurn 
      #
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   out$computation.time<-comp.time
    #
       #class(out) <- "spAR"
      #
       out
      #
      #
}
##
## initials checking for AR
##
 initials.checking.ar<-function(x,o,X,Xsp,n,r,T,d){
   #
   # x = initial values
   #
   if(is.null(Xsp)){
   #
   ## for non-spatial beta
     if(is.null(x)){
       x$phi=NULL; x$sig2eps=NULL; x$sig2eta=NULL; x$sig_l0=NULL;
       x$rho=NULL; x$beta=NULL; x$mu_l=NULL;
       x$phi<-(3/max(c(d)))#(-log(0.05)/max(c(d)))
       x$sig2eps <- 0.01; x$sig2eta <- 0.1
       x$mu_l<-rep(NA,r); x$sig_l0<-rep(0.1,r)
       if(length(T)>1){ 
         rT <- sum(T)
       }
       else{ 
         rT <- r*T
       }
       om <- matrix(o,rT,n)
       for(i in 1:r){
         x$mu_l[i] <- mean(om[1+c(0,cumsum(T))[[i]],],na.rm=TRUE)
       }
       o1 <- matrix(NA,rT,n)
       for(i in 1:n){
        for(j in 1:r){
         o1[c((1+c(0,cumsum(T))[[j]]):(c(0,cumsum(T))[[j+1]])),i] <- om[c(0,cumsum(T))[[j]]+c(2:T[[j]],NA),i]
        }
       }
       lm.coef<-lm(c(o) ~ c(o1) + X-1)$coef
       lm.coef[is.na(lm.coef)]<-0
       x$rho<-lm.coef[[1]]
       x$beta<-lm.coef[2:(dim(X)[[2]]+1)]
       x
     }
     else{
       if(is.null(x$phi)){ 
        x$phi<-(-log(0.05)/max(c(d)))      
       }
       if(is.null(x$sig2eps)){
        x$sig2eps<-0.01
       }
       if(is.null(x$sig2eta)){
        x$sig2eta<-0.1
       }
       if(is.null(x$mu_l)){
        x$mu_l<-rep(NA,r); 
        if(length(T)>1){ 
          rT <- sum(T)
        }
        else{ 
          rT <- r*T
        }
        om <- matrix(o,rT,n)
        for(i in 1:r){
         x$mu_l[i] <- mean(om[1+c(0,cumsum(T))[[i]],],na.rm=TRUE)#mean(om[1+T*(i-1),],na.rm=T)
        }
       }
       if(is.null(x$sig_l0)){
        x$sig_l0<-rep(0.1,r)
       }
       if(is.null(x$rho)|is.null(x$beta)){
        if(length(T)>1){ 
          rT <- sum(T)
        }
        else{ 
          rT <- r*T
        }
        om <- matrix(o,rT,n)
        o1 <- matrix(NA,rT,n)
        for(i in 1:n){
         for(j in 1:r){
          o1[c((1+c(0,cumsum(T))[[j]]):(c(0,cumsum(T))[[j+1]])),i] <- om[c(0,cumsum(T))[[j]]+c(2:T[[j]],NA),i]#om[(j-1)*T+c(2:T,NA),i]
         }
        }
        lm.coef<-lm(c(o) ~ c(o1) + X -1)$coef
        lm.coef[is.na(lm.coef)]<-0
        x$rho<-lm.coef[[1]]
        x$beta<-lm.coef[2:(dim(X)[[2]]+1)]
       }
       x
     }
   }
   else{
   #
   ## for spatial beta
     if(is.null(x)){
       x$phi=NULL; x$sig2eps=NULL; x$sig2eta=NULL; x$sig_l0=NULL;
       x$rho=NULL; x$beta=NULL; x$mu_l=NULL;
       x$phi<-(3/max(c(d)))#(-log(0.05)/max(c(d)))
       x$sig2eps <- 0.01; x$sig2eta <- 0.1
       x$mu_l<-rep(NA,r); x$sig_l0<-rep(0.1,r)
       om <- matrix(o,r*T,n)
       for(i in 1:r){
         x$mu_l[i] <- mean(om[1+T*(i-1),],na.rm=TRUE)
       }
       o1 <- matrix(NA,r*T,n)
       for(i in 1:n){
        for(j in 1:r){
         o1[,i] <- om[(j-1)*T+c(2:T,NA),i]
        }
       }
       q<-length(c(Xsp))/(n*r*T)
       dump<-sort(rep(1:n,r*T)) 
       dump<-model.matrix(~factor(dump)-1)
       tdump<-NULL; for(i in 1:q){ tdump <- cbind(tdump,dump) }
       dump<-array(c(Xsp),dim=c(r*T,n,q))
       for(j in 1:q){for(i in 1:n){ tdump[,i+(j-1)*n]<-tdump[,i+(j-1)*n]*dump[,i,j] }}
       lm.coef<-lm(c(o) ~ c(o1) + X + tdump - 1)$coef
       lm.coef[is.na(lm.coef)]<-0
       x$rho<-lm.coef[[1]]
       x$beta<-lm.coef[2:(dim(X)[[2]]+1)]
       x$betasp<-lm.coef[(dim(X)[[2]]+2):(dim(X)[[2]]+1+(n*q))]
       x
     }
     else{
       if(is.null(x$phi)){ 
        x$phi<-(-log(0.05)/max(c(d)))      
       }
       if(is.null(x$sig2eps)){
        x$sig2eps<-0.01
       }
       if(is.null(x$sig2eta)){
        x$sig2eta<-0.1
       }
       if(is.null(x$mu_l)){
        x$mu_l<-rep(NA,r); om <- matrix(o,r*T,n)
        for(i in 1:r){
         x$mu_l[i] <- mean(om[1+T*(i-1),],na.rm=TRUE)
        }
       }
       if(is.null(x$sig_l0)){
        x$sig_l0<-rep(0.1,r)
       }
       if(is.null(x$rho)|is.null(x$beta)){
        o1 <- matrix(NA,r*T,n); om <- matrix(o,r*T,n)
        for(i in 1:n){
         for(j in 1:r){
          o1[,i] <- om[(j-1)*T+c(2:T,NA),i]
         }
        }
        q<-length(c(Xsp))/(n*r*T)
        dump<-sort(rep(1:n,r*T)) 
        dump<-model.matrix(~factor(dump)-1)
        tdump<-NULL; for(i in 1:q){ tdump <- cbind(tdump,dump) }
        dump<-array(c(Xsp),dim=c(r*T,n,q))
        for(j in 1:q){for(i in 1:n){ tdump[,i+(j-1)*n]<-tdump[,i+(j-1)*n]*dump[,i,j] }}
        lm.coef<-lm(c(o) ~ c(o1) + X + tdump - 1)$coef
        lm.coef[is.na(lm.coef)]<-0
        x$rho<-lm.coef[[1]]
        x$beta<-lm.coef[2:(dim(X)[[2]]+1)]
        x$betasp<-lm.coef[(dim(X)[[2]]+2):(dim(X)[[2]]+1+(n*q))]
       }
       x
     }
   }
}
##
## priors checking for AR
##
 priors.checking.ar<-function(x){
     #
     if(is.null(x)){
       x$prior_a<-NULL; x$prior_b<-NULL;
       x$prior_sig<-NULL 
       x$prior_a<-2; x$prior_b<-1;
       x$prior_sig<-10^10
       x
     }
     else{
     if(is.null(x$prior_a)){
       x$prior_a<-2
     } 
     if(is.null(x$prior_b)){
       x$prior_b<-1
     } 
     if(is.null(x$prior_sig)){
       x$prior_sig<-10^10
     }
     x 
     }
}
##
## 
##