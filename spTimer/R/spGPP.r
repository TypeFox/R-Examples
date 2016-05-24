##
## priors checking gpp
##
 priors.checking.gpp<-function(x=NULL,r,p){
     #
     if(is.null(x)){
       x$prior_a=NULL; x$prior_b=NULL;
       x$mu_beta=NULL; x$delta2_beta=NULL;
       x$mu_rho=NULL; x$delta2_rho=NULL;
       x$alpha_l=NULL; x$delta2_l=NULL;
       x$prior_a<-2; x$prior_b<-1;
       x$mu_beta<-rep(0,p); x$delta2_beta<-matrix(0,p,p);
       diag(x$delta2_beta)<-10^10;
       x$mu_rho<-0; x$delta2_rho<-10^10;
       x$alpha_l<-rep(0,r); x$delta2_l<-rep(10^10,r)
       x
     }
     else{
     if(is.null(x$prior_a)){
       x$prior_a<-2
     } 
     if(is.null(x$prior_b)){
       x$prior_b<-1
     } 
     if(is.null(x$mu_beta)){
       x$mu_beta<-rep(0,p)
     }
     if(length(c(x$mu_beta)) != p){
      x$mu_beta<-rep(c(x$mu_beta),p)
     }
     if(is.null(x$delta2_beta)){
       x$delta2_beta<-matrix(0,p,p);
       diag(x$delta2_beta)<-10^10;
     }
     if(length(c(x$delta2_beta)) != (p*p)){
       delta2_beta<-matrix(0,p,p);
       diag(delta2_beta)<-x$delta2_beta;
       x$delta2_beta<-delta2_beta
     }
     if(is.null(x$mu_rho)){
       x$mu_rho<-0; 
     }
     if(is.null(x$delta2_rho)){
       x$delta2_rho<-10^10;
     }
     if(is.null(x$alpha_l)){
       x$alpha_l<-rep(0,r); 
     }
     if(length(c(x$alpha_l)) != r){
       x$alpha_l<-rep(c(x$alpha_l),r)
     }
     if(is.null(x$delta2_l)){
       x$delta2_l<-rep(10^10,r);
     }
     if(length(c(x$delta2_l)) != r){
       x$delta2_l<-rep(c(x$delta2_l),r)
     }
     x 
     }
}
##
## initials checking for the GPP
##
 initials.checking.gpp<-function(x,z,X,n,r,T,d){
     # checking unequal T
     if(length(T)>1){
       rT <- sum(T)
     }
     else{
       rT <- r*T
     } 
     #
     # x = initial values
     #
     if(is.null(x)){
      x$phi=NULL; x$sig2eps=NULL; x$sig2eta=NULL; 
      x$sig2l=NULL; x$rho=NULL; x$beta=NULL; x$mu_l=NULL; 
      x$phi<-(-log(0.05)/max(c(d)))#(3/max(c(d)))
      x$sig2eps <- 0.01; x$sig2eta <- 0.1; x$rho<-0.5
      x$sig2l <- rep(0.1, r);
      x$mu_l<-rep(NA,r); 
      zm <- matrix(z,rT,n)
      for(i in 1:r){
        x$mu_l[i] <- 0
      }
      lm.coef<-lm(c(z) ~ X-1)$coef
      lm.coef[is.na(lm.coef)]<-0
      x$beta<-lm.coef[1:(dim(X)[[2]])]
      x
     }
     else{
      if(is.null(x$phi)){ 
      x$phi<-(3/max(c(d)))#(-log(0.05)/max(c(d)))
      }
      if(is.null(x$sig2eps)){
      x$sig2eps<-0.01
      }
      if(is.null(x$sig2eta)){
      x$sig2eta<-0.1
      }
      if(is.null(x$rho)){
      x$rho<-0.5
      }
      if(is.null(x$sig2l)){
      x$sig2l <- rep(0.1, r);
      }
      if(is.null(x$mu_l)){
      x$mu_l<-rep(NA,r); zm <- matrix(z,rT,n)
      for(i in 1:r){
       x$mu_l[i] <- 0
      }
      }
      if(is.null(x$beta)){
       lm.coef<-lm(c(z) ~ X-1)$coef
       lm.coef[is.na(lm.coef)]<-0
       x$beta<-lm.coef[1:(dim(X)[[2]])]
      }
      x
     }
}
##
## MCMC sampling with only one phi parameter
##
spGPP.Gibbs<-function(formula, data=parent.frame(), time.data,
         knots.coords, coords, priors=NULL, initials=NULL, 
         nItr, nBurn=0, report=1, tol.dist=2, distance.method="geodetic:km", 
         cov.fnc="exponential", scale.transform="NONE", 
         spatial.decay, X.out=TRUE, Y.out=FALSE)
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
    if (class(formula) == "formula") {
         XY <- Formula.matrix(formula, data)
         Y <- XY[[1]]
         X <- as.matrix(XY[[2]])
         x.names <- XY[[3]]
         Xsp <- XY[[4]]
         x.names.sp <- XY[[5]]
         Xtp <- XY[[6]]
         x.names.tp <- XY[[7]]
         #if((!is.null(x.names.sp)) | (!is.null(x.names.tp))){
         #  stop("\n## \n# Error: spatially and/or temporally varying approach is not available for the GPP model\n ##\n")
         #}
         if((!is.null(x.names.tp))){
           stop("\n## \n# Error: dynamic temporally varying approach is not available for the GPP model\n ##\n")
         }
    }
   #
   #
    if (missing(coords)) {
         stop("Error: need to specify the coords")
    }
    if ( !is.matrix(coords) ) {
         stop("Error: coords must be a (n x 2) matrix of xy-coordinate locations")
    }
    if ( (!is.numeric(coords[,1])) | (!is.numeric(coords[,2]))) {
         stop("\n Error: coords columns should be numeric \n")
    }
   #
    if (missing(knots.coords)) {
         stop("Error: need to specify the knots.coords")
    }
    if ( !is.matrix(knots.coords) ) {
         stop("Error: knots.coords must be a (n x 2) matrix of xy-coordinate locations")
    }
    if ( (!is.numeric(knots.coords[,1])) | (!is.numeric(knots.coords[,2]))) {
         stop("\n Error: knots.coords columns should be numeric \n")
    }
   #
   # check time.data
   if(is.null(time.data)){
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
         #
         # checking unequal T
         if(length(T) > 1){
           rT <- sum(T)
         }
         else{
           rT <- r*T
         }
         N <- n*rT
         p <- length(x.names)
         knots <- length(knots.coords[,1]) 
   #
    if(n <= knots){
         stop("Error: n must not be smaller than or equal to knots")
    }
   #
      method <- distance.method
      spT.check.sites.inside(knots.coords, method, tol=tol.dist)
   #
    if(method=="geodetic:km"){
      coords.D.knots <- as.matrix(spT.geodist(Lon=knots.coords[,1],Lat=knots.coords[,2], KM=TRUE))
    }
    else if(method=="geodetic:mile"){
      coords.D.knots <- as.matrix(spT.geodist(Lon=knots.coords[,1],Lat=knots.coords[,2], KM=FALSE))
    }
    else {
      coords.D.knots <- as.matrix(dist(knots.coords, method, diag=TRUE, upper=TRUE))
    }
   #
    if(knots != length(coords.D.knots[,1])){
         stop("Error: knots must be equal to knots.coords")
    }
    #
      all.coords <- rbind(coords,knots.coords)
      spT.check.sites.inside(all.coords, method, tol=tol.dist)
    #
    if(method=="geodetic:km"){
       coords.D.all <- as.matrix(spT.geodist(Lon=all.coords[,1],Lat=all.coords[,2], KM=TRUE))
       coords.D.obs.knots <- coords.D.all[1:length(coords[,1]),(length(coords[,1])+1):length(all.coords[,1])]
    }
    else if(method=="geodetic:mile"){
       coords.D.all <- as.matrix(spT.geodist(Lon=all.coords[,1],Lat=all.coords[,2], KM=FALSE))
       coords.D.obs.knots <- coords.D.all[1:length(coords[,1]),(length(coords[,1])+1):length(all.coords[,1])]
    }
    else {
       coords.D.all <- as.matrix(dist(all.coords, method, diag=TRUE, upper=TRUE))
       coords.D.obs.knots <- coords.D.all[1:length(coords[,1]),(length(coords[,1])+1):length(all.coords[,1])]
    }
    #
    #
    if (N != length(Y)) {
         stop("#\n# Error: Years, Months, and Days are misspecified,\n     i.e., total number of observations in the data set \n     should be equal to N\n  : N = n * r * T \n   where, N = total number of observations in the data,\n          n = total number of sites,\n          r = total number of years,\n          T = total number of days. \n# Check the function spT.time.\n#\n")
    }
    # 
          priors<-priors.checking.gpp(priors,r,p)
    #
         shape_e <- N/2+priors$prior_a
         shape_eta <- (knots*rT)/2+priors$prior_a
         shape_l <- (knots*r)/2+priors$prior_a
    #
         zm <- matrix(Y,rT,n)
         zm <- apply(zm,1,median,na.rm=TRUE)
         zm <- rep(zm,n)
         zm <- cbind(Y,zm)
         zm[is.na(zm[,1]),1] <- zm[is.na(zm[,1]),2]
         zm[is.na(zm[,1]),1] <- median(zm[,2],na.rm=TRUE)
    #
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
    if(scale.transform=="NONE"){
         zm <- zm
         trans <- 0
    }
    else if(scale.transform=="SQRT"){
         zm <- sqrt(zm)
         trans <- 1
    }
    else if(scale.transform=="LOG"){
         zm <- log(zm)
         trans <- 2
    }
    else{
         stop("\n Error: scale.transform is not correctly specified \n")
    }
    #
    #
      initials<-initials.checking.gpp(initials,zm,X,n,r,T,coords.D.all)
    #
    #
    if(spatial.decay$type=="FIXED"){
         spdecay <- 1
		 if(is.null(spatial.decay$value)){
		 spatial.decay$value <- (3/max(c(coords.D.all))) 
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
    if (length(initials$mu_l) != r){
         stop("Error: need to specify correct number of years (r) for mu_l initials.")
    }
    if (length(initials$sig2l) != r){
         stop("Error: need to specify correct number of years (r) for sig2l initials.")
    }
    if (length(initials$beta) != p){
         stop("Error: need to specify correct number (p) of initial parameters for beta.")
    }
    #
    if (length(priors$mu_beta) != p){
         stop("Error: need to specify correct number (p) of mu_beta priors.")
    }
    if (length(priors$delta2_beta) != p*p){
         stop("Error: need to specify correct number (p) of delta2_beta priors.")
    }
    if (length(priors$alpha_l) != r){
         stop("Error: need to specify correct number (r) of alpha_l priors.")
    }
    if (length(priors$delta2_l) != r){
         stop("Error: need to specify correct number (r) of delta2_l priors.")
    }
    #
     tmp <- zm-median(zm)
     tmp<-matrix(tmp,rT,n)
     tmp<-apply(tmp,1,mean)
     w<-rep(tmp,knots)
     rm(tmp)
     w<-matrix(w,rT,knots)
     w<-t(w); w<-c(w);
     w0 <- rep(0,knots*r)
    # 
    if((length(x.names.sp) == 0) & (length(x.names.tp) == 0)){
      # non-spatial and non-temporal beta
      # check for T
        if(r > 1){ 
         if(length(T) != r){         
           T<-rep(T,r) 
         }
        }
        #  
        out <- NULL
        out <- .C("GIBBS_zfitsum_onephi_gpp", as.integer(cov),
            as.integer(spdecay), as.double(flag), as.integer(nItr), 
            as.integer(nBurn), as.integer(n), as.integer(knots),
            as.integer(T), as.integer(r), as.integer(rT), 
            as.integer(p), as.integer(N), as.integer(report), 
            as.double(shape_e), as.double(shape_eta),as.double(shape_l), 
			as.double(phi_a), as.double(phi_b),
			as.double(priors$prior_a),as.double(priors$prior_b), as.double(priors$mu_beta),
            as.double(priors$delta2_beta), as.double(priors$mu_rho), 
            as.double(priors$delta2_rho), as.double(priors$alpha_l), 
            as.double(priors$delta2_l), as.double(init.phi),
            as.double(tuning), as.double(phis), as.integer(phik),
            as.double(coords.D.knots), as.double(coords.D.obs.knots), 
            as.integer(1), 
            as.double(initials$sig2eps), as.double(initials$sig2eta),
            as.double(initials$sig2l), as.double(initials$beta),
            as.double(initials$rho), as.double(initials$mu_l),
            as.double(X), as.double(zm), as.double(w0), as.double(w),
            as.integer(trans),
            phip = double(nItr), accept = double(1), nup = double(nItr),
            sig2eps = double(nItr), sig2etap = double(nItr), 
            betap = matrix(double(p*nItr), p, nItr), rhop = double(nItr),
            mu_lp = matrix(double(r*nItr), r, nItr), 
            sig2lp = matrix(double(r*nItr),r, nItr), 
            w0p = matrix(double(knots * r * nItr),knots * r, nItr), 
            wp = matrix(double(knots*rT*nItr), knots * rT, nItr), 
            gof = as.double(1), penalty = as.double(1),
            fitted = matrix(double(2*n*rT),n*rT,2))[45:58]

    } # end for fixed beta
    #  
    else if((length(x.names.sp) > 0) & (length(x.names.tp) == 0)){
      # for spatial beta
      if(length(T)> 1){ stop("\n ## Error: Unequal T is currently not possible for spatially varying model. ##\n")}
      q <- length(x.names.sp)          # number of spatial covariates
      shape_beta<-(knots*q)/2+priors$prior_a
      sig2beta<-0.1
      betas<-rep(0,knots)
      if("(Intercept)" %in% x.names){
         if(sum(c(X[,1])) == 0){
         intercept <- 0
         }     
         else{ intercept <- 1 }
      }
      else{ 
         intercept <- 1 
      }
      #
        out <- .C("GIBBSsp_zfitsum_onephi_gpp", as.integer(intercept),
            as.integer(cov),
            as.integer(spdecay), as.double(flag), as.integer(nItr), 
            as.integer(nBurn), as.integer(n), as.integer(knots),
            as.integer(T), as.integer(r), as.integer(rT), 
            as.integer(p), as.integer(q), as.integer(N), as.integer(report), 
            as.double(shape_e), as.double(shape_eta), as.double(shape_beta),
            as.double(shape_l), as.double(priors$prior_a), 
            as.double(priors$prior_b), as.double(priors$mu_beta),
            as.double(priors$delta2_beta), as.double(priors$mu_rho), 
            as.double(priors$delta2_rho), as.double(priors$alpha_l), 
            as.double(priors$delta2_l), as.double(init.phi),
            as.double(tuning), as.double(phis), as.integer(phik),
            as.double(coords.D.knots), as.double(coords.D.obs.knots), 
            as.integer(1), 
            as.double(initials$sig2eps), as.double(initials$sig2eta),
            as.double(sig2beta),
            as.double(initials$sig2l), as.double(initials$beta), as.double(betas),
            as.double(initials$rho), as.double(initials$mu_l),
            as.double(X), as.double(Xsp), as.double(zm), as.double(w0), as.double(w),
            as.integer(trans),
            phip = double(nItr), accept = double(1), nup = double(nItr),
            sig2eps = double(nItr), sig2etap = double(nItr), sig2betap = double(nItr),
            betap = matrix(double(p*nItr), p, nItr), betasp = matrix(double(knots*q*nItr), knots*q, nItr), 
            rhop = double(nItr),
            mu_lp = matrix(double(r*nItr), r, nItr), 
            sig2lp = matrix(double(r*nItr),r, nItr), 
            w0p = matrix(double(knots * r * nItr),knots * r, nItr), 
            wp = matrix(double(knots*rT*nItr), knots * r * T, nItr), 
            gof = as.double(1), penalty = as.double(1),
            fitted = matrix(double(2*n*rT),n*rT,2))[49:64]
    } # end for spatial beta
    else{
         stop("\n#\n## Error: Dynamic temporally varying models are not available for the GPP based models. \n#")
    } # end
    #
    ###
    #
     if(X.out==TRUE){
        if((length(x.names.sp) == 0) & (length(x.names.tp) == 0)){
        # non-spatial and temporal beta
          out$X <- X
        }
        else if((length(x.names.sp) > 0) & (length(x.names.tp) == 0)){
        # for spatial beta
          out$X <- X
          out$Xsp <- Xsp
          dimnames(out$X)[[2]] <- as.list(x.names)
        }  
        else{
          stop("\n#\n## Error: \n#")
        }
     }
     if(Y.out==TRUE){
        out$Y <- Y
     }
    #
      out$accept <- round(out$accept/nItr*100,2)
      out$call<-formula
     #
      cat("##","\n")
      cat("# nBurn = ",nBurn,". Iterations = ",nItr,".", "\n")
      cat("# Acceptance rate: (phi) = ",out$accept,"%", "\n")
      cat("##","\n")
     #
           out$phip <- as.matrix(out$phip[(nBurn+1):nItr])
           if(cov==4){
           out$nup <- as.matrix(out$nup[(nBurn+1):nItr])
           }
           out$sig2eps <- as.matrix(out$sig2eps[(nBurn+1):nItr])
           out$sig2etap <- as.matrix(out$sig2etap[(nBurn+1):nItr])
           if((length(x.names.sp) > 0)){
           out$sig2betap <- as.matrix(out$sig2betap[(nBurn+1):nItr])
           } 
           out$sig2lp <- matrix(out$sig2lp[1:r,(nBurn+1):nItr],r,length((nBurn+1):nItr))
           out$rhop <- as.matrix(out$rhop[(nBurn+1):nItr])
           out$betap <- matrix(out$betap[1:p,(nBurn+1):nItr],p,length((nBurn+1):nItr))
           if(length(x.names.sp) != 0){          
           out$betasp <- matrix(out$betasp[1:(knots*q),(nBurn+1):nItr],knots*q,length((nBurn+1):nItr))
           }
           out$mu_lp <- matrix(out$mu_lp[1:r,(nBurn+1):nItr],r,length((nBurn+1):nItr))	
           out$w0p <- out$w0p[1:(knots*r),(nBurn+1):nItr]	
           out$wp <- out$wp[1:(knots*rT),(nBurn+1):nItr]
           dimnames(out$fitted)[[2]] <- c("Mean","SD")
           out$tol.dist<-tol.dist
           out$distance.method<-method
           out$cov.fnc<-cov.fnc
           out$scale.transform<-scale.transform
           out$sampling.sp.decay<-spatial.decay
           out$covariate.names<-c(x.names)
           if(length(x.names.sp) != 0){          
           out$sp.covariate.names<-c(x.names.sp)
           }
           out$Distance.matrix.knots<-coords.D.knots
           out$knots.coords<-knots.coords
           out$coords<-coords
           if(method=="geodetic:km"){
           out$KM <- TRUE
           }
           if(method=="geodetic:mile"){
           out$KM <- FALSE
           }
           out$n <- n
           out$r <- r
           out$T <- T
           out$p <- p
           if(length(x.names.sp) != 0){          
           out$q<-q
           }
           out$knots <- knots
           out$initials <- initials	
           out$priors <- priors	
           out$gof <- round(out$gof,2)
           out$penalty <- round(out$penalty,2)
           tmp <- matrix(c(out$gof,out$penalty,out$gof+out$penalty),1,3)
           dimnames(tmp)[[2]]<-c("Goodness.of.fit","Penalty","PMCC")
           dimnames(tmp)[[1]]<-c("values:")
           out$PMCC <- tmp
           tmp <- NULL
           out$gof <- NULL
           out$penalty <- NULL
     #
          out$iterations <- nItr	
          out$nBurn <- nBurn	
     #
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   out$computation.time<-comp.time
    #
     #
     #class(out) <- "spGPP"
     #
     out
}
##
## Prediction of Z_lt for the GPP models using MCMC samples
##
spGPP.prediction<-function(nBurn, pred.data, pred.coords, 
               posteriors, tol.dist, Summary=TRUE, fitted.X=NULL)
{
    start.time<-proc.time()[3]
     #
        if (missing(posteriors)) {
           stop("Error: need to provide the posterior MCMC samples.")
        }
     #
          if (is.null(posteriors$X)) {
           stop("Error: need to provide the fitted covariate values.")
          }
        fitted.X <- posteriors$X
      # 
        knots.coords<-posteriors$knots.coords 	
        coords<-posteriors$coords 	
      #
        if (missing(knots.coords)) {
           stop("Error: need to specify the knots coords.")
        }
        if (!is.matrix(knots.coords)) {
           stop("Error: knots coords must be a (knots x 2) matrix of xy-coordinate locations.")
        }
        if ( (!is.numeric(knots.coords[,1])) | (!is.numeric(knots.coords[,2]))) {
           stop("\n Error: knots coords columns should be numeric \n")
        }
        if (missing(coords)) {
            stop("Error: need to specify the coords")
        }
        if ( !is.matrix(coords) ) {
            stop("Error: coords must be a (n x 2) matrix of xy-coordinate locations")
        }
        if ( (!is.numeric(coords[,1])) | (!is.numeric(coords[,2]))) {
            stop("\n Error: coords columns should be numeric \n")
        }
        if (missing(pred.coords)) {
           stop("Error: need to specify the prediction coords.")
        }
        if (!is.matrix(pred.coords)) {
           stop("Error: prediction coords must be a (n.pred.site x 2) matrix of xy-coordinate locations.")
        }
        if ( (!is.numeric(pred.coords[,1])) | (!is.numeric(pred.coords[,2]))) {
           stop("\n Error: prediction coords columns should be numeric \n")
        }
      #
           coords.all <- rbind(knots.coords,pred.coords)
           tn.knotsites <- length(knots.coords[, 1])
           nknot.sites <- 1:tn.knotsites
           tn.predsites <- length(coords.all[, 1]) - tn.knotsites
           npred.sites <- (tn.knotsites + 1):(length(coords.all[, 1]))
      #    
           method <- posteriors$distance.method   
           spT.check.sites.inside(coords.all, method, tol=tol.dist)
      #                                                           
        if(method == "geodetic:km"){
           coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=TRUE))
        }
        else if(method == "geodetic:mile"){
           coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=FALSE))
        }
        else {
           coords.D <- as.matrix(dist(coords.all, method, diag=TRUE, upper=TRUE))
        }
      #
      #
           dmns <- coords.D[nknot.sites, npred.sites] # m x ns
           dnsm <- t(dmns) # ns x m
           dm <- coords.D[1:tn.knotsites, 1:tn.knotsites] # m x m
      #
           all.coords <- rbind(coords,knots.coords)
           spT.check.sites.inside(all.coords, method, tol=tol.dist)
      #
        if(method == "geodetic:km"){
           coords.D.all <- as.matrix(spT.geodist(Lon=all.coords[,1],Lat=all.coords[,2], KM=TRUE))
        }
        else if(method == "geodetic:mile"){
           coords.D.all <- as.matrix(spT.geodist(Lon=all.coords[,1],Lat=all.coords[,2], KM=FALSE))
        }
        else {
           coords.D.all <- as.matrix(dist(all.coords, method, diag=TRUE, upper=TRUE))
        }
      #
           dnm <- coords.D.all[1:length(coords[,1]),(length(coords[,1])+1):length(all.coords[,1])]
      #
           al.coords <- rbind(pred.coords,coords)
           spT.check.sites.inside(al.coords, method, tol=tol.dist)
        if(method == "geodetic:km"){
           coords.D.al <- as.matrix(spT.geodist(Lon=al.coords[,1],Lat=al.coords[,2], KM=TRUE))
           dnsn <- coords.D.al[1:length(pred.coords[,1]),(length(pred.coords[,1])+1):length(al.coords[,1])]
        }
        else if(method == "geodetic:mile"){
           coords.D.al <- as.matrix(spT.geodist(Lon=al.coords[,1],Lat=al.coords[,2], KM=FALSE))
           dnsn <- coords.D.al[1:length(pred.coords[,1]),(length(pred.coords[,1])+1):length(al.coords[,1])]
        }
        else {
           coords.D.al <- as.matrix(dist(al.coords, method, diag=TRUE, upper=TRUE))
           dnsn <- coords.D.al[1:length(pred.coords[,1]),(length(pred.coords[,1])+1):length(al.coords[,1])]
        }
      #
           nsite <- tn.predsites
      #
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
           nItr <- posteriors$iterations
           if((nBurn+posteriors$nBurn) >= nItr){
              stop(paste("\n Error: burn-in >= iterations\n Here, nBurn = ",nBurn+posteriors$nBurn," and iterations = ",nItr,"."))
           }
           nItr <-(nItr-posteriors$nBurn) 
           itt <-(nItr-nBurn)
       # 
           m <- tn.knotsites
           r <- posteriors$r
           T <- posteriors$T
         # checking unequal T
         if(length(T) > 1){
           rT <- sum(T)
         }
         else{
           rT <- r*T
         }
           p <- length(c(posteriors$betap))/nItr
           n <- posteriors$n
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
      pred.xsp<-Formula.matrix(call.f,data=pred.data)[[4]]
      pred.xtp<-Formula.matrix(call.f,data=pred.data)[[6]]
    #
    #
           phip<-posteriors$phip[(nBurn+1):nItr,]
           if(cov==4){
           nup<-posteriors$nup[(nBurn+1):nItr,]
           }
           else{
           nup<-0
           } 
           sig2ep<-posteriors$sig2ep[(nBurn+1):nItr,]
           sig2etap<-posteriors$sig2etap[(nBurn+1):nItr,]
           sig2lp<-matrix(posteriors$sig2lp,r,nItr)
           sig2lp<-sig2lp[,(nBurn+1):nItr]
           rhop<-posteriors$rhop[(nBurn+1):nItr,]
           mu_lp<-matrix(posteriors$mu_lp,r,nItr)
           mu_lp<-mu_lp[,(nBurn+1):nItr]
           betap<-matrix(posteriors$betap,p,nItr)
           betap<-betap[,(nBurn+1):nItr]
           w0p<-matrix(posteriors$w0p,m*r,nItr)
           w0p<-w0p[,(nBurn+1):nItr]	
           wp<-matrix(posteriors$wp,m*rT,nItr)
           wp<-wp[,(nBurn+1):nItr]	
      #
          #
      if(posteriors$scale.transform=="NONE"){ 
        scale.transform <- 1
      }
      else if(posteriors$scale.transform=="SQRT"){ 
        scale.transform <- 2
      }
      else if(posteriors$scale.transform=="LOG"){ 
        scale.transform <- 3
      }
      else {
         stop("scale.transform is not correctly defined.")
      }
      #
      if((is.null(pred.xsp)) & (is.null(pred.xtp))){
      # fixed beta
      out<-matrix(.C('z_pr_its_gpp1', as.integer(cov), 
           as.integer(scale.transform), as.integer(itt),
           as.integer(nsite),as.integer(n),as.integer(m),
           as.integer(r),as.integer(T),as.integer(rT),as.integer(p), 
           as.integer(nsite*rT),as.double(phip),as.double(nup),as.double(dm),
           as.double(dnsm),as.double(wp),as.double(sig2ep),
           as.double(betap),as.double(pred.x),as.integer(1),
           zpred=double(itt*nsite*rT))$zpred,rT*nsite,itt)
      #
      }
      else if((!is.null(pred.xsp)) & (is.null(pred.xtp))){
      # spatial beta
        sig2betap<-posteriors$sig2betap[(nBurn+1):nItr,]
        betasp<-posteriors$betasp[,(nBurn+1):nItr]
        q<-posteriors$q 
        out<-.C("z_pr_its_gpp1_sp",as.integer(cov),as.integer(scale.transform),
             as.integer(itt),as.integer(nsite),as.integer(n),as.integer(m),
             as.integer(r),as.integer(T),as.integer(rT),as.integer(p),as.integer(q),
             as.integer(nsite*rT),as.double(phip),as.double(nup),as.double(dm),
             as.double(dnsm),as.double(wp),as.double(sig2ep),as.double(sig2betap),
             as.double(betap),as.double(betasp),as.double(pred.x),as.double(pred.xsp),
             as.integer(1),betapred=double(q*nsite*itt),#matrix(double(q*nsite*itt),q*nsite,itt), #array(double(q*nsite*itt),dim=c(q,nsite,itt)),
             zpred=matrix(double(rT*nsite*itt),rT*nsite,itt))[25:26]
      }
      #
      else if((is.null(pred.xsp)) & (!is.null(pred.xtp))){
      # temporal beta
         stop("Error: dynamic modellings is not currently available for the GPP based models")
      }
      else if((!is.null(pred.xsp)) & (!is.null(pred.xtp))){
      # both spatial and temporal beta
         stop("Error: spatio-dynamic modellings is not currently available for the GPP based models")
      }
      else{
         stop("Error: correctly define X variables")
      } 
      #
      output<-NULL
      if((!is.null(pred.xsp)) & (is.null(pred.xtp))){
      # spatial beta
        output$pred.samples<-out$zpred
        output$pred.spbeta.samples<-array(c(out$betapred),dim=c(q,nsite,itt))
      }
      else{
        output$pred.samples<-out
      }
      #
      out<-NULL
      output$knots.coords <- posteriors$knots.coords
      output$pred.coords <- pred.coords
      output$distance.method <- posteriors$distance.method
      output$cov.fnc <- posteriors$cov.fnc
      output$scale.transform <- posteriors$scale.transform
      output$n <- n
      output$r <- r
      output$T <- T
      output$knots <- m
      #
        if(Summary == TRUE){
         if(itt < 40){
          cat("##", "\n")
          cat("# Summary statistics are not given, because the number of predicted samples (",nItr,") are too small.\n# nBurn = ",nBurn+posteriors$nBurn,". Iterations = ",posteriors$iterations,".", "\n")
          cat("##", "\n")
    #
   end.time<-proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   output$computation.time<-comp.time
    #
          #class(output) <- "spGPP"
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
          szp <- NULL
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   output$computation.time<-comp.time
    #
          #class(output) <- "spGPP"
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
          #class(output) <- "spGPP"
          output
      }
}
##
## K-step forecasts for the GPP models 
##
 spGPP.forecast<-function(nBurn, K, fore.data=NULL, 
       fore.coords=NULL, posteriors, tol.dist, Summary=TRUE)
{
      start.time<-proc.time()[3]
      #
      if (missing(posteriors)) {
           stop("Error: need to provide the posterior MCMC samples.")
      }
      #
           nItr <- posteriors$iterations
           if((nBurn+posteriors$nBurn) >= nItr){
              stop(paste("\n Error: burn-in >= iterations\n Here, nBurn = ",nBurn+posteriors$nBurn," and iterations = ",nItr,"."))
           }
           nItr <-(nItr-posteriors$nBurn) 
           itt <-(nItr-nBurn)
      #
           scale.transform<-posteriors$scale.transform
           coords <- posteriors$coords
           knots.coords <- posteriors$knots.coords
      #
       if(is.null(fore.coords)){
              stop("Error: need to provide fore.coords ")
       }
       if(missing(fore.coords)) {
         stop("Error: need to specify the fore.coords")
       }
       if ( !is.matrix(fore.coords) ) {
         stop("Error: fore.coords must be a matrix of xy-coordinate locations")
       }
       if ( (!is.numeric(fore.coords[,1])) | (!is.numeric(fore.coords[,2]))) {
         stop("\n Error: fore.coords columns should be numeric \n")
       }
      #
           nsite <- dim(fore.coords)[[1]]
           n <- posteriors$n
           r <- posteriors$r
           T <- posteriors$T
           m <- posteriors$knots  
           if(length(T)>1){
             rT <- sum(T)
           } 
           else{
             rT <- r*T
           }
           N <- n*rT
           p <- posteriors$p
      #
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
    # adding the formula in to the forecast dataset
      if(!is.data.frame(fore.data) & !is.null(fore.data)){
        stop("#\n# Error: fore.data should be in data format\n#")
      }
      call.f<-posteriors$call  
      call.f<-as.formula(paste("tmp~",paste(call.f,sep="")[[3]]))
      if(is.data.frame(fore.data)){
      if((nsite*r*K)!=dim(fore.data)[[1]]){
        stop("\n# Check the fore.data and/or fore.coords and/or spT.time function#\n")
        #print("#\n # Check the fore.data \n#\n")
      }
      fore.data$tmp<-rep(1,nsite*r*K)
      }
      if(is.null(fore.data)){
        fore.data<-data.frame(tmp=rep(1,nsite*r*K))
      }
      fore.x<-Formula.matrix(call.f,data=fore.data)[[2]]
    #
    #
      #
           dm <- posteriors$Distance.matrix.knots
      #
           coords.all <- rbind(knots.coords,fore.coords)
           tn.knotsites <- length(knots.coords[, 1])
           nknot.sites <- 1:tn.knotsites
           tn.predsites <- length(coords.all[, 1]) - tn.knotsites
           npred.sites <- (tn.knotsites + 1):(length(coords.all[, 1]))
      #    
           method <- posteriors$distance.method   
           spT.check.sites.inside(coords.all, method, tol=tol.dist)
      #                                                           
        if(method == "geodetic:km"){
           coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=TRUE))
        }
        else if(method == "geodetic:mile"){
           coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=FALSE))
        }
        else {
           coords.D <- as.matrix(dist(coords.all, method, diag=TRUE, upper=TRUE))
        }
      #
      #
           dmns <- coords.D[nknot.sites, npred.sites] # m x ns
           dnsitem <- t(dmns) # ns x m
           dmns <- NULL
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
           wp<-posteriors$wp[,(nBurn+1):nItr]
      #
      out<-NULL
      out<-matrix(.C("zlt_fore_gpp_its",as.integer(cov),
       as.integer(itt),as.integer(K),as.integer(nsite),
       as.integer(m),as.integer(r),as.integer(p),as.integer(rT),
       as.integer(T),as.integer(r*K),as.integer(nsite*r*K),
       as.double(dnsitem),as.double(dm),as.double(phip),as.double(nup),
       as.double(sig_ep),as.double(sig_etap),as.double(betap),
       as.double(rhop),as.double(wp),as.double(fore.x),
       as.integer(1),fz=double(nsite*r*K*itt))$fz,nsite*r*K,itt)
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
      output$knots.coords <- knots.coords
      output$fore.coords <- fore.coords
      output$distance.method<-posteriors$distance.method  
      output$cov.fnc<-posteriors$cov.fnc  
      #output$scale.transform <- posteriors$scale.transform

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
        out <- NULL
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   output$computation.time<-comp.time
    #
          #class(output) <- "spGPP"
          output
         }
         else { 
          cat("##", "\n")
          cat("# Forecast samples and summary statistics are given.\n# nBurn = ",nBurn,". Iterations = ",nItr,".", "\n")
          cat("##", "\n")
          #
          szp<-spT.Summary.Stat(output$fore.samples[,])
          # 
          out <- NULL
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
          #class(output) <- "spGPP"
          output
         }
      }
      else {
          cat("##", "\n")
          cat("# Forecast samples are given.\n# nBurn = ",nBurn,", Iterations = ",itt,".", "\n")
          cat("##", "\n")
      #
       out <- NULL
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   output$computation.time<-comp.time
    #
          #class(output) <- "spGPP"
          output
      }
}
##
## up to this point
##

##
## MCMC fit and predictions
##
spGPP.MCMC.Pred<-function(formula, data=parent.frame(), pred.data,
             time.data, knots.coords, coords, pred.coords, priors, 
             initials, nItr, nBurn=0, report=1, tol.dist=2,
             distance.method="geodetic:km", cov.fnc="exponential", 
             scale.transform="NONE", spatial.decay, annual.aggregation="NONE")
{
   #
    start.time<-proc.time()[3]
   #
    if (nBurn >= nItr) {
         stop(paste("\n Error: iterations =< nBurn\n Here, nBurn = ",nBurn," and iterations = ",nItr,"."))
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
         Xtp <- XY[[6]]
         x.names.tp <- XY[[7]]
   #
   if(length(x.names.sp)>0){ stop("\n Error: currently not available for spatially varying model. \n")} 
   if(length(x.names.tp)>0){ stop("\n Error: currently not available for temporally varying model. \n")} 
   #
    if (missing(coords)) {
         stop("Error: need to specify the coords")
    }
    if ( !is.matrix(coords) ) {
         stop("Error: coords must be a (n x 2) matrix of xy-coordinate locations")
    }
    if ( (!is.numeric(coords[,1])) | (!is.numeric(coords[,2]))) {
         stop("\n Error: coords columns should be numeric \n")
    }
   #
    if (missing(knots.coords)) {
         stop("Error: need to specify the knots.coords")
    }
    if ( !is.matrix(knots.coords) ) {
         stop("Error: knots.coords must be a (n x 2) matrix of xy-coordinate locations")
    }
    if ( (!is.numeric(knots.coords[,1])) | (!is.numeric(knots.coords[,2]))) {
         stop("\n Error: knots.coords columns should be numeric \n")
    }
   #
   #
    if (missing(pred.coords)) {
         stop("Error: need to specify the pred.coords")
    }
    if ( !is.matrix(pred.coords) ) {
         stop("Error: pred.coords must be a (pred.site x 2) matrix of xy-coordinate locations")
    }
    if ( (!is.numeric(pred.coords[,1])) | (!is.numeric(pred.coords[,2]))) {
         stop("\n Error: pred.coords columns should be numeric \n")
    }
   #
   #
    if (missing(pred.data)) {
         stop("Error: need to specify the pred.data")
    }
    if ( !is.data.frame(pred.data) ) {
         stop("Error: pred.data must be a data frame")
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
         knots <- length(knots.coords[,1]) 
   #
    if(n <= knots){
         stop("Error: n must not be smaller than or equal to knots")
    }
   #
    #
      method <- distance.method
      spT.check.sites.inside(knots.coords, method, tol=tol.dist)
   #
         all.coords <- rbind(coords,knots.coords)
   #
      if(method=="geodetic:km"){
         coords.D.knots <- as.matrix(spT.geodist(Lon=knots.coords[,1],Lat=knots.coords[,2], KM=TRUE))
         coords.D.all <- as.matrix(spT.geodist(Lon=all.coords[,1],Lat=all.coords[,2], KM=TRUE))
     }
     else if(method=="geodetic:mile"){
         coords.D.knots <- as.matrix(spT.geodist(Lon=knots.coords[,1],Lat=knots.coords[,2], KM=FALSE))
         coords.D.all <- as.matrix(spT.geodist(Lon=all.coords[,1],Lat=all.coords[,2], KM=FALSE))
     }
     else {
         coords.D.knots<- as.matrix(dist(knots.coords, method, diag = T, upper = T))
         coords.D.all <- as.matrix(dist(all.coords, method, diag = T, upper = T))
     }
   #
    if(knots != length(coords.D.knots[,1])){
         stop("Error: knots must be equal to knots.coords")
    }
    #
         coords.D.obs.knots <- coords.D.all[1:length(coords[,1]),(length(coords[,1])+1):length(all.coords[,1])]
    #
    #
    if (N != length(Y)) {
         stop("#\n# Error: Years, Months, and Days are misspecified,\n     i.e., total number of observations in the data set \n     should be equal to N\n  : N = n * r * T \n   where, N = total number of observations in the data,\n          n = total number of sites,\n          r = total number of years,\n          T = total number of days. \n# Check the function spT.time.\n#\n")
    }
    # 
         p <- length(x.names)          # number of covariates
    #
         priors <- priors.checking.gpp(priors,r,p)
    #
         shape_e <- N/2+priors$prior_a
         shape_eta <- (knots*rT)/2+priors$prior_a
         shape_l <- (knots*r)/2+priors$prior_a
    #
         zm <- matrix(Y,rT,n)
         zm <- apply(zm,1,median,na.rm=TRUE)
         zm <- rep(zm,n)
         zm <- cbind(Y,zm)
         zm[is.na(zm[,1]),1]=zm[is.na(zm[,1]),2]
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
         zm <- zm[,1]
         trans <- 0
         scale.transform = "NONE"
    }
    else if (scale.transform == "SQRT") {
         zm <- sqrt(zm[,1])
         trans <- 1
         scale.transform = "SQRT"
    }
    else if (scale.transform == "LOG") {
         zm <- log(zm[,1])
         trans <- 2
         scale.transform = "LOG"
    }
    else{
         stop("\n Error: scale.transform is not correctly specified \n")
    }
    #
    #
      initials<-initials.checking.gpp(initials,zm,X,n,r,T,coords.D.knots)
    #
    #
    if(spatial.decay$type=="FIXED"){
         spdecay <- 1
		 if(is.null(spatial.decay$value)){
		 spatial.decay$value <- (3/max(c(coords.D.knots))) 
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
    if (length(initials$mu_l) != r){
         stop("Error: need to specify correct number of years (r) for mu_l initials.")
    }
    if (length(initials$sig2l) != r){
         stop("Error: need to specify correct number of years (r) for sig2l initials.")
    }
    if (length(initials$beta) != p){
         stop("Error: need to specify correct number (p) of initial parameters for beta.")
    }
    #
    #
    if (length(priors$mu_beta) != p){
         stop("Error: need to specify correct number (p) of mu_beta priors.")
    }
    if (length(priors$delta2_beta) != p*p){
         stop("Error: need to specify correct number (p) of delta2_beta priors.")
    }
    if (length(priors$alpha_l) != r){
         stop("Error: need to specify correct number (r) of alpha_l priors.")
    }
    if (length(priors$delta2_l) != r){
         stop("Error: need to specify correct number (r) of delta2_l priors.")
    }
    #
    if (missing(nBurn)) {
         stop("Error: nBurn must be specified")
    }
    #
     w <- rep(0,rT*knots)
     w0 <- rep(0,knots*r)
    #
    # prediction part
    #
    #  	
         coords.all <- rbind(knots.coords,pred.coords)
         tn.knotsites <- length(knots.coords[, 1])
         nknot.sites <- 1:tn.knotsites
         tn.predsites <- length(coords.all[, 1]) - tn.knotsites
         npred.sites <- (tn.knotsites + 1):(length(coords.all[, 1]))
  #                                                                         
      if(method=="geodetic:km"){
         coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=TRUE))
      }
      else if(method=="geodetic:mile"){
         coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=FALSE))
      }
      else {
           coords.D<- as.matrix(dist(coords.all, method, diag = T, upper = T))
      }
    #
         dmns <- coords.D[nknot.sites, npred.sites] # m x ns
         dnsm <- t(dmns) # ns x m
    #
         nsite <- tn.predsites
         m <- tn.knotsites
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
        out <- NULL
        out <- .C("GIBBS_sumpred_gpp",as.integer(aggtype),as.integer(cov), as.integer(spdecay), 
            as.double(flag), as.integer(nItr),
            as.integer(nBurn), as.integer(n), as.integer(m), as.integer(T),
            as.integer(r), as.integer(rT), as.integer(p), as.integer(N), as.integer(report),
            as.double(shape_e), as.double(shape_eta), as.double(shape_l),
            as.double(phi_a), as.double(phi_b),			
            as.double(priors$prior_a), as.double(priors$prior_b), as.double(priors$mu_beta),
            as.double(priors$delta2_beta), as.double(priors$mu_rho), as.double(priors$delta2_rho),
            as.double(priors$alpha_l), as.double(priors$delta2_l), as.double(init.phi),
            as.double(tuning), as.double(phis), as.integer(phik),
            as.double(coords.D.knots), as.double(coords.D.obs.knots),
            as.integer(1), as.double(initials$sig2eps), as.double(initials$sig2eta),
            as.double(initials$sig2l), as.double(initials$beta), as.double(initials$rho), 
            as.double(initials$mu_l), as.double(X), as.double(zm), as.double(w0),
            as.double(w), as.integer(nsite), as.integer(nsite*rT), as.double(dnsm), 
            as.double(pred.x), as.integer(trans), accept=double(1), 
            gof=double(1),penalty=double(1))[50:52]    
     #
      out$accept <- round(out$accept/nItr*100,2)
      out$call<-formula
     #
      cat("##","\n")
      cat("# nBurn = ",nBurn,". Iterations = ",nItr,".", "\n")
      cat("# Acceptance rate: (phi) = ",out$accept,"%", "\n")
      cat("##","\n")
      cat("# Text output is given: \n##\n")
     #
           tmp<-read.table('OutGPP_Values_Parameter.txt', sep='',header=FALSE)
           tmp<-tmp[(nBurn+1):nItr,]
           tmp<-spT.Summary.Stat(t(tmp))
           if(cov==4){
           row.names(tmp)<- c(x.names,"rho","sig2eps","sig2eta","phi","nu")
           }
           else{
           row.names(tmp)<- c(x.names,"rho","sig2eps","sig2eta","phi")
           }
           out$parameters<-round(tmp, 4)
           tmp<-NULL

           tmp<-read.table('OutGPP_Stats_FittedValue.txt', sep=',',header=FALSE)
           names(tmp)<- c("Mean","SD")
           out$fitted<-round(tmp, 4)
           tmp<-NULL

           tmp<-read.table('OutGPP_Stats_PredValue.txt', sep=',',header=FALSE)
           names(tmp)<- c("Mean","SD")
           out$prediction<-round(tmp, 4)
           tmp<-NULL
#
       if(annual.aggregation=="ave"){
           tmp<-read.table('OutGPP_Annual_Average_Prediction.txt', sep='',header=FALSE)
           tmp<-spT.Summary.Stat(t(tmp))
           out$an.agr.pred.average<-round(tmp, 4)
           tmp<-NULL
       }   
       if(annual.aggregation=="an4th"){
           tmp<-read.table('OutGPP_Annual_4th_Highest_Prediction.txt', sep='',header=FALSE)
           tmp<-spT.Summary.Stat(t(tmp))
           out$an.agr.pred.an4th<-round(tmp, 4)
           tmp<-NULL
       }   
       if(annual.aggregation=="w126"){
           tmp<-read.table('OutGPP_Annual_w126_Prediction.txt', sep='',header=FALSE)
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
           out$knots <- knots
           out$pred.n<-nsite
           out$r <- r
           out$T <- T
           out$T <- Y
           out$initials <- initials	
           out$priors <- priors 
           out$iterations <- nItr 
           out$nBurn <- nBurn 
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   out$computation.time<-comp.time
    #

     #
     #class(out) <- "spGPP"
     #
     out
}
##
##
##