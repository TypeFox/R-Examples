##
## initials checking for GP
##
 initials.checking.gp<-function(x,o,X,Xsp,Xtp,n,r,T,d){
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
     if((is.null(Xsp)) & (is.null(Xtp))){
      # 
      ## for non spatial and temporal beta
      if(is.null(x)){
       x$phi=NULL; x$sig2eps=NULL; x$sig2eta=NULL; x$beta=NULL;
       x$phi<-(3/max(c(d)))#(-log(0.05)/max(c(d)))
       x$sig2eps <- 0.01; x$sig2eta <- 0.1
       lm.coef<-lm(c(o) ~ X-1)$coef
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
       if(is.null(x$beta)){
        lm.coef<-lm(c(o) ~ X -1)$coef
        lm.coef[is.na(lm.coef)]<-0
        x$beta<-lm.coef[1:(dim(X)[[2]])]
       }
      x
      }
     }
     else if((!is.null(Xsp)) & (is.null(Xtp))){
     #
     ## for spatial beta
      if(is.null(x)){
       x$phi=NULL; x$sig2eps=NULL; x$sig2eta=NULL; x$beta=NULL;
       x$phi<-(3/max(c(d)))#(-log(0.05)/max(c(d)))
       x$sig2eps <- 0.01; x$sig2eta <- 0.1; x$sig2beta <- 0.1
       q<-length(c(Xsp))/(n*rT)
       dump<-sort(rep(1:n,rT)) 
       dump<-model.matrix(~factor(dump)-1)
       tdump<-NULL; for(i in 1:q){ tdump <- cbind(tdump,dump) }
       dump<-array(c(Xsp),dim=c(rT,n,q))
       for(j in 1:q){for(i in 1:n){ tdump[,i+(j-1)*n]<-tdump[,i+(j-1)*n]*dump[,i,j] }}
       lm.coef<-lm(c(o) ~ X+tdump-1)$coef
       lm.coef[is.na(lm.coef)]<-0
       x$beta<-lm.coef[1:(dim(X)[[2]])]
       x$betasp<-lm.coef[(dim(X)[[2]]+1):((dim(X)[[2]])+(n*q))]
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
       if(is.null(x$sig2beta)){
        x$sig2beta <- 0.1
       }
       if((is.null(x$beta)) | (is.null(x$betasp))){
        q<-length(c(Xsp))/(n*rT)
        dump<-sort(rep(1:n,rT)) 
        dump<-model.matrix(~factor(dump)-1)
        tdump<-NULL; for(i in 1:q){ tdump <- cbind(tdump,dump) }
        dump<-array(c(Xsp),dim=c(rT,n,q))
        for(j in 1:q){for(i in 1:n){ tdump[,i+(j-1)*n]<-tdump[,i+(j-1)*n]*dump[,i,j] }}
        lm.coef<-lm(c(o) ~ X+tdump-1)$coef
        lm.coef[is.na(lm.coef)]<-0
        x$beta<-lm.coef[1:(dim(X)[[2]])]
        x$betasp<-lm.coef[(dim(X)[[2]]+1):((dim(X)[[2]])+(n*q))]
       }
       x
      }
     }
     else if((is.null(Xsp)) & (!is.null(Xtp))){
     #
     ## for temporal beta
      if(is.null(x)){
       x$phi=NULL; x$sig2eps=NULL; x$sig2eta=NULL; x$beta=NULL; x$rhotp <- 1;
       x$phi<-(3/max(c(d)))#(-log(0.05)/max(c(d)))
       x$sig2eps <- 0.01; x$sig2eta <- 0.1; x$sig2delta <- 0.1
       u<-length(c(Xtp))/(n*rT)
       dump<-rep(1:T,n*r) 
       dump<-model.matrix(~factor(dump)-1)
       tdump<-NULL; for(i in 1:u){ tdump <- cbind(tdump,dump) }
       dump<-array(c(Xtp),dim=c(rT*n,u))
       for(j in 1:u){ for(i in 1:T){ tdump[,i+(j-1)*T]<-tdump[,i+(j-1)*T]*dump[,j]}}
       lm.coef<-lm(c(o) ~ X+tdump-1)$coef
       lm.coef[is.na(lm.coef)]<-0
       x$beta<-lm.coef[1:(dim(X)[[2]])]
       x$betatp<-t(matrix(lm.coef[(dim(X)[[2]]+1):((dim(X)[[2]])+(T*u))],T,u))
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
       if(is.null(x$sig2delta)){
        x$sig2delta <- 0.1
       }
       if(is.null(x$rhotp)){
        x$rhotp <- 1;
       }
        if((is.null(x$beta)) | (is.null(x$betatp))){
        u<-length(c(Xtp))/(n*rT)
        dump<-rep(1:T,n*r) 
        dump<-model.matrix(~factor(dump)-1)
        tdump<-NULL; for(i in 1:u){ tdump <- cbind(tdump,dump) }
        dump<-array(c(Xtp),dim=c(rT*n,u))
        for(j in 1:u){ for(i in 1:T){ tdump[,i+(j-1)*T]<-tdump[,i+(j-1)*T]*dump[,j]}}
        lm.coef<-lm(c(o) ~ X+tdump-1)$coef
        lm.coef[is.na(lm.coef)]<-0
        x$beta<-lm.coef[1:(dim(X)[[2]])]
        x$betatp<-t(matrix(lm.coef[(dim(X)[[2]]+1):((dim(X)[[2]])+(T*u))],T,u))
       }
       x
      }
     }
     else if((!is.null(Xsp)) & (!is.null(Xtp))){
     #
     ## for both spatial and temporal beta
      if(is.null(x)){
       x$phi=NULL; x$sig2eps=NULL; x$sig2eta=NULL; x$sig2beta=NULL; x$beta=NULL;
       x$phi<-(3/max(c(d)))#(-log(0.05)/max(c(d)))
       x$sig2eps <- 0.01; x$sig2eta <- 0.1; x$sig2beta <- 0.1; x$sig2delta <- 0.1; x$rhotp <- 1;

       q<-length(c(Xsp))/(n*rT)
       dump<-sort(rep(1:n,rT)) 
       dump<-model.matrix(~factor(dump)-1)
       tdump.sp<-NULL; for(i in 1:q){ tdump.sp <- cbind(tdump.sp,dump) }
       dump<-array(c(Xsp),dim=c(rT,n,q))
       for(j in 1:q){for(i in 1:n){ tdump.sp[,i+(j-1)*n]<-tdump.sp[,i+(j-1)*n]*dump[,i,j] }}

       u<-length(c(Xtp))/(n*rT)
       dump<-rep(1:T,n*r) 
       dump<-model.matrix(~factor(dump)-1)
       tdump.tp<-NULL; for(i in 1:u){ tdump.tp <- cbind(tdump.tp,dump) }
       dump<-array(c(Xtp),dim=c(rT*n,u))
       for(j in 1:u){ for(i in 1:T){ tdump.tp[,i+(j-1)*T]<-tdump.tp[,i+(j-1)*T]*dump[,j]}}

       lm.coef<-lm(c(o) ~ X+tdump.sp+tdump.tp-1)$coef
       lm.coef[is.na(lm.coef)]<-0
       x$beta<-lm.coef[1:(dim(X)[[2]])]
       x$betasp<-lm.coef[(dim(X)[[2]]+1):((dim(X)[[2]])+(n*q))]
       x$betatp<-t(matrix(lm.coef[(dim(X)[[2]]+(n*q)+1):((dim(X)[[2]])+(n*q)+(T*u))],T,u))
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
       if(is.null(x$sig2beta)){
        x$sig2beta <- 0.1
       }
       if(is.null(x$sig2delta)){
        x$sig2delta <- 0.1
       }
       if(is.null(x$rhotp)){
        x$rhotp <- 1;
       }
       if((is.null(x$beta)) | (is.null(x$betasp))){
        q<-length(c(Xsp))/(n*rT)
        dump<-sort(rep(1:n,rT)) 
        dump<-model.matrix(~factor(dump)-1)
        tdump.sp<-NULL; for(i in 1:q){ tdump.sp <- cbind(tdump.sp,dump) }
        dump<-array(c(Xsp),dim=c(rT,n,q))
        for(j in 1:q){for(i in 1:n){ tdump.sp[,i+(j-1)*n]<-tdump.sp[,i+(j-1)*n]*dump[,i,j] }}

        u<-length(c(Xtp))/(n*rT)
        dump<-rep(1:T,n*r) 
        dump<-model.matrix(~factor(dump)-1)
        tdump.tp<-NULL; for(i in 1:u){ tdump.tp <- cbind(tdump.tp,dump) }
        dump<-array(c(Xtp),dim=c(rT*n,u))
        for(j in 1:u){ for(i in 1:T){ tdump.tp[,i+(j-1)*T]<-tdump.tp[,i+(j-1)*T]*dump[,j]}}

        lm.coef<-lm(c(o) ~ X+tdump.sp+tdump.tp-1)$coef
        lm.coef[is.na(lm.coef)]<-0
        x$beta<-lm.coef[1:(dim(X)[[2]])]
        x$betasp<-lm.coef[(dim(X)[[2]]+1):((dim(X)[[2]])+(n*q))]
        x$betatp<-lm.coef[(dim(X)[[2]]+(n*q)+1):((dim(X)[[2]])+(n*q)+(T*u))]
       }
       x
      }
     }
     else{
       stop("Error:")
     }
     ##
}
##
## priors checking for GP
##
 priors.checking.gp<-function(x){
     #
     if(is.null(x)){
       x$prior_a<-NULL; x$prior_b<-NULL;
       x$prior_mubeta<-NULL; x$prior_sigbeta<-NULL 
       x$prior_omu<-NULL; x$prior_osig<-NULL 
       x$prior_a<-2; x$prior_b<-1;
       x$prior_mubeta<-0
       x$prior_sigbeta<-10^10
       x$prior_omu<-0
       x$prior_osig<-10^10
       x
     }
     else{
     if(is.null(x$prior_a)){
       x$prior_a<-2
     } 
     if(is.null(x$prior_b)){
       x$prior_b<-1
     } 
     if(is.null(x$prior_mubeta)){
       x$prior_mubeta<-0
     }
     if(is.null(x$prior_sigbeta)){
       x$prior_sigbeta<-10^10
     }
     if(is.null(x$prior_omu)){
       x$prior_omu<-0
     }
     if(is.null(x$prior_osig)){
       x$prior_osig<-10^10
     }
     x 
     }
}
##
## MCMC sampling for the GP models
## time.data format: col-1: year, col-2: month, col-3: day
##
spGP.Gibbs<-function(formula, data=parent.frame(), time.data, coords, 
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
    #if(is.wholenumber(nItr/report)==FALSE) {
    #     stop(paste("Error: ITERATION/",report," should be an integer number"))
    #}     
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
    if (missing(coords)) {
         stop("\n Error: need to specify the coords \n")
    }
    if ( !is.matrix(coords) & !is.data.frame(coords)) {
         stop("\n Error: coords must be a (n x 2) matrix or data frame of xy-coordinate locations \n")
    }
    if ( dim(coords)[[2]] !=2) {
         stop("\n Error: coords should have 2 columns \n")
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
   else{
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
         n <- length(coords[,1])              # number of sites
         r <- time.data[[1]]                  # number of years
         T <- time.data[[2]]                  # number of days
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
         stop(" Error: Years, Months, and Days are misspecified,\n i.e., total number of observations in the data set should be equal to N\n  : N = n * r * T \n   where, N = total number of observations in the data,\n          n = total number of sites,\n          r = total number of years,\n          T = total number of days. \n## Check spT.time function.")
    }
    #
         priors<-priors.checking.gp(priors)
    #
         shape_e<-N/2+priors$prior_a
         shape_eta<-N/2+priors$prior_a
    #
         zm <- matrix(Y,rT,n)
         zm <- apply(zm,1,median,na.rm=TRUE)
         zm <- rep(zm,n)
         zm <- cbind(Y,zm)
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
      initials<-initials.checking.gp(initials,zm[,1],X,Xsp,Xtp,n,r,T,coords.D)
    #
         o <- zm[,1]
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
         p <- length(x.names)          # number of covariates
    #
    if (length(initials$beta) != p){
         stop("Error: need to specify correct number of parameters for beta.")
    }
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
    out<-.C('GIBBS_gp',as.double(flag[,2]),as.integer(nItr),
         as.integer(nBurn),as.integer(n),as.integer(T),
         as.integer(r),as.integer(rT),as.integer(p),
         as.integer(N),as.integer(report),as.integer(cov),
         as.integer(spdecay),as.double(shape_e),
         as.double(shape_eta),
		 as.double(phi_a), as.double(phi_b),
		 as.double(priors$prior_a),
         as.double(priors$prior_b),as.double(priors$prior_mubeta),    
         as.double(priors$prior_sigbeta),as.double(priors$prior_omu), 
         as.double(priors$prior_osig),as.double(init.phi),
         as.double(tuning),as.double(phis),as.integer(phik),
         as.double(coords.D),
         as.double(initials$sig2eps), as.double(initials$sig2eta),
         as.double(initials$beta),as.double(X),as.double(zm[,1]),  
         as.double(o),as.integer(1),phip=double(nItr),
         accept=double(1),nup=double(nItr),sig2ep=double(nItr),
         sig2etap=double(nItr),betap=matrix(double(p*nItr),p,nItr), 
         op=matrix(double(N*nItr),N,nItr),fit=matrix(double(N*2),N,2), 
         gof=double(1),penalty=double(1))[35:44]
    }
    else if((length(x.names.sp) > 0) & (length(x.names.tp) == 0)){
    # for spatial beta
    stop("Error: Currently not available for this version.")
	if(length(T)> 1){ stop("\n ## Error: Unequal T is currently not possible for spatially varying model. ##\n")}
    q <- length(x.names.sp)          # number of spatial covariates
    shape_beta<-(n*q)/2+priors$prior_a
    if("(Intercept)" %in% x.names){
       if(sum(c(X[,1])) == 0){
       intercept <- 0
       }     
       else{ intercept <- 1 }
    }
    else{ 
       intercept <- 1 
    }
    out<-.C("GIBBSsp_gp",as.integer(intercept),as.double(flag[,2]),as.integer(nItr),as.integer(nBurn),
         as.integer(n),as.integer(T),as.integer(r),as.integer(rT),as.integer(p),
         as.integer(q),as.integer(N),as.integer(report),as.integer(cov),
         as.integer(spdecay),as.double(shape_e),as.double(shape_eta),as.double(shape_beta),   
         as.double(priors$prior_a),as.double(priors$prior_b),as.double(priors$prior_mubeta), 
         as.double(priors$prior_sigbeta),as.double(priors$prior_omu),as.double(priors$prior_osig),
         as.double(init.phi),as.double(tuning),as.double(phis),as.integer(phik),
         as.double(coords.D),as.double(initials$sig2eps),as.double(initials$sig2eta),as.double(initials$sig2beta), 
         as.double(initials$beta),as.double(initials$betas),as.double(X),as.double(Xsp),
         as.double(zm[,1]),as.double(o),as.integer(1), 
         phip=double(nItr), accept=double(1),nup=double(nItr),sig2ep=double(nItr),sig2etap=double(nItr),  
         sig2betasp=double(nItr),betap=matrix(double(nItr*p),p,nItr),betasp=matrix(double(nItr*n*q),n*q,nItr),
         op=matrix(double(N*nItr),N,nItr),fit=matrix(double(N*2),N,2),gof=double(1),penalty=double(1))[39:50]
    } 
    else if((length(x.names.sp) == 0) & (length(x.names.tp) > 0)){
    # for temporal beta
    stop("Error: Currently not available for this version.")
    if(length(T)> 1){ stop("\n ## Error: Unequal T is currently not possible for temporally varying model. ##\n")}
    u <- length(x.names.tp)          # number of temporal covariates
    shape_del<-(u*T)/2+priors$prior_a
    shape_0<-u/2+priors$prior_a
    if("(Intercept)" %in% x.names){
       if(sum(c(X[,1])) == 0){
       intercept <- 0
       }
       else{ intercept <- 1 }     
    }
    else{ 
       intercept <- 1 
    }
    if(is.null(initials$rhotp)){
      stop("Error: need to provide indication for rho sampling")
    }
    out<-.C("GIBBStp_gp",as.integer(intercept),as.double(flag[,2]),as.integer(nItr),as.integer(nBurn),
         as.integer(n),as.integer(T),as.integer(r),as.integer(rT),as.integer(p),
         as.integer(u),as.integer(N),as.integer(report),as.integer(cov),as.integer(spdecay),as.integer(initials$rhotp),
         as.double(shape_e),as.double(shape_eta),as.double(shape_del),as.double(shape_0),
         as.double(priors$prior_a),as.double(priors$prior_b), 
         as.double(priors$prior_mubeta),as.double(priors$prior_sigbeta),as.double(priors$prior_omu),
         as.double(priors$prior_osig),
         as.double(init.phi),as.double(tuning),as.double(phis),as.integer(phik),
         as.double(coords.D),as.double(initials$sig2eps),as.double(initials$sig2eta),as.double(initials$sig2delta),
         as.double(0.01),as.double(initials$beta),as.double(initials$betat),as.double(rep(1,u)),
         as.double(X),as.double(Xtp),as.double(zm[,1]),as.double(o),as.integer(1), 
         phip=double(nItr),accept=double(1),nup=double(nItr),sig2ep=double(nItr),sig2etap=double(nItr),  
         sig2deltap=double(nItr),sig2op=double(nItr),rhotp=matrix(double(u*nItr),u,nItr),
         betap=matrix(double(nItr*p),p,nItr),betat0p=matrix(double(u*nItr),u,nItr),betatp=matrix(double(nItr*u*T),u*T,nItr), 
         op=matrix(double(N*nItr),N,nItr),fit=matrix(double(N*2),N,2),gof=double(1),penalty=double(1))[43:57]
    }
    else if((length(x.names.sp) > 0) & (length(x.names.tp) > 0)){
    # for both spatial and temporal beta
    stop("Error: Currently not available for this version.")
    if(length(T)> 1){ stop("\n ## Error: Unequal T is currently not possible for spatially and temporally varying model.## \n")}
    q <- length(x.names.sp)          # number of spatial covariates
    u <- length(x.names.tp)          # number of temporal covariates
    shape_beta<-(n*q)/2+priors$prior_a
    shape_del<-(u*T)/2+priors$prior_a
    shape_0<-u/2+priors$prior_a
    if(is.null(initials$rhotp)){
      stop("Error: need to provide indication for rho sampling")
    }
    if("(Intercept)" %in% x.names){
       if(sum(c(X[,1])) == 0){
       intercept <- 0
       }
       else{ intercept <- 1 }     
    }
    else{ 
       intercept <- 1 
    }
    out<-.C("GIBBSsptp_gp",as.integer(intercept), as.double(flag[,2]),as.integer(nItr),as.integer(nBurn),
         as.integer(n),as.integer(T),as.integer(r),as.integer(rT),as.integer(p),
         as.integer(q),as.integer(u),as.integer(N),as.integer(report),as.integer(cov),as.integer(spdecay),as.integer(initials$rhotp),
         as.double(shape_e),as.double(shape_eta),as.double(shape_beta),as.double(shape_del),as.double(shape_0),
         as.double(priors$prior_a),as.double(priors$prior_b),
         as.double(priors$prior_mubeta),as.double(priors$prior_sigbeta),as.double(priors$prior_omu),as.double(priors$prior_osig), 
         as.double(init.phi),as.double(tuning),as.double(phis),as.integer(phik),
         as.double(coords.D),as.double(initials$sig2eps),as.double(initials$sig2eta),as.double(initials$sig2beta),
         as.double(initials$sig2delta),as.double(0.1),
         as.double(initials$beta),as.double(initials$betas),
         as.double(initials$betat),as.double(rep(1,u)),as.double(X),as.double(Xsp),as.double(Xtp), 
         as.double(zm[,1]),as.double(o),as.integer(1),
         phip=double(nItr),accept=double(1),nup=double(nItr),sig2ep=double(nItr),sig2etap=double(nItr),  
         sig2betasp=double(nItr),sig2deltap=double(nItr),sig2op=double(nItr),rhotp=matrix(double(u*nItr),u,nItr),
         betap=matrix(double(nItr*p),p,nItr),betasp=matrix(double(nItr*n*q),n*q,nItr),
         betat0p=matrix(double(u*nItr),u,nItr),betatp=matrix(double(nItr*u*T),u*T,nItr), 
         op=matrix(double(N*nItr),N,nItr),fit=matrix(double(N*2),N,2),gof=double(1),penalty=double(1))[48:64]
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
        if((length(x.names.sp) == 0) & (length(x.names.tp) == 0)){
        # non-spatial and temporal beta
          output$X <- X
          #dimnames(output$X)[[2]] <- x.names
        }
        else if((length(x.names.sp) > 0) & (length(x.names.tp) == 0)){
        # for spatial beta
          output$X <- X
          output$Xsp <- Xsp
          dimnames(output$X)[[2]] <- as.list(x.names)
          #dimnames(output$Xsp)[[2]] <- x.names.sp
        }  
        else if((length(x.names.sp) == 0) & (length(x.names.tp) > 0)){
        # for temporal beta
          output$X <- X
          output$Xtp <- Xtp
          dimnames(output$X)[[2]] <- as.list(x.names)
          #dimnames(output$Xtp)[[2]] <- x.names.tp
        }
        else if((length(x.names.sp) > 0) & (length(x.names.tp) > 0)){
        # for both spatial and temporal beta
          output$X <- X
          output$Xsp <- Xsp
          output$Xtp <- Xtp
          dimnames(output$X)[[2]] <- as.list(x.names)
          #dimnames(output$Xsp)[[2]] <- x.names.sp
          #dimnames(output$Xtp)[[2]] <- x.names.tp
        }
        else{
          stop("\n#\n## Error: \n#")
        }
     }
     if(Y.out==TRUE){
        output$Y <- Y
     }
    #
      output$accept <- accept
      output$call <- formula
    #
           output$phip <- as.matrix(out$phip[(nBurn+1):nItr])
           if(cov==4){
           output$nup <- as.matrix(out$nup[(nBurn+1):nItr])
           }
           output$sig2ep <- as.matrix(out$sig2ep[(nBurn+1):nItr])
           output$sig2etap <- as.matrix(out$sig2etap[(nBurn+1):nItr])
           if(length(x.names.sp) != 0){          
           output$sig2betap <- as.matrix(out$sig2betasp[(nBurn+1):nItr])
           }
           if(length(x.names.tp) != 0){          
           output$sig2deltap <- as.matrix(out$sig2deltap[(nBurn+1):nItr])
           output$sig2op <- as.matrix(out$sig2op[(nBurn+1):nItr])
           }
           output$betap <- matrix(out$betap[1:p,(nBurn+1):nItr],p,length((nBurn+1):nItr))
           if(length(x.names.sp) != 0){          
           output$betasp <- matrix(out$betasp[1:(n*q),(nBurn+1):nItr],n*q,length((nBurn+1):nItr))
           }
           if(length(x.names.tp) != 0){    
           output$betat0p <- matrix(out$betat0p[1:u,(nBurn+1):nItr],u,length((nBurn+1):nItr))
           output$betatp <- matrix(out$betatp[1:(u*T),(nBurn+1):nItr],u*T,length((nBurn+1):nItr))
           output$rhotp <- matrix(out$rhotp[1:u,(nBurn+1):nItr],u,length((nBurn+1):nItr))
             #if(initials$rhotp == 0){
             #  output$rhotp <- matrix(out$rhotp[1:u,(nBurn+1):nItr],u,length((nBurn+1):nItr))
             #}
           }
           output$op <- out$op[1:N,(nBurn+1):nItr]
           output$wp <- output$op-output$X%*%output$betap
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
           if(length(x.names.tp) != 0){          
           output$tp.covariate.names<-c(x.names.tp)
           }
           output$Distance.matrix <- coords.D
           output$coords <- coords
           output$n <- n
           output$r <- r
           output$T <- T
           output$p <- p
           if(length(x.names.sp) != 0){          
           output$q <- q
           }
           if(length(x.names.tp) != 0){          
           output$u <- u
           }
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
     #class(output) <- "spGP"
     output
    #
    #
}
##
## Prediction of Z_lt for the GP models using MCMC samples
##
 spGP.prediction<-function(nBurn, pred.data, pred.coords, posteriors,
                  tol.dist=2, Summary=TRUE)
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
      #
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
      if(method=="geodetic:km"){
           coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=TRUE))
      }
      else if(method=="geodetic:mile"){
           coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=FALSE))
      }
      else {
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
      pred.xsp<-Formula.matrix(call.f,data=pred.data)[[4]]
      pred.xtp<-Formula.matrix(call.f,data=pred.data)[[6]]
    #
       #
           phip<-posteriors$phip[(nBurn+1):nItr,]
           if(cov==4){
           nup<-posteriors$nup[(nBurn+1):nItr,]
           }
           else{
           nup<-0 # nup<-rep(0,itt)
           } 
           sig_ep<-posteriors$sig2ep[(nBurn+1):nItr,]
           sig_etap<-posteriors$sig2etap[(nBurn+1):nItr,]
           betap<-posteriors$betap[,(nBurn+1):nItr]
           op<-posteriors$op[,(nBurn+1):nItr] 	
      #
      if((is.null(pred.xsp)) & (is.null(pred.xtp))){
        out<-matrix(.C('z_pr_its_gp',as.integer(cov),as.integer(itt),as.integer(nsite), 
            as.integer(n),as.integer(r),as.integer(rT),as.integer(T),
            as.integer(p),as.integer(N),as.integer(predN),as.double(d), 
            as.double(d12),as.double(phip),as.double(nup),as.double(sig_ep),as.double(sig_etap),
            as.double(betap),as.double(posteriors$X),
            as.double(pred.x),as.double(op),as.integer(1),
            out=double(itt*predN))$out,predN,itt)
      }
      else if((!is.null(pred.xsp)) & (is.null(pred.xtp))){
      # spatial beta
        sig_betap<-posteriors$sig2betap[(nBurn+1):nItr,]
        betasp<-posteriors$betasp[,(nBurn+1):nItr]
        q<-posteriors$q 
        if (is.null(posteriors$Xsp)) {
           stop("Error: need to provide the fitted spatial covariate values.")
        }
        out<-.C('z_pr_its_gp_sp',as.integer(cov),as.integer(itt),as.integer(nsite), 
            as.integer(n),as.integer(r),as.integer(rT),as.integer(T),as.integer(p),
            as.integer(q),as.integer(N),as.integer(predN),as.double(d),as.double(d12),
            as.double(phip),as.double(nup),as.double(sig_ep),as.double(sig_etap),as.double(sig_betap),
            as.double(betap),as.double(betasp),as.double(posteriors$X),as.double(pred.x),
            as.double(posteriors$Xsp),as.double(pred.xsp),as.double(op),as.integer(1),
            betapred=double(q*nsite*itt),zpred=matrix(double(itt*predN),predN,itt))[27:28]
      } 
      else if((is.null(pred.xsp)) & (!is.null(pred.xtp))){
      # temporal beta
        sig_deltap<-posteriors$sig2deltap[(nBurn+1):nItr,]
        sig_op<-posteriors$sig2op[(nBurn+1):nItr,]
        betatp<-posteriors$betatp[,(nBurn+1):nItr]
        betat0p<-posteriors$betat0p[,(nBurn+1):nItr]
        u<-posteriors$u 
        if(is.null(posteriors$Xtp)) {
           stop("Error: need to provide the fitted temporal covariate values.")
        }
        if(is.null(posteriors$rhotp)){
           rhotp<-rep(1,u*itt)
        }
        else{
           rhotp<-posteriors$rhotp[,(nBurn+1):nItr]
        } 
        out<-matrix(.C('z_pr_its_gp_tp',as.integer(cov),as.integer(itt),as.integer(nsite), 
            as.integer(n),as.integer(r),as.integer(rT),as.integer(T),
            as.integer(p),as.integer(u),as.integer(N),as.integer(predN),as.double(d), 
            as.double(d12),as.double(phip),as.double(nup),as.double(sig_ep),as.double(sig_etap),
            as.double(sig_deltap),as.double(sig_op),as.double(betap),
            as.double(rhotp),as.double(betat0p),as.double(betatp),
            as.double(posteriors$X),as.double(pred.x),as.double(posteriors$Xtp),as.double(pred.xtp),
            as.double(op),as.integer(1),out=double(itt*predN))$out,predN,itt)
      }
      else if((!is.null(pred.xsp)) & (!is.null(pred.xtp))){
      # both spatial and temporal beta
        sig_betap<-posteriors$sig2betap[(nBurn+1):nItr,]
        betasp<-posteriors$betasp[,(nBurn+1):nItr]
        q<-posteriors$q 
        if (is.null(posteriors$Xsp)) {
           stop("Error: need to provide the fitted spatial covariate values.")
        }
        sig_deltap<-posteriors$sig2deltap[(nBurn+1):nItr,]
        sig_op<-posteriors$sig2op[(nBurn+1):nItr,]
        betatp<-posteriors$betatp[,(nBurn+1):nItr]
        betat0p<-posteriors$betat0p[,(nBurn+1):nItr]
        u<-posteriors$u 
        if(is.null(posteriors$Xtp)) {
           stop("Error: need to provide the fitted temporal covariate values.")
        }
        if(is.null(posteriors$rhotp)){
           rhotp<-rep(1,u*itt)
        }
        else{
           rhotp<-posteriors$rhotp[,(nBurn+1):nItr]
        } 
        out<-.C('z_pr_its_gp_sptp',as.integer(cov),as.integer(itt),as.integer(nsite), 
            as.integer(n),as.integer(r),as.integer(rT),as.integer(T),as.integer(p),as.integer(q),
            as.integer(u),as.integer(N),as.integer(predN),as.double(d),as.double(d12),
            as.double(phip),as.double(nup),as.double(sig_ep),as.double(sig_etap),as.double(sig_betap),
            as.double(sig_deltap),as.double(sig_op),as.double(betap),as.double(betasp),as.double(rhotp),
            as.double(betat0p),as.double(betatp),as.double(posteriors$X),as.double(pred.x),
            as.double(posteriors$Xsp),as.double(pred.xsp),as.double(posteriors$Xtp),as.double(pred.xtp),
            as.double(op),as.integer(1),
            betapred=double(q*nsite*itt),zpred=matrix(double(itt*predN),predN,itt))[35:36]
      }
      else{
         stop("Error: correctly define X variables")
      } 
      #
      output <- NULL
      #
      if((!is.null(pred.xsp)) & (is.null(pred.xtp))){
      # spatial beta
        output$pred.samples<-out$zpred
        output$pred.spbeta.samples<-array(c(out$betapred),dim=c(q,nsite,itt))
      }
      else if((!is.null(pred.xsp)) & (!is.null(pred.xtp))){
      # both spatial and temporal beta
        output$pred.samples<-out$zpred
        output$pred.spbeta.samples<-array(c(out$betapred),dim=c(q,nsite,itt))
      }
      else{
      # for fixed and dynamic beta
        output$pred.samples <- out
      } 
      ##
      out<-NULL
          #
          if(scale.transform=="NONE"){ 
              output$pred.samples <- output$pred.samples           
          }
          if(scale.transform=="SQRT"){ 
              output$pred.samples <- output$pred.samples^2
          }
          if(scale.transform=="LOG"){ 
              output$pred.samples <- exp(output$pred.samples)
          }
          # 
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
          #class(output) <- "spGP"
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
          #class(output) <- "spGP"
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
          #class(output) <- "spGP"
          output
      }
}
##
## K-step forecasts for the GP models 
##
 spGP.forecast <- function(nBurn, K, fore.data=NULL,
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
      #
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

      newr<-r

      if(is.data.frame(fore.data)){
       if((nsite*newr*K)!=dim(fore.data)[[1]]){
        stop("\n# Check the fore.data and/or fore.coords and/or spT.time function#\n")
       }
       fore.data$tmp<-rep(1,nsite*newr*K)
      }
      if(is.null(fore.data)){
        fore.data<-data.frame(tmp=rep(1,nsite*newr*K))
      }
      fore.x<-Formula.matrix(call.f,data=fore.data)[[2]]
      fore.xsp<-Formula.matrix(call.f,data=fore.data)[[4]]
      fore.xtp<-Formula.matrix(call.f,data=fore.data)[[6]]
    #
    #
      #
           method <- posteriors$distance.method
           coords <- posteriors$coords             
           spT.check.sites.inside(fore.coords, method, tol=tol.dist)
      #                                                           
        if(method == "geodetic:km"){
           #coords.D <- as.matrix(spT.geodist(Lon=fore.coords[,1],Lat=fore.coords[,2],KM=TRUE))
           coords.D <- posteriors$Distance.matrix
           coords.f.D<-as.matrix(spT.geodist(Lon=c(fore.coords[,1],coords[,1]),
                                 Lat=c(fore.coords[,2],coords[,2]),KM=TRUE))
           coords.f.D<-coords.f.D[1:nsite,(nsite+1):(nsite+n)] 
        }
        else if(method == "geodetic:mile"){
           #coords.D <- as.matrix(spT.geodist(Lon=fore.coords[,1],Lat=fore.coords[,2],KM=FALSE))
           coords.D <- posteriors$Distance.matrix
           coords.f.D<-as.matrix(spT.geodist(Lon=c(fore.coords[,1],coords[,1]),
                                 Lat=c(fore.coords[,2],coords[,2]),KM=FALSE))
           coords.f.D<-coords.f.D[1:nsite,(nsite+1):(nsite+n)] 
        }
        else {
           #coords.D <- as.matrix(dist(fore.coords, method, diag=TRUE, upper=TRUE))
           coords.D <- posteriors$Distance.matrix
           coords.f.D <- as.matrix(dist(rbind(fore.coords,coords), method, diag=TRUE, upper=TRUE))
           coords.f.D<-coords.f.D[1:nsite,(nsite+1):(nsite+n)]
        }
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
           betap<-matrix(posteriors$betap,p,nItr)
           betap<-betap[,(nBurn+1):nItr]
      #
      #
        w<-apply(posteriors$wp,1,median)
        w<-matrix(w,rT,n)
        w<-apply(w,2,median)
        #w<-apply(posteriors$wp[,(nBurn+1):nItr],2,function(x)apply(matrix(x,rT,n),2,median))
      #
      if((is.null(fore.xsp)) & (is.null(fore.xtp))){
      out<-matrix(.C('zlt_fore_gp_its',as.integer(cov),
           as.integer(itt),as.integer(K),as.integer(nsite),as.integer(n),
           as.integer(newr),as.integer(p),as.integer(newr*T),as.integer(T),
           as.integer(newr*K),as.integer(nsite*newr*K), 
           as.double(coords.D),as.double(t(coords.f.D)),as.double(phip),as.double(nup),
           as.double(sig_ep),as.double(sig_etap),as.double(fore.x),
           as.double(betap),as.double(w),as.integer(1),
           foreZ=double(nsite*newr*K*itt))$foreZ,nsite*newr*K,itt)
      }
      else if((!is.null(fore.xsp)) & (is.null(fore.xtp))){
      # forecast for spatially varying beta
      sig_betap<-posteriors$sig2betap[(nBurn+1):nItr,]
      betasp<-posteriors$betasp[,(nBurn+1):nItr]
      q<-posteriors$q 
      if (is.null(posteriors$Xsp)) {
           stop("Error: need to provide the fitted spatial covariate values.")
      }
      out<-matrix(.C('zlt_fore_gp_sp_its',as.integer(cov),
           as.integer(itt),as.integer(K),as.integer(nsite),as.integer(n),
           as.integer(newr),as.integer(p),as.integer(q),as.integer(newr*T),as.integer(T),
           as.integer(newr*K),as.integer(nsite*newr*K), 
           as.double(coords.D),as.double(t(coords.f.D)),as.double(phip),as.double(nup),
           as.double(sig_ep),as.double(sig_etap),as.double(sig_betap),as.double(fore.x),
           as.double(fore.xsp),as.double(betap),as.double(betasp),as.double(w),as.integer(1),
           foreZ=double(nsite*newr*K*itt))$foreZ,nsite*newr*K,itt)
      }
      else if((is.null(fore.xsp)) & (!is.null(fore.xtp))){
      # forecast for temporally varying beta
      u<-posteriors$u 
      sig_deltap<-posteriors$sig2deltap[(nBurn+1):nItr,]      
      sig_op<-posteriors$sig2op[(nBurn+1):nItr,]      
      betat0p<-posteriors$betat0p[,(nBurn+1):nItr]
      betatp<-posteriors$betatp[,(nBurn+1):nItr]
      if(is.null(posteriors$Xtp)) {
           stop("Error: need to provide the fitted temporal covariate values.")
      }
      if(is.null(posteriors$rhotp)){
           rhotp<-rep(1,u*itt)
      }
      else{
           rhotp<-posteriors$rhotp[,(nBurn+1):nItr]
      } 
      out<-matrix(.C('zlt_fore_gp_tp_its',as.integer(cov),
           as.integer(itt),as.integer(K),as.integer(nsite),as.integer(n),
           as.integer(newr),as.integer(p),as.integer(u),as.integer(newr*T),as.integer(T),
           as.integer(newr*K),as.integer(nsite*newr*K), 
           as.double(coords.D),as.double(t(coords.f.D)),as.double(phip),as.double(nup),
           as.double(sig_ep),as.double(sig_etap),as.double(sig_deltap),as.double(sig_op),
           as.double(fore.x),as.double(fore.xtp),
           as.double(betap),as.double(rhotp),as.double(betat0p),as.double(betatp),
           as.double(w),as.integer(1),
           foreZ=double(nsite*newr*K*itt))$foreZ,nsite*newr*K,itt)
      }
      else if((!is.null(fore.xsp)) & (!is.null(fore.xtp))){
      # forecast for spatially and temporally varying beta
      sig_betap<-posteriors$sig2betap[(nBurn+1):nItr,]
      betasp<-posteriors$betasp[,(nBurn+1):nItr]
      q<-posteriors$q 
      if (is.null(posteriors$Xsp)) {
           stop("Error: need to provide the fitted spatial covariate values.")
      }
      u<-posteriors$u 
      sig_deltap<-posteriors$sig2deltap[(nBurn+1):nItr,]      
      sig_op<-posteriors$sig2op[(nBurn+1):nItr,]      
      betat0p<-posteriors$betat0p[,(nBurn+1):nItr]
      betatp<-posteriors$betatp[,(nBurn+1):nItr]
      if(is.null(posteriors$Xtp)) {
           stop("Error: need to provide the fitted temporal covariate values.")
      }
      if(is.null(posteriors$rhotp)){
           rhotp<-rep(1,u*itt)
      }
      else{
           rhotp<-posteriors$rhotp[,(nBurn+1):nItr]
      } 
      out<-matrix(.C('zlt_fore_gp_sptp_its',as.integer(cov),as.integer(itt),as.integer(K),
           as.integer(nsite),as.integer(n),as.integer(newr),as.integer(p),as.integer(q),
           as.integer(u),as.integer(newr*T),as.integer(T),as.integer(newr*K),as.integer(nsite*newr*K), 
           as.double(coords.D),as.double(t(coords.f.D)),as.double(phip),as.double(nup),
           as.double(sig_ep),as.double(sig_etap),as.double(sig_betap),as.double(sig_deltap),as.double(sig_op),
           as.double(fore.x),as.double(fore.xsp),as.double(fore.xtp),
           as.double(betap),as.double(betasp),as.double(rhotp),as.double(betat0p),as.double(betatp),
           as.double(w),as.integer(1),foreZ=double(nsite*newr*K*itt))$foreZ,nsite*newr*K,itt)
      }
      else{
        stop("Error")
      }
      #
      output <- NULL
          #
          if(scale.transform=="NONE"){ 
      output$fore.samples <- out
          }
          if(scale.transform=="SQRT"){ 
      output$fore.samples <- out^2
          }
          if(scale.transform=="LOG"){ 
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
          cat("# Summary statistics are not given,\n#   because the number of forecast samples (",itt,") are too small.\n# nBurn = ",nBurn,". Iterations = ",nItr,".", "\n")
          #cat("# Columns = iterations, Rows = kth step forecasts in each site.\n")
          cat("##", "\n")
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   output$computation.time<-comp.time
    #
          #class(output) <- "spGP"
          output
         }
         else { 
          cat("##", "\n")
          cat("# Forecast samples and summary statistics are given.\n# nBurn = ",nBurn,". Iterations = ",nItr,".", "\n")
          cat("##", "\n")
          #
          szp<-spT.Summary.Stat(output$fore.samples[,])
          # 
          output$Mean <- matrix(szp$Mean,newr*K, nsite)
          output$Median <- matrix(szp$Median,newr*K, nsite)
          output$SD <- matrix(szp$SD,newr*K, nsite)
          output$Low <- matrix(szp[,4],newr*K, nsite)
          output$Up <- matrix(szp[,5],newr*K, nsite)
    #
   end.time <- proc.time()[3]
   comp.time<-end.time-start.time
   comp.time<-fnc.time(comp.time)
   output$computation.time<-comp.time
    #
          #class(output) <- "spGP"
          output
         }
      }
}
##
## MCMC fit and predictions for the GP models
##
spGP.MCMC.Pred<-function(formula, data=parent.frame(), time.data, 
            coords, pred.coords, priors, initials, pred.data, nItr, nBurn, 
			report=1, tol.dist=2, distance.method="geodetic:km", cov.fnc="exponential",
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
         priors<-priors.checking.gp(priors)
    #
         shape_e<-N/2+priors$prior_a
         shape_eta<-N/2+priors$prior_a
    #
         zm <- matrix(Y,rT,n)
         zm <- apply(zm,1,median,na.rm=TRUE)
         zm <- rep(zm,n)
         zm <- cbind(Y,zm)
         zm[is.na(zm[,1]),1]=zm[is.na(zm[,1]),2]
         zm <- zm[,1]
         names(zm)<-NULL
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
      initials<-initials.checking.gp(initials,zm,X,Xsp,Xtp,n,r,T,coords.D)
    #
         o <- zm
    #
    if(spatial.decay$type=="FIXED"){
         spdecay <- 1
		 if(is.null(spatial.decay$value)){
		 spatial.decay$value <- (3/max(c(coords.D))) 
		 }
         init.phi <- spatial.decay$value 
         tuning <- 0; phis<-0; phik<-0; 
		 phi_a<-0; phi_b<-0
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
      if (length(initials$beta) != p){
         stop("Error: need to specify correct number of parameters for beta.")
      }
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
        out <- .C("GIBBS_sumpred_txt_gp", as.integer(aggtype),as.double(flag), as.integer(nItr),
            as.integer(nBurn), as.integer(n), as.integer(T), as.integer(r), 
            as.integer(rT), as.integer(p), as.integer(N), as.integer(report),
            as.integer(cov), as.integer(spdecay), as.double(shape_e), as.double(shape_eta),
            as.double(phi_a), as.double(phi_b),			
            as.double(priors$prior_a), as.double(priors$prior_b), 
            as.double(priors$prior_mubeta), as.double(priors$prior_sigbeta),
            as.double(priors$prior_omu), as.double(priors$prior_osig),
            as.double(init.phi), as.double(tuning),
            as.double(phis), as.integer(phik), as.double(coords.D), 
            as.double(initials$sig2eps), as.double(initials$sig2eta),
            as.double(initials$beta), as.double(X), as.double(zm), as.double(o),
            as.integer(1), 
            as.integer(nsite), as.integer(predN), as.double(d12), as.double(pred.x), 
            as.integer(trans), accept=double(1), gof=double(1),penalty=double(1))[41:43] 
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
           tmp<-read.table('OutGP_Values_Parameter.txt',sep='',header=FALSE)
           tmp<-tmp[(nBurn+1):nItr,]
           tmp<-spT.Summary.Stat(t(tmp))
           if(cov==4){
           row.names(tmp)<- c(x.names,"sig2eps","sig2eta","phi","nu")
           }
           else{
           row.names(tmp)<- c(x.names,"sig2eps","sig2eta","phi")
           } 
           out$parameters<-round(tmp, 4)
           tmp<-NULL
     #
           tmp<-read.table('OutGP_Stats_FittedValue.txt', sep=',',header=FALSE)
           names(tmp)<- c("Mean","SD")
           out$fitted<-round(tmp, 4)
           tmp<-NULL
     #
           tmp<-read.table('OutGP_Stats_PredValue.txt', sep=',',header=FALSE)
           names(tmp)<- c("Mean","SD")
           out$prediction<-round(tmp, 4)
           tmp<-NULL
     #
       if(annual.aggregation=="ave"){
           tmp<-read.table('OutGP_Annual_Average_Prediction.txt', sep='',header=FALSE)
           tmp<-spT.Summary.Stat(t(tmp))
           out$an.agr.pred.average<-round(tmp, 4)
           tmp<-NULL
       }   
       if(annual.aggregation=="an4th"){
           tmp<-read.table('OutGP_Annual_4th_Highest_Prediction.txt', sep='',header=FALSE)
           tmp<-spT.Summary.Stat(t(tmp))
           out$an.agr.pred.an4th<-round(tmp, 4)
           tmp<-NULL
       }   
       if(annual.aggregation=="w126"){
           tmp<-read.table('OutGP_Annual_w126_Prediction.txt', sep='',header=FALSE)
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
       #class(out) <- "spGP"
      #
       out
      #
      #
}
##
##
##