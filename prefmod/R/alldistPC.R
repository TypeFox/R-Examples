"alldistPC" <-
 function(formula,
                    random=~1,
                    family = gaussian(),
                    data ,
                    k = 1,
                    random.distribution="np",
                    tol = 0.5,
                    startp = NULL,
                    offset,
                    weights,
                    pluginz,
                    na.action,
                    EMmaxit=500,
                    EMdev.change=0.001,
                    lambda=0,
                    damp=FALSE,
                    damp.power=1,
                    spike.protect=0,
                    sdev,
                    shape,
                    plot.opt=0,
                    verbose=FALSE,
                    pr.it = FALSE,
                    ...)
{
  # R function alldist in package npmlreg.  NPML/GQ for overdispersed GLMs.
  # Type  ?npmlreg for licence and copyright information.
  # Version 0.30, 06-03-2006


  call <- match.call()
    if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    stop("`family' not recognized")
  }

  data <- as.data.frame(data)   #0.34-1
  ddim <- dim(data)
  mf   <- match.call(expand.dots = FALSE)


  # Test for incorrect offset specification in formula object
  testoffset<-try(is.null(attr(terms(formula(mf)),"offset")),silent=TRUE)
  if (!(class(testoffset)=="try-error" || testoffset)){
      stop("Please specify offset as separate argument outside the model formula.")
  }

  m    <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0)
  mf   <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- try(eval(mf,  parent.frame()), silent=TRUE)
  ######if (class(mf)=="try-error"){
  ######    if (!missing(offset) && length(offset) != ddim[1]) {
  ######    stop("Number of offsets is ", length(offset), ", should equal ", ddim[1], " (number of observations)")
  ######    }
  ######    if (!missing(weights) && length(weights) != ddim[1]){
  ######    stop("Number of weights is ", length(weights), ", should equal ", ddim[1], " (number of observations)")
  ######    }
  ######    stop(geterrmessage())
  ######}

  ######offset  <- model.offset(mf)
  ######weights <- model.weights(mf)
  X  <- Y <- XZ <- NULL
#####  Y  <- model.response(mf, "numeric") # response
  Y  <- data$y  # response
  Ym <- is.matrix(Y)

  N  <- NROW(Y)  #corresponds to ddim[1] if there are no missing values
  if (is.null(offset)){offset<-rep(0,N) }
  if (is.null(weights)){weights<-rep(1,N)}
#####  data$offset <- numeric(1); data$pweights<-numeric(1)
#####  data  <- if (is.matrix(Y)) data[dimnames(Y)[[1]],] else  data[names(Y),] #omit missing values
  data$offset<-offset;  data$pweights<-weights

  # initial fit and simple glm for k=1
  fit    <- glm(formula, family=family, weights=pweights, offset=offset, data=data,...)

  names0 <- dimnames(data)[[1]]
  w0     <- fit$prior.weights  #store prior weights for output
  off0   <- fit$offset;   names(off0)<-names0   #store offset for output
  Y      <- fit$y
  l0     <- length(fit$coef)

  ######if(family$family=="binomial"){
  ######         data$pweights<- data$pweights^Ym         # these are the actual case weights
  ######         YP<- binomial.expand(Y,1,w0/data$pweights)   # w0/data$pweights are the actual pop. sizes
  ######         Y<- YP[[1]]; PY<-YP[[2]];r<-YP[[3]]; n<-YP[[4]]
  ######}

  #w    <- data$pweights #  removed with version 0.32


  sdev.miss <- missing(sdev)
  shape.miss <-missing(shape)


  if (family$family=="gaussian"){
        #sdev <- ifelse(sdev.miss, sqrt(sum(w*(Y-fitted(fit))^2)/(sum(w)-l0)), sdev)  # identical to:
        sdev  <-  ifelse(sdev.miss, sqrt(summary(fit)$dispersion), sdev)
        shape <- 0
  } else if (family$family =="Gamma") {
         # Estimate sdev from residuals on linear predictor scale, see Einbeck & Hinde (2006):
        sdev  <- ifelse(sdev.miss, sqrt(switch(family$link,
                        "log"= sum(data$pweights*(log(Y)-log(fitted(fit)))^2)/(sum(data$pweights)),
                        "inverse"= sum(data$pweights*(1/Y-1/(fitted(fit)))^2)/(sum(data$pweights)),
                        "identity"= sum(data$pweights*(Y-fitted(fit))^2)/(sum(data$pweights)),
                        )), sdev)
        #shape <- ifelse(shape.miss,(sum(w)-l0)/sum(w*((Y-fitted(fit))/fitted(fit))^2), shape)#ident. to:
        shape <- ifelse(shape.miss,1/summary(fit)$dispersion, shape)
  } else {
        sdev  <- 1
        shape <- 0
  }


  ML.dev0 <- -2*sum(data$pweights*switch(family$family,
                "gaussian"= dnorm(fit$y, fitted(fit), sdev, log=TRUE),
                "poisson" = dpois(fit$y, fitted(fit), log=TRUE),
                "binomial"= dbinom(Y[,1],Y[,1]+Y[,2], fitted(fit),log=TRUE),
                "Gamma"   = dgamma(fit$y,  shape=shape, scale=fitted(fit)/shape, log=TRUE),
                ))


  if (k == 1) {

      fit <- c( fit[c(1,2,3,8,9)],
              disparity = ML.dev0,
              deviance = fit$dev,
              fit[c(12,16,17, 18)],
              call = call,
              formula = formula,
              random = "none",
              data = list(subset(data,select=1:ddim[2])),
              model = list(model.matrix(fit)),
              weights = list(w0),
              offset = list(off0),
              mass.points = list(fit$coef[1]),
              masses = list(c("MASS1"=1)),
              sdev = list(list(sdev=sdev, sdevk=sdev)),
              shape = list(list(shape=shape,shapek=shape)),
              post.prob = list(matrix(1,N,1,dimnames=list(names0,"") )),
              post.int =  list(fit$coef[1]),
              ebp = family$linkfun(fit$fitted),
              EMiter = 0,
              EMconverged = "none",
              lastglm = list(fit),
              Misc = list(list(lambda=lambda))
              )
      class(fit)<-'pattNPML'
      return(fit)
  } else if (!(k %in% 1:21)){
      stop("This choice of k is not supported.")
     }


  # Expand the response
  if(family$family=="binomial"){
#####      YP<- binomial.expand(Y,k,rep(1,N)); Y<- YP[[1]]; PY<-YP[[2]];r<-YP[[3]]; n<-YP[[4]]
  }  else Y <- rep(Y,k)


  datak  <- expand(data,k)# expand data.
  kindex <- rep(1:k,rep(N,k))# index for the mixtures, for version 0.31 or higher only used for allvc
  #tmp   <- hermite(k) # grab weights and abscissas
  tmp    <- gqz(k,minweight=1e-14)  #from version 0.32-2
  #z0    <- tmp$z    # for GQ - masspoints
  z0     <- -tmp$l
  #z     <- rep(tmp$z,rep(N,k))
  z      <- rep(-tmp$l,rep(N,k))
  p      <- tmp$w
  rform  <- random

  #  Generate the random design matrix and append to fixed matrix
  mform <- strsplit(as.character(random)[2],'\\|')[[1]]
  mform <- gsub(' ', '',mform)

  if (length(mform)==2){stop("Please use function allvc for two-level models")}

  # Expand offset and weights and attach them to expanded data frame
  offset   <- datak$offset
  pweights <- datak$pweights

  if (random.distribution=='np'){
      # Nonparametric random effect
#      X <- model.matrix(formula,datak)[,-1,drop=FALSE]
      X <- model.matrix(formula,datak)[,,drop=FALSE]      ### intercept included
      datak$MASS <- gl(k,N)
#      if (mform=='1') random <- formula(~MASS-1) else {
#          # Nonparametric random coefficient
#          random <- formula(paste('~ MASS + ',paste(mform,'MASS',sep=":",collapse='+'), '-1',sep=''))
      if (mform=='1') random <- formula(1) else {
          # Nonparametric random coefficient
#          random <- formula(paste('~ ',paste(mform,'MASS',sep=":",collapse='+'), '-1',sep=''))
          random <- formula(paste('~ ',paste(paste("(",mform,")",sep=''),'MASS',sep=":",collapse='+'), '-1',sep=''))
#          random <-formula("~-1+a:MASS+b:MASS+c:MASS+d:MASS+e:MASS")
      }
  } else {
      X <- model.matrix(formula,datak)
      if (mform=='1') {
          random <- formula(~ z - 1 )
      } else {
          stop("random coefficients only allowed with option 'np'")
      }
  }


  Z <- model.matrix(random,datak)
  if (dim(X)[1]!= dim(Z)[1]){ cat("The missing value routine cannot cope with this model. Please specify the random term also as fixed term and try again. " )}
  XZ <- cbind(X,Z)

#####  if (missing(pluginz)){sz<-tol* sdev*z} else {sz<-rep(pluginz-fit$coef[[1]],rep(N,k))}
  if (is.null(pluginz)){sz<-tol* sdev*z} else {sz<-rep(pluginz-fit$coef[[1]],rep(N,k))}
  Eta <- fit$linear.predictor + sz

  if (random.distribution=="np"){
      tol<- max(min(tol,1),1-damp)  #For tol >  1 or damp=F no Damping
      if(length(fit$coef)==1){
          followmass<-matrix(Eta[(1:k)*N],1,k)-offset[(1:k)*N]
      } else {
          followmass<-matrix(fit$coef[1]+sz[(1:k)*N],1,k)
  }
  } else {
      followmass<-NULL; tol<-1
  }

  Mu <- family$linkinv(Eta) # expanded fitted values

  # Calculate loglikelihood for fixed model
  f <- switch(family$family,
######              "gaussian" = dnorm(Y,Mu,tol*sdev,log=TRUE),
              "poisson"  = dpois(Y,Mu,log=TRUE),
######              "binomial" = dbinom(r,n,Mu,log=TRUE),
######              "Gamma"    = dgamma(Y,shape=shape/tol^2,scale=Mu*tol^2/shape ,log=TRUE),
              )


  # Calculate the weights from initial model

################# inserted #################################
if (!is.null(startp)) p<-startp
#cat("RANDOM STARTING PROPORTIONS:", p, "\n")
############################################################
  tmp <-weightslogl.calc.w(p,matrix(f,ncol=k),pweights[1:N])
  w <- tmp$w
#### new randow starts ############## 2011-11-11
  w[,]<-runif(length(w))
  w<-t(apply(w,1,function(x)x/sum(x)))
#### new randow starts ############## 2011-11-11

 # Initialize for EM loop
  ML.dev <- ML.dev0
  iter <- ml<- 1
  converged <- FALSE
  sdevk<-rep(sdev,k);  shapek<-rep(shape,k)

  p.iter<-0

  while (iter <= EMmaxit && (!converged || (iter<=9 && random.distribution=='np' && damp && (family$family=="gaussian" && sdev.miss || family$family=="Gamma"&& shape.miss)  ))){
  ##########Start of EM ##########

      # print dot at each iteration (max 50/line)
      if(pr.it) {
          if (p.iter>50){
            cat("\n")
            p.iter<-0
          }
          p.iter<-p.iter+1
          cat(".")
          flush.console()
      }

   #####   if (verbose){ cat(iter,"..") }

      fit <- try(glm.fit(x=XZ, y=Y, weights = as.vector(w)*pweights, family = family, offset=offset, ...), silent=TRUE)
      if (class(fit)=="try-error"){
              stop("Singularity or Likelihood-Spike at iteration #",iter,".
              Check model specification,  enable spike protection or smooth among components.")
      }

      if (random.distribution=="np"){   #Save EM Trajectories
          masspoint<- fit$coef[l0:(l0+k-1)]
          followmass<-rbind(followmass, masspoint)
      }

      Mu <- fitted(fit)

      #Unequal component dispersion parameters

      if (family$family=="gaussian"){
 ######## if (sdev.miss){ sdev<- sqrt(sum((as.vector(w)*pweights)*(Y-Mu)^2)/sum(as.vector(w)*pweights))}
 ######## sdevk<-rep(sdev,k)
 ######## if (lambda!=0){
 ########     for (l in 1:k){
 ########       wk<-matrix(1,k,N); wk[1:k,]<-dkern(1:k,l,k,lambda);wk<-t(wk)
 ########       sdevk[l] <-  sqrt(sum(wk* as.vector(w)*pweights *(Y-Mu)^2)/sum(wk*as.vector(w)*pweights))
 ########     }
 ########     sk<-rep(sdevk,rep(N,k))
 ######## } else {
 ########     sk<-sdev
 ######## }
      }  else {
          sdevk <-rep(NA,k)
      }

      if (family$family=="Gamma"){
   ########if (shape.miss) { shape<-(sum(as.vector(w)*pweights))*1/sum(as.vector(w)*pweights*((Y-fitted(fit))/fitted(fit))^2)}
   ########shapek<-rep(shape,k)
   ########if (lambda!=0){
   ########     for (l in 1:k){
   ########       wk<-matrix(1,k,N); wk[1:k,]<-dkern(1:k,l,k,lambda);wk<-t(wk)
   ########       shapek[l] <- sum(wk*as.vector(w)*pweights)/ sum(wk* as.vector(w)*pweights*((Y-Mu)/Mu)^2)
   ########     }
   ########     shk<-rep(shapek,rep(N,k))
   ########} else {
   ########     shk<-shape
   ########}
      } else {shapek<-rep(NA,k)}

      # Calculate loglikelihood for expanded model for this iteration
      f <- matrix(switch(family$family,
   ######           "gaussian"=dnorm(Y,Mu,(1-(1-tol)^(damp.power*iter+1))*sk,log=TRUE),
              "poisson" =dpois(Y,Mu,log=TRUE),
   ######           "binomial"=dbinom(r,n,Mu,log=TRUE),
   ######           "Gamma"=dgamma(Y,shape=shk/(1-(1-tol)^(damp.power*iter+1))^2,scale=Mu*(1-(1-tol)^(damp.power*iter+1))^2/shk,log=TRUE),
              ),nrow=N, ncol=k)


      # Calculate the component proportions from the weights
      if (random.distribution=='np') p <- as.vector(apply(w*pweights,2,sum))/sum(pweights[1:N]) #16-03-06

      # Calculate updated weights and loglikehood
      tmp <- weightslogl.calc.w(p,f,pweights[1:N])
      w   <- tmp$w
      ML.dev[iter+1] <- ifelse(is.na(tmp$ML.dev), Inf, tmp$ML.dev)
      if (ML.dev[iter+1]>ML.dev0) {ml<-ml+1}  #only relevant for graphical output
      converged <- abs(ML.dev[iter+1] - ML.dev[iter])< EMdev.change
      iter <- iter + 1

      #Likelihood Spike Protection
      if (random.distribution != 'gq' && spike.protect!=0){
          if (family$family=='gaussian' && abs(min(sdevk/(fit$coef[(length(fit$coef)-k+1):length(fit$coef)]))) <0.000001*spike.protect){break}
          if (family$family=='Gamma' && abs(max(shapek/(fit$coef[(length(fit$coef)-k+1):length(fit$coef)])))> 10^6*spike.protect){break}
    }
  }###########################End of EM loop#############

#####  if (verbose){
      cat("\n")
      if (converged){
#####          cat("EM algorithm met convergence criteria at iteration # ", iter-1,"\n")
      } else {
        cat("EM algorithm failed to meet convergence criteria at iteration # ",
        iter-1,"\n")
      }
#####

  Deviance<- switch(family$family,
########               "gaussian"= sdev^2*ML.dev[iter]-sdev^2* sum(data$pweights[1:N] * log(2*pi*sdev^2)),
              "poisson" = ML.dev[iter] +2*sum(data$pweights[1:N]*(-Y[1:N]+Y[1:N]*log(Y[1:N]+ (Y[1:N]==0))-lfactorial(Y[1:N]))),
########               "binomial"= ML.dev[iter] +2*sum(data$pweights[1:N]*(lfactorial(n)-lfactorial(r)-lfactorial(n-r) - n*log(n) + r*log(r+(r==0))+(n-r)*log(n-r+((n-r)==0)))[1:N]),
              "Gamma"   = 1/shape*ML.dev[iter]+2/shape*(sum(data$pweights[1:N])*shape*(log(shape)-1)-sum(data$pweights[1:N])*lgamma(shape)-sum(data$pweights[1:N]*log(Y[1:N]))),
              )

  mass.points   <- masses <- NULL
  np            <- length(fit$coef)
  ebp           <- apply(w*matrix(fit$linear.predictor,N,k,byrow=FALSE),1,sum)  #Emp. Bayes Pred. (Aitkin, 96)
  ebp.fitted    <- family$linkinv(ebp)
  ebp.residuals <- Y[1:N]- ebp.fitted
  names(ebp)    <- names(ebp.fitted) <- names(ebp.residuals) <- names0
  if (is.na(fit$coefficients[np])){length(fit$coefficients)<-np<-np-1}# if one variable is random and fixed
  m <- seq(1,np)[substr(attr(fit$coefficients,'names'),1,4)=='MASS']
  mass.points   <- fit$coefficients[m]
  post.prob     <- matrix(w, nrow=N, byrow=FALSE, dimnames=list(names0, 1:k) )

   if ((plot.opt==1 || plot.opt==2) && par("mfrow")[1]>2) { #Write tol values as plot title if alldist is called from tolfind:
           plot.main <- substitute("tol"== tol, list(tol=tol))
      } else {
           plot.main <- c("")
  }

  if (plot.opt==3 && random.distribution=="np"){
      par(mfrow=c(2,1),cex=0.5,cex.axis=1.5,cex.lab=1.5)
  }

  if (plot.opt==1|| plot.opt==3){
      if  ((family$family=="gaussian" && sdev.miss|| family$family=="Gamma" && shape.miss) && damp  && random.distribution=='np' && iter>=max(8,ml+1)){
          #Linear interpolation for initial cycles
          ML.dev[2: max(7,ml)]<-ML.dev0+ 1:max(6,ml-1)/ max(7,ml)*(ML.dev[max(8,ml+1)]-ML.dev0)
      }
      plot(0:(iter-1),ML.dev, col=1,type="l",xlab='EM iterations',ylab='-2logL', main= plot.main )
      if (verbose){ cat("Disparity trend plotted.\n")}
  }


  if (random.distribution=="np") {
      #masses <- apply(w,2,mean)
      masses <- as.vector(apply(w*pweights,2,sum))/sum(pweights[1:N]) #16-03-05
      names(masses) <- paste('MASS',1:k,sep='')

      if (family$family=="binomial"){
 ######       R0 <- family$linkfun(PY[1:N])
      } else  {
          R0 <- family$linkfun(Y[1:N])
      }

      #Residuals R= y-Xb
      #if (family$family=="binomial"){
      #    R0<- switch(family$link, "log"= log(PY[1:N]),
      #                                    log(PY[1:N]/(1-PY[1:N])))
      #}  else if (family$family=="poisson") {
      #    R0<- switch(family$link, "log"= log(Y[1:N]),
      #                             "sqrt"=sqrt(Y[1:N]),
      #                             "identity"=Y[1:N])
      #} else if (family$family=="Gamma") {
      #    R0<- switch(family$link, "log"= log(Y[1:N]),
      #                             "inverse"=1/Y[1:N],
      #                             "identity"=Y[1:N] )
      #} else {
      #    R0<- switch(family$link,  "log"= log(Y[1:N]),
      #                              "inverse"=1/Y[1:N],
      #                              "identity"=Y[1:N])
      #}



      if(dim(X)[2]>0){
        R <- R0 - X[1:N,]%*%matrix(fit$coef[1:dim(X)[2]])-offset[1:N]
      } else {
        R <- R0 - offset[1:N]
      }
      R<-as.vector(R);  names(R)<-names0


      if (mform=='1' && any(is.finite(R))){
           ylim <- c(min(R[is.finite(R)]), max(R[is.finite(R)]))  #29/06/2006
           # ylim<- c(min(na.omit(R)),max(na.omit(R)))
           #if (abs(ylim[1]) == Inf){ylim[1] <- min(followmass[,])}
           # if (abs(ylim[2]) == Inf){ylim[2] <- max(followmass[,])}
      } else  {
            ylim <- c(min(followmass[,]),max(followmass[,]))
      }

      if(any(is.na(ylim)) &  plot.opt >1 ){
              cat("Singularity: EM Trajectory plot not available.", "\n");
              plot.opt<-min(plot.opt,1)
      }

      if (plot.opt==2|| plot.opt==3){
            plot(0:(iter-1),followmass[,1],col=1,type='l',ylim=ylim,ylab='mass points',xlab='EM iterations',  main=plot.main )
            for (i in 1:k){ lines(0:(iter-1), followmass[,i],col=i)
                        if (mform=='1'){ points(rep(iter-1,length(R)),R)}}
            if (verbose){ cat("EM Trajectories plotted.\n")}
      }

      fit <- c( fit[1],
                residuals = list(ebp.residuals),
                fitted.values = list(ebp.fitted),
                fit[c(8,9)],
                disparity = ML.dev[iter],
                deviance = Deviance,
                fit[12],
                df.residual = N-np-k+1,
                df.null= N-1,
                fit[18],
                call = call,
                formula = formula,
                random = rform,
              ##  data = list(subset(data,select=1:ddim[2])),
                data = list(data),
                model = list(XZ),
                weights = list(w0),
                offset = list(off0),
                mass.points = list(mass.points),
                masses = list(masses),
                sdev = list(list(sdev=sdev, sdevk=sdevk)),
                shape = list(list(shape=shape,shapek=shapek)),
                post.prob = list(post.prob),
                post.int =  list(post.prob %*% mass.points[1:k]),
                ebp = list(ebp),
                EMiter = iter - 1,
                EMconverged = converged,
                lastglm = list(fit),
                Misc = list(list(Disparity.trend=ML.dev,EMTrajectories=followmass,res=R,ylim=ylim, lambda=lambda,mform=mform))
            )
      class(fit) <- 'pattNPML'

  ##### else {
  #####   mass.points<-fit$coef[1]+fit$coef[np]*z0
  #####   fit <- c( fit[1],
  #####             residuals = list(ebp.residuals),
  #####             fitted.values = list(ebp.fitted),
  #####             fit[c(8,9)],
  #####             disparity = ML.dev[iter],
  #####             deviance = Deviance,
  #####             fit[12],
  #####             df.residual = N-np,
  #####             df.null= N-1,
  #####             fit[18],
  #####             call = call,
  #####             formula = formula,
  #####             random = rform,
  #####             data = list(subset(data,select=1:ddim[2])),
  #####             model = list(XZ),
  #####             weights = list(w0),
  #####             offset = list(off0),
  #####             mass.points=list(mass.points),
  #####             masses=list(gqz(k, minweight=1e-14)$w),
  #####             sdev=list(list(sdev=sdev, sdevk=sdevk)),
  #####             shape=list(list(shape=shape,shapek=shapek)),
  #####             post.prob=list(post.prob),
  #####             post.int =  list(post.prob %*% mass.points[1:k]),
  #####             ebp=list(ebp),
  #####             EMiter=iter - 1,
  #####             EMconverged=converged,
  #####             lastglm =list(fit),
  #####             Misc=list(list(Disparity.trend=ML.dev, lambda=lambda,mform=mform))
  #####        )
  #####   class(fit) <- 'glmmGQ'
  }
  fit
}
