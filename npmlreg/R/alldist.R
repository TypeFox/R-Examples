"alldist" <- 
 function(formula,
                    random=~1,
                    family = gaussian(),
                    data ,
                    k = 4,
                    random.distribution="np",
                    tol = 0.5,
                    offset,
                    weights,
                    pluginz,
                    na.action,
                    EMmaxit = 500,
                    EMdev.change = 0.001,
                    lambda = 0,
                    damp = TRUE,
                    damp.power = 1,
                    spike.protect = 0,
                    sdev,
                    shape,
                    plot.opt = 3,
                    verbose = TRUE,
                    ...)
{

  # R function alldist in package npmlreg.  NPML/GQ for overdispersed GLMs.
  # Type  ?npmlreg for licence, copyright, and version information.
  
  call <- match.call()
    if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("`family' not recognized")
  }

  data <- as.data.frame(data)   #0.34-1
  ddim <- dim(data)
  mf   <- match.call(expand.dots = FALSE)
 
  # intercept term removed?
  int.removed<-FALSE
  if (random.distribution=='np' && max(length(grep('- 1', deparse(formula(mf)))),length(grep('-1', deparse(formula(mf))))) >0 ){
     # stop(" term '-1'  in model formula not supported for k>1 & random.distribution='np'. ")
     int.removed<-TRUE
  }

 
  # Test for incorrect offset specification in formula object
  testoffset<-try(is.null(attr(terms(formula(mf)),"offset")),silent=TRUE)
  if (!(class(testoffset)=="try-error" || testoffset)){
      stop("Please specify offset as separate argument outside the model formula.")
  }
  
  # Extract variables from call and set up initial (fixed effect) model
  m    <- match(c("formula", "data", "subset", "weights", "na.action", 
               "etastart", "mustart", "offset"), names(mf), 0)
  mf   <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- try(eval(mf,  parent.frame()), silent=TRUE)
  if (class(mf)=="try-error"){
      if (!missing(offset) && length(offset) != ddim[1]) {
      stop("Number of offsets is ", length(offset), ", should equal ", ddim[1], " (number of observations)")
      }
      if (!missing(weights) && length(weights) != ddim[1]){
      stop("Number of weights is ", length(weights), ", should equal ", ddim[1], " (number of observations)")
      }
      stop(geterrmessage())
  }
  offset  <- model.offset(mf) 
  weights <- model.weights(mf)
  X  <- Y <- XZ <- NULL
  Y  <- model.response(mf, "any") # response  # for 0.42 changed from "numeric to "any"
  Ym <- is.matrix(Y)
  N  <- NROW(Y)  # corresponds to ddim[1] if there are no missing values
  
  # Set up weights and offset for initial glm
  if (is.null(offset)){offset  <-rep(0,N) }
  if (is.null(weights)){weights<-rep(1,N)}   
  data$offset <- numeric(1); data$pweights<-numeric(1)
  data  <- if (is.matrix(Y)) data[dimnames(Y)[[1]],] else  data[names(Y),] # omit missing values   
  data$offset<-offset;  data$pweights<-weights
   
  # Extract variable names from random part 
  rform   <- random
  mform   <- strsplit(as.character(random)[2],'\\|')[[1]]
  mform   <- gsub(' ', '',mform)
  
  if (length(mform)==2){stop("Please use function allvc for two-level models")}
  if (random.distribution=='gq' && mform!="1"){stop("Random coefficient models are only supported for random.distribution='np'.")}
  
 if(int.removed && mform=="1"){
     k<-1
     cat("k was set equal to 1, since this model cannot have a random component. \n")
  }   
  
  # Initial fit and simple glm for k=1
  fit    <- glm(formula, family=family, weights=pweights, offset=offset, data=data,...)
  names0 <- dimnames(data)[[1]] 
  w0     <- fit$prior.weights  #store prior weights for output
  off0   <- fit$offset;   names(off0)<-names0   #store offset for output
  Y      <- fit$y  
  l0     <- length(fit$coef)
  tol0   <- tol # for main title of graphical output 
  
  # For binomial models, check weights: case weights, or number of trials?
  if(family$family=="binomial"){
           data$pweights<- data$pweights^Ym         # these are the actual case weights
           YP<- binomial.expand(Y,1,w0/data$pweights)   # w0/data$pweights are the actual pop. sizes 
           Y<- YP[[1]]; PY<-YP[[2]];r<-YP[[3]]; n<-YP[[4]]                 
  } 
  
  # Initial estimates of  standard deviation and/or shape parameter (response distr.)       
  sdev.miss <- missing(sdev)
  shape.miss <-missing(shape)
  if (family$family=="gaussian"){
        sdev  <-  ifelse(sdev.miss, sqrt(summary(fit)$dispersion), sdev)
        shape <- 0
  } else if (family$family =="Gamma" ||family$family =="inverse.gaussian" ) {
         # Estimate sdev from residuals on linear predictor scale, see Einbeck & Hinde (2006):
        sdev  <- ifelse(sdev.miss, sqrt(switch(family$link,
                        "log"= sum(data$pweights*(log(Y)-log(fitted(fit)))^2)/(sum(data$pweights)),
                        "inverse"= sum(data$pweights*(1/Y-1/(fitted(fit)))^2)/(sum(data$pweights)),
                        "identity"= sum(data$pweights*(Y-fitted(fit))^2)/(sum(data$pweights)),
                        "1/mu^2"=   sum(data$pweights*(1/(Y^2)-1/(fitted(fit))^2)^2)/(sum(data$pweights)),                    
                        )), sdev)
        shape <- ifelse(shape.miss,1/summary(fit)$dispersion, shape)                   
  } else  {
        sdev  <- 1
        shape <- 0
  } 



  
  # Initial disparity (-2logL)
  ML.dev0 <- -2*sum(data$pweights*switch(family$family,
                "gaussian"= dnorm(fit$y, fitted(fit), sdev, log=TRUE),
                "poisson" = dpois(fit$y, fitted(fit), log=TRUE),
                "binomial"= dbinom(Y[,1],Y[,1]+Y[,2], fitted(fit),log=TRUE),
                "Gamma"   = dgamma(fit$y,  shape=shape, scale=fitted(fit)/shape, log=TRUE),
                "inverse.gaussian"= dinvgauss(fit$y,fitted(fit), shape=shape, log=TRUE),              
                ))


  
  # Return (glm) output and terminate if k=1  
  if (k == 1) {
      if (random.distribution=="np"){
            names(fit$coefficients) <- ifelse(names(fit$coefficients)=="(Intercept)", "MASS1", names(fit$coefficients))
      }
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
              rsdev = 0,
              post.prob = list(matrix(1,N,1,dimnames=list(names0,"") )),
              post.int =  list(fit$coef[1]),
              ebp = list(family$linkfun(fit$fitted)),
              EMiter = 0,
              EMconverged = "none",
              lastglm = list(fit),
              Misc = list(list(lambda=lambda))
              )          
      if (random.distribution =="np"){
              class(fit) <- 'glmmNPML'
      } else {
             class(fit) <- 'glmmGQ'
      }
      return(fit)
  } else if (!(k %in% 1:600)){
      stop("This choice of k is not supported.")
  }

  # Omit integration point if GH weights are too small  
  tmp       <- gqz(k, minweight=1e-50)  # from version 0.39-1; changed to 1e-50 in v0.42
  k0        <- k  # from 0.42: save for glmmGQ output
  k         <- min(k, dim(tmp)[1])
  
  # Expand the data
  if(family$family=="binomial"){
      YP    <- binomial.expand(Y,k,rep(1,N)); Y<- YP[[1]]; PY<-YP[[2]];r<-YP[[3]]; n<-YP[[4]]
  }  else { 
      Y     <- rep(Y,k)
  }          
  datak     <- expand(data,k)
  kindex    <- rep(1:k,rep(N,k))# index for the mixtures, for version 0.31 or higher only used for allvc 
  #tmp       <- gqz(k,minweight=1e-55)  # omitted from version 0.39-1
  z0        <- -tmp$location
  z         <- rep(-tmp$location,rep(N,k))
  p         <- tmp$weight
  offset    <- datak$offset
  pweights  <- datak$pweights

  # Generate the random design matrix and append to fixed matrix
  if (random.distribution=='np'){
      # Nonparametric random effect
      X <- model.matrix(formula,datak)[,-1,drop=FALSE]
      datak$MASS <- gl(k,N) 
      if (mform=='1') {
        random <- formula(~MASS-1)
      } else {  
          # Nonparametric random coefficient
          if (!int.removed){
             random <- formula(paste('~ MASS + ',paste(mform,'MASS',sep=":",collapse='+'), '-1',sep=''))
         } else{
             X <- model.matrix(formula,datak)[,,drop=FALSE]
             random <- formula(paste('~ ',paste(mform,'MASS',sep=":",collapse='+'), '-1',sep=''))
          }  
      }
  } else {
      X <- model.matrix(formula,datak)
      #if (mform=='1') {
      random <- formula(~ z - 1 )
      #} else { 
      #    stop("Random coefficients only supported with option random.distribution='np'")
      #} # 05/11/07
  }
  Z <- model.matrix(random,datak)
  if (dim(X)[1]!= dim(Z)[1]){ cat("The missing value routine cannot cope with this model. Please specify the random term also as fixed term and try again. " )}
  XZ <- cbind(X,Z)
  
  # Extend the linear predictor:
  if (missing(pluginz)){
      sz <- tol* sdev*z
  } else {
      if (length(pluginz)!=k){
        stop("pluginz needs to be a vector of length k.")  #30/09/09
      } else {
        sz <- rep(pluginz-fit$coef[[1]],rep(N,k))
      }
  }      
  
  Eta <- fit$linear.predictor + sz
           # The extra term stops unrelated regressions
                                        
  # Initial EM trajectory values 
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

  
  # Expanded fitted values     
  Mu <- family$linkinv(Eta) 
  
  if (sum(is.na(Mu))>0){
     if (family$link=="1/mu^2"){                                                                      
        warning("The squared reciprocal link will often fail. Try family=...(link=log) instead.")
     }
     stop("Unable to transform extended linear predictor to response scale.")
   }  
  
  # Calculate loglikelihood for fixed model
  f <- switch(family$family,
              "gaussian" = dnorm(Y,Mu,tol*sdev,log=TRUE),
              "poisson"  = dpois(Y,Mu,log=TRUE),
              "binomial" = dbinom(r,n,Mu,log=TRUE),
              "Gamma"    = dgamma(Y,shape=shape/tol^2,scale=Mu*tol^2/shape ,log=TRUE),
              "inverse.gaussian"= dinvgauss(Y, Mu, shape= shape/tol^2, log=TRUE),
              )
   

  # Calculate the weights from initial model
  tmp <-weightslogl.calc.w(p,matrix(f,ncol=k),pweights[1:N])
  w <- tmp$w
  
  # Initialize for EM loop
  ML.dev <- ML.dev0
  iter <- ml<- 1
  converged <- FALSE
  sdevk<-rep(sdev,k);  shapek<-rep(shape,k)


  
  ##########Start of EM ##########
  while (iter <= EMmaxit && (!converged || (iter<=9 && random.distribution=='np' && damp && (family$family=="gaussian" && sdev.miss || (family$family=="Gamma"|| family$family=="inverse.gaussian")  && shape.miss)  ))){
      if (verbose){ cat(iter,"..") }      
       
      # M-Step: Weighted GLM
      fit <- try(glm.fit(x=XZ, y=Y, weights = as.vector(w)*pweights, family = family, offset=offset, ...), silent=TRUE)
      if (class(fit)=="try-error"){
              stop("Singularity or Likelihood-Spike at iteration #",iter,". 
              Check model specification,  enable spike protection or smooth among components.")
      }

      # EM Trajectories   
      if (random.distribution=="np" && !int.removed ){
          masspoint<- fit$coef[l0:(l0+k-1)]
          followmass<-rbind(followmass, masspoint)
      }
   
      # Fitted response from current model
      Mu <- fitted(fit)
      
      # Unequal component dispersion parameters  
      if (family$family=="gaussian"){
          if (sdev.miss){ sdev<- sqrt(sum((as.vector(w)*pweights)*(Y-Mu)^2)/sum(as.vector(w)*pweights))}
          sdevk<-rep(sdev,k) 
          if (lambda!=0){
              for (l in 1:k){
                wk<-matrix(1,k,N); wk[1:k,]<-dkern(1:k,l,k,lambda);wk<-t(wk)
                sdevk[l] <-  sqrt(sum(wk* as.vector(w)*pweights *(Y-Mu)^2)/sum(wk*as.vector(w)*pweights))
              }
              sk<-rep(sdevk,rep(N,k))
          } else {
              sk<-sdev
          }
      }  else {
          sdevk <-rep(NA,k)
      }
      if (family$family=="Gamma"){
           if (shape.miss) { shape<-(sum(as.vector(w)*pweights))*1/sum(as.vector(w)*pweights*((Y-fitted(fit))/fitted(fit))^2)}
           shapek<-rep(shape,k) 
           if (lambda!=0){
                for (l in 1:k){
                  wk<-matrix(1,k,N); wk[1:k,]<-dkern(1:k,l,k,lambda);wk<-t(wk)
                  shapek[l] <- sum(wk*as.vector(w)*pweights)/ sum(wk* as.vector(w)*pweights*((Y-Mu)/Mu)^2)
                }
                shk<-rep(shapek,rep(N,k))
           } else {
                shk<-shape
           }
      } else {
           shapek<-rep(NA,k)
      }
      

     if (family$family=="inverse.gaussian"){
           if (shape.miss) { shape<-(sum(as.vector(w)*pweights))*1/sum(as.vector(w)*pweights*(Y-fitted(fit))^2/(fitted(fit))^3)}
           shapek<-rep(shape,k) 
           if (lambda!=0){
                for (l in 1:k){
                  wk<-matrix(1,k,N); wk[1:k,]<-dkern(1:k,l,k,lambda);wk<-t(wk)
                  shapek[l] <- sum(wk*as.vector(w)*pweights)/ sum(wk* as.vector(w)*pweights*(Y-Mu)^2/Mu^3)
                }
                shk<-rep(shapek,rep(N,k))
           } else {
                shk<-shape
           }
      } else {
           shapek<-rep(NA,k)
      }

      
      
      # Calculate loglikelihood for expanded model for this iteration
      f <- matrix(switch(family$family,
              "gaussian"=dnorm(Y,Mu,(1-(1-tol)^(damp.power*iter+1))*sk,log=TRUE),
              "poisson" =dpois(Y,Mu,log=TRUE),
              "binomial"=dbinom(r,n,Mu,log=TRUE),
              "Gamma"=dgamma(Y,shape=shk/(1-(1-tol)^(damp.power*iter+1))^2,scale=Mu*(1-(1-tol)^(damp.power*iter+1))^2/shk,log=TRUE),
              "inverse.gaussian"=dinvgauss(Y, Mu, shape= shk/(1-(1-tol)^(damp.power*iter+1))^2, log=TRUE),
           ),       
           nrow=N, ncol=k)
              
      # Calculate the component proportions from the weights
      if (random.distribution=='np') {
          p <- as.vector(apply(w*pweights,2,sum))/sum(pweights[1:N]) #16-03-06
      }
             
      # E-Step: Update weights 
      tmp <- weightslogl.calc.w(p,f,pweights[1:N])
      w   <- tmp$w
      
      # Update disparity and check for convergence
      ML.dev[iter+1] <- ifelse(is.na(tmp$ML.dev), Inf, tmp$ML.dev)
      if (ML.dev[iter+1] > ML.dev0) {ml<-ml+1}  #only relevant for graphical output
      converged <- abs(ML.dev[iter+1] - ML.dev[iter])< EMdev.change
      iter <- iter + 1
      
      # Likelihood Spike Protection
      if (random.distribution != 'gq' && spike.protect!=0){
          if (family$family=='gaussian' && abs(min(sdevk/masspoint)) < 0.000001*spike.protect){break}  
          if ((family$family=='Gamma'||family$family=='inverse.gaussian')  && abs(max(shapek/masspoint)) > 10^6*spike.protect){
            break
          }
      }  
  }###########################End of EM loop#############

   
  # Print on screen information on EM convergence
  if (verbose){
      cat("\n")
      if (converged){
          cat("EM algorithm met convergence criteria at iteration # ", iter-1,"\n")
      } else {
        cat("EM algorithm failed to meet convergence criteria at iteration # ",
        iter-1,"\n")  
      }
  }

  
  # Compute model deviance
  Deviance <- switch(family$family,
              "gaussian"= sdev^2*ML.dev[iter]-sdev^2* log(2*pi*sdev^2)* sum(data$pweights[1:N]),
              "poisson" = ML.dev[iter] +2*sum(data$pweights[1:N]*(-Y[1:N]+Y[1:N]*log(Y[1:N]+ (Y[1:N]==0))-lfactorial(Y[1:N]))),
              "binomial"= ML.dev[iter] +2*sum(data$pweights[1:N]*(lfactorial(n)-lfactorial(r)-lfactorial(n-r) - n*log(n) + r*log(r+(r==0))+(n-r)*log(n-r+((n-r)==0)))[1:N]),
              "Gamma"   = 1/shape*ML.dev[iter]+2/shape*(sum(data$pweights[1:N])*shape*(log(shape)-1)-sum(data$pweights[1:N])*lgamma(shape)-sum(data$pweights[1:N]*log(Y[1:N]))),
              "inverse.gaussian"=   1/shape*( ML.dev[iter]+ sum(data$pweights[1:N])*log(shape/(2*pi)) - 3*sum(data$pweights[1:N]*log(Y[1:N])) ),
              )

  # Compute  posterior prob. etc. 
  mass.points   <- masses <- NULL
  np            <- length(fit$coef)
  ebp           <- apply(w*matrix(fit$linear.predictor,N,k,byrow=FALSE),1,sum)  # Emp. Bayes Pred. (Aitkin, 96)
  ebp.fitted    <- family$linkinv(ebp)
  ebp.residuals <- Y[1:N]- ebp.fitted
  names(ebp)    <- names(ebp.fitted) <- names(ebp.residuals) <- names0
  if (mform %in% substring(names(fit$coef),1, nchar(mform)) &&!int.removed ){length(fit$coefficients) <- np <- np-1}# if one variable is random *and* fixed 
  # if (is.na(fit$coefficients[np])){length(fit$coefficients)<-np<-np-1}# replaced by the line above from 0.42 on
  m <- seq(1,np)[substr(attr(fit$coefficients,'names'),1,4)=='MASS']
  if (random.distribution=="np"){
      mass.points   <-  fit$coefficients[m] # from 0.42
   #   if (int.removed){          
   #       print(length(fit$coef))
   #       fit$coef<-fit$coef[-length(fit$coef)]
   #       print(fit$coef)
   #    }    
  } else {  
      a <- ifelse(names(fit$coef[1])== "(Intercept)", fit$coef[1], 0) #02-08-06
      mass.points <- a + fit$coef["z"]*z0           # from 0.42, np replaced by "z"
    }
  post.prob     <- matrix(w, nrow=N, byrow=FALSE, dimnames=list(names0, 1:k) )
  post.int      <- as.vector(post.prob %*% mass.points[1:k]); names(post.int) <- names0

  # Write tol values as plot title if alldist is called from tolfind:
  if ((plot.opt==1 || plot.opt==2) && par("mfrow")[1]>2) { 
        plot.main <- substitute("tol"== tol0, list(tol0=tol0))
  } else {
        plot.main <- c("")
  }

  # Set up graphics device and plot disparity trend 
  if (plot.opt==3 && random.distribution=="np"){
      par(mfrow=c(2,1), cex=0.5, cex.axis=1.5, cex.lab=1.5)
  }
  if (plot.opt==1|| plot.opt==3){
      if  ((family$family=="gaussian" && sdev.miss|| (family$family=="Gamma" ||family$family=="inverse.gaussian") && shape.miss) && damp  && random.distribution=='np' && iter>=max(8,ml+1)){
          # Linear interpolation for initial cycles
          ML.dev[2: max(7,ml)]<-ML.dev0+ 1:max(6,ml-1)/ max(7,ml)*(ML.dev[max(8,ml+1)]-ML.dev0) 
      }  
      plot(0:(iter-1),ML.dev, col=1,type="l",xlab='EM iterations',ylab='-2logL', main= plot.main )  
      if (verbose){ cat("Disparity trend plotted.\n")}
  }

  # Prepare output for glmmNPML objects
  if (random.distribution=="np") {
      
      # Mixture proportions 
      masses <- as.vector(apply(w*pweights,2,sum))/sum(pweights[1:N]) # 16-03-05
      names(masses) <- paste('MASS',1:k,sep='')
       
      # Estimate random effect standard deviation                  # from 0.42
      rsdev <- sqrt(sum(masses * (mass.points[1:length(masses)] - sum(masses*mass.points[1:length(masses)]) )^2))

      # Compute fixed part residuals
      if (family$family=="binomial"){
          R0 <- family$linkfun(PY[1:N])
      } else  {
          R0 <- family$linkfun(Y[1:N])
      }
      if(dim(X)[2]>0){
          R <- R0 - X[1:N,]%*%matrix(fit$coef[1:dim(X)[2]])-offset[1:N]
      } else {
          R <- R0 - offset[1:N]
      }
      R <-as.vector(R);  names(R) <- names0  
      
      # EM trajectory plot  
      if (mform=='1' && any(is.finite(R))){
           ylim <- c(min(R[is.finite(R)]), max(R[is.finite(R)]))  #29/06/2006
      } else  { 
            ylim <- c(min(followmass[,]),max(followmass[,]))
      }
      if(any(is.na(ylim)) &  plot.opt >1 ){
              cat("Singularity: EM Trajectory plot not available.", "\n");
              plot.opt<-min(plot.opt,1)
      }
      if ((plot.opt==2|| plot.opt==3)&&!int.removed ){
            plot(0:(iter-1),followmass[,1],col=1,type='l',ylim=ylim,ylab='mass points',xlab='EM iterations',  main=plot.main )
            for (i in 1:k){ lines(0:(iter-1), followmass[,i],col=i)
                        if (mform=='1'){ points(rep(iter-1,length(R)),R)}}
            if (verbose){ cat("EM Trajectories plotted.\n")}
      }
      
      # glmmNPML output    
      fit <- c( fit[1],
                residuals = list(ebp.residuals),
                fitted.values = list(ebp.fitted),
                fit[c(8,9)],
                disparity = ML.dev[iter],
                deviance = Deviance,
                fit[12],  
                df.residual = N-np-k+1,
                df.null = N-1, 
                fit[18],
                call = call,
                formula = formula,
                random = rform,
                data = list(subset(data,select=1:ddim[2])),
                model = list(XZ),
                weights = list(w0),
                offset = list(off0),
                mass.points = list(mass.points),
                masses = list(masses),               
                sdev = list(list(sdev=sdev, sdevk=sdevk)),
                shape = list(list(shape=shape,shapek=shapek)),
                rsdev = list(rsdev),
                post.prob = list(post.prob),
                post.int = list(post.int),
                ebp = list(ebp),
                EMiter = iter - 1,
                EMconverged = converged,
                lastglm = list(fit),
                Misc = list(list(Disparity.trend=ML.dev,EMTrajectories=followmass,res=R,ylim=ylim, lambda=lambda,mform=mform))
            )
      class(fit) <- 'glmmNPML'
  
  } else {
      # glmmGQ output
      fit <- c( fit[1],
                residuals = list(ebp.residuals),
                fitted.values = list(ebp.fitted),
                fit[c(8,9)],
                disparity = ML.dev[iter],
                deviance = Deviance,
                fit[12],  
                df.residual = N - np,
                df.null = N - 1, 
                fit[18],
                call = call,
                formula = formula,
                random = rform,
                data = list(subset(data,select=1:ddim[2])),
                model = list(XZ),
                weights = list(w0),
                offset = list(off0), 
                mass.points = list(mass.points),
                masses = list(gqz(k0, minweight=1e-50)$weight),          
                sdev = list(list(sdev=sdev, sdevk=sdevk)),
                shape = list(list(shape=shape,shapek=shapek)),
                rsdev = fit$coef[["z"]],
                post.prob = list(post.prob),
                post.int =  list(post.int),
                ebp = list(ebp),
                EMiter = iter - 1,
                EMconverged = converged,
                lastglm = list(fit),
                Misc = list(list(Disparity.trend=ML.dev, lambda=lambda, mform=mform))
           )
      class(fit) <- 'glmmGQ'       
  }
  fit
}

