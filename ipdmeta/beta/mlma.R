mlma <- function(
   ipd.formula,
   meta.formula,
   random.formula,
   ipd.data,
   meta.data,
   ipd.groups,
   meta.groups,
   sigma2,
   study.group.interaction,
   max.iter=5,
   min.sample=300,
   est.delta=.1,
   mc.step=2,
   df=30,
   echo = TRUE,
   max.sample = NULL,
   converge.index,
   error=1/100,
   fixed=FALSE
){

   ### MISSING REMOVAL

   vars <- all.vars(ipd.formula)
   if(!fixed) vars <- c(vars,all.vars(random.formula))
   
   missing.cases <- sapply(vars,function(x){any(is.na(ipd.data[,x]))})

   if(any(missing.cases)){
     i <- which(missing.cases)
     missing.index <- sapply(vars[i],function(x){which(is.na(ipd.data[,x]))})
     missing.index <- unique(as.vector(unlist(missing.index)))
     ipd.data <- ipd.data[-missing.index,]
     print(paste(length(missing.index),"Missing cases excluded"))
   }
   
if(fixed){

fit <- coxmeta.fixed(
              ipd.formula=ipd.formula,
              meta.formula=meta.formula,
              ipd.data=ipd.data,
              meta.data=meta.data,
              sigma2=sigma2,
              study.group=study.group.interaction,
              error=error)

 }
else{
fit <- coxmeta.mixed(
              ipd.formula,
              meta.formula,
              random.formula,
              ipd.data,
              meta.data,
              ipd.groups,
              meta.groups,
              sigma2,
              study.group.interaction,
              max.iter,
              min.sample,
              est.delta,
              mc.step,
              df,
              echo,
              max.sample,
              converge.index,
              error)
 }

return(fit)

}

###DEPENDENCIES

coxmeta.fixed <- function(
              ipd.formula,
              meta.formula,
              ipd.data,
              meta.data,
              sigma2,
              study.group.interaction,
              error=1/1000)
{

    beta.index <- beta.indices(ipd.formula,meta.formula)
    beta.ad = runif(length(beta.index$ad))
    beta.ipd = runif(max(beta.index$ipd))

    #KNOWN COVARIANCE-VARIANCE FOR SURVIVAL OUTCOMES
    s <- meta.data[,all.vars(meta.formula)[[1]]]
    Sigma <- surv.covariance.variance(s,sigma2,study.group.interaction)

####NEWTON-RAPHSON MAXIMIZATION

    converge <- error*2
    loglik <- 0
    last <- FALSE

    monitor <- list(coef=c(),loglik=c())

    while(converge>error|!last){
    
    if(converge<error) last = TRUE

    #PATIENT LEVEL SCORE/INFO/LIKELIHOOD
    NR.cox <- coxph.with.offset(ipd.formula,ipd.data,beta=beta.ipd)        

    #STUDY LEVEL SCORE/INFO/LIKELIHOOD
    NR.surv <- surv.with.offset(meta.formula,meta.data,beta=beta.ad,Sigma)

    #COMBINED
    NR <- combined.score.info(
          beta.index,
	  U.cox=NR.cox$score,
	  I.cox=NR.cox$info,
	  loglik.cox=NR.cox$loglik,
	  U.surv=NR.surv$score,
	  I.surv=NR.surv$info,
	  loglik.surv=NR.surv$loglik)

    #UPDATE
    update <- as.vector(solve(NR$info)%*%NR$score)
 
    beta.ipd <- beta.ipd+update[beta.index$ipd]
    beta.ad <- beta.ad+update[beta.index$ad]
    beta <- beta.construct(beta.ipd,beta.ad,beta.index$shared)

    converge = abs(NR$loglik-loglik)
    loglik = NR$loglik
  
    monitor$coef <- cbind(monitor$coef,beta)
    monitor$loglik <- c(monitor$loglik,loglik)
    }

    beta <- array(beta)
    names(beta) <- beta.index$labels

    return(
           list(
           coef=beta,
           var=solve(NR$info),
           loglik=loglik,
           monitor=monitor)
           )
}


combined.score.info <- function(
                    beta.index,
                    U.cox,
                    I.cox,
                    loglik.cox,
                    U.surv,
                    I.surv,
                    loglik.surv)
{

   p <- max(beta.index$ad)
   U <- rep(0,p)
   U[beta.index$ipd] <- U.cox
   U[beta.index$ad] <- U[beta.index$ad]+U.surv

   #COMBINING INFO
   I <- matrix(0,p,p)
   I[beta.index$ipd,beta.index$ipd] <- I.cox
   I[beta.index$ad,beta.index$ad] <- I[beta.index$ad,beta.index$ad]+I.surv

   loglik <- loglik.cox+loglik.surv

   return(
           list(
           score=U,
           info=I,
           loglik=loglik)
           )
}


coxph.with.offset <- function(
                  formula,
                  data,
                  beta,
                  offset)
{

					#FORMULA CONSTRUCTION
   names <- all.vars(terms(formula))
   time <- data[,names[1]]
   event <- data[,names[2]]

   S <- Surv(time,event)
   covariates <- paste(attr(terms(formula),"term.labels"),collapse="+")

   f <- as.formula(paste(c("S~",covariates),c=""))

   if(!missing(offset)) f <- as.formula(paste(c(f,"offset(offset)"),
							collapse="+"))

  
					#MODEL ESTIMATION AT BETA/OFFSET
   fit <- coxph(f,data=data,init=beta,iter=0)
   fit.detail <- coxph.detail(fit)

					#SCORE/INFO/LOGLIK

   p <- length(beta)

   U <- fit.detail$score
   if(is.matrix(U)) U <- colSums(U) else U <- sum(U)

   I <- matrix(as.vector(fit.detail$imat),nrow=p^2)
   I <- apply(I,1,sum)   
   I <- matrix(I,p,p)

   return(
          list(
          score=U,
          info=I,
          loglik=fit$loglik[2])
          )
}


surv.with.offset <- function(
                 formula,
                 data,
                 beta,
                 Sigma,
                 offset)
{

   X <- model.matrix(terms(formula),data)

   S <- data[,all.vars(terms(formula))[1]]
   
   R <- log(-log(S))-X%*%beta
   
   if(!missing(offset)) R <- R-offset

   U <- as.vector(t(X)%*%solve(Sigma)%*%R)
   I <- as.vector(t(X)%*%solve(Sigma)%*%X)
   I <- matrix(I,length(beta),length(beta))

   loglik <- -1/2*t(R)%*%solve(Sigma)%*%R

   return(
          list(
          score=U,
          info=I,
          loglik=loglik)
          )
}


coxmeta.mixed <- function(
              ipd.formula,
              meta.formula,
              random.formula,
              ipd.data,
              meta.data,
              ipd.groups,
              meta.groups,
              sigma2,
              study.group.interaction,
              max.iter,
              min.sample,
              est.delta,
              mc.step,
              df,
              echo = TRUE,
              max.sample=NULL,
              converge.index,
              error)
{

###INITIALIZATION
       
       beta.index <- beta.indices(ipd.formula,meta.formula)

	coxmeta.objects <- 
	   coxmeta.initialize(
	      ipd.formula,
   	      meta.formula,
   	      random.formula,
   	      ipd.data,
   	      meta.data,
   	      ipd.groups,
   	      meta.groups,
   	      sigma2,
   	      study.group.interaction,
              beta.index
	   )

   D = coxmeta.objects$params$D

   fixed.coxmeta <- coxmeta.fixed(
 	      ipd.formula,
   	      meta.formula,
   	      ipd.data,
   	      meta.data,
   	      sigma2,
   	      study.group.interaction,
              error=error
	  )

   D = coxmeta.objects$params$D
   beta.ipd = fixed.coxmeta$coef[beta.index$ipd]
   beta.ad = fixed.coxmeta$coef[beta.index$ad]
   loglik.null = fixed.coxmeta$loglik
   p.D = dim(D)[1]

   if(missing(converge.index)) converge.index = 1:(max(beta.index$ad)+p.D)

###PROPOSAL SAMPLING

    rmt.cox <- function(n){rmt(n,coxmeta.objects$params$cox.loc,
				coxmeta.objects$params$cox.scale,df=df)}

    rmt.surv <- function(n){rmt(n,coxmeta.objects$params$surv.loc,
				coxmeta.objects$params$surv.scale,df=df)}

    weights <- function(B.cox,B.surv){
       
        coxmeta.importance.weights(
	coxmeta.objects$model[[1]],coxmeta.objects$model[[2]],
	B.cox,B.surv,D,beta.cox=beta.ipd,beta.surv=beta.ad,
   	coxmeta.objects$params$cox.loc,coxmeta.objects$params$cox.scale,
	coxmeta.objects$params$surv.loc,coxmeta.objects$params$surv.scale,df)

}

#CONVERGENCE/SAMPLE SIZE INITIALIZATION
      N = min.sample
      if(!is.null(max.sample)) fixed.escalation = seq(min.sample,max.sample,length=max.iter)

      n.iter <- 0
      converged = FALSE

#COXMCEM OUTPUT

      coxmcem.object <- 
                     list(
                     max.weight=c(),
                     mc.samples=c(),
                     est.converge=c(),
                     sd.loglik=c(),
                     loglik=c()
                     )
   
while(!converged&n.iter<max.iter){

   n.iter <- n.iter+1

#SAMPLING RANDOM EFFECTS

   B.cox <- rmt.cox(N)
   B.surv <- rmt.surv(N)

#IMPORTANCE WEIGHTS

   weight.object <- weights(B.cox,B.surv)
   w <- weight.object$weights

#ESTEP

   Estep <- 
      coxmeta.Estep(coxmeta.objects$model[[1]],coxmeta.objects$model[[2]],
		B.cox,B.surv,weights=w)

#MSTEP

    Mstep <- coxmeta.Mstep(beta.index,beta.ipd,beta.ad,
		coxmeta.objects$model[[1]],
		coxmeta.objects$model[[2]],
		Estep)


###CONVERGENCE CRITERIA
###RELATIVE CHANGE

    beta = beta.construct(beta.ipd,beta.ad,beta.index$shared)
    zeta.now = c(Mstep$fit$coef,diag(Mstep$D))
    zeta.last = c(beta,get.diag(D))

    est.converge <-  relative.criteria(zeta.now[converge.index],
                                       zeta.last[converge.index])

###UPDATE COX OBJECT

    coxmcem.object$max.weight <- append(coxmcem.object$max.weight,max(w))

    tries <- length(coxmcem.object$max.weight)
    m <- length(coxmcem.object$est.converge)

    if(tries>3) cv.last <- sd(coxmcem.object$est.converge[(m-2):m])/mean(coxmcem.object$est.converge[(m-2):m])

    coxmcem.object$est.converge <- append(coxmcem.object$est.converge,
                                          est.converge)

    if(tries>3) cv.now <- sd(coxmcem.object$est.converge[(m-1):(m+1)])/mean(coxmcem.object$est.converge[(m-1):(m+1)])


    coxmcem.object$mc.samples <- append(coxmcem.object$mc.samples,N)
    coxmcem.object$sd.loglik <- append(coxmcem.object$sd.loglik,
                                       weight.object$sd.loglik)
    coxmcem.object$loglik <- append(coxmcem.object$loglik,
                                    weight.object$loglik)

###CHECK IF CONVERGENCE IS MET; UPDATE SAMPLE SIZE IF NOT


     if(tries>2){

        est.convergence <- max(coxmcem.object$est.converge[(m-1):(m+1)])
  
          if(est.convergence<est.delta){
             converged = TRUE
	          }
              }

     if(is.null(max.sample)){

	    if(tries>3){
		if(cv.now>cv.last) N = ceiling(N+N/mc.step)
       	       	 }
  	    }
     else{
	  N = ceiling(fixed.escalation[n.iter])
         }

#UPDATE ESTIMATES; PRINT CONVERGENCE
      if(tries>3&echo) print(paste(c("N","Max Weight","% Relative Change"),
         round(c(coxmcem.object$mc.samples[length(coxmcem.object$mc.samples)],
	 max(w),est.converge*100),5)))

###UPDATE

    beta.ipd = Mstep$fit$coef[beta.index$ipd]
    beta.ad = Mstep$fit$coef[beta.index$ad]
    D = Mstep$D

  }

#FINAL ESTIMATES

   coxmcem.object$iterations <- n.iter
   coxmcem.object$coef <- beta.construct(beta.ipd,beta.ad,beta.index$shared)
   coxmcem.object$vcov <- D
   coxmcem.object$cluster <- rbind(Estep$cox.mean,Estep$surv.mean)

#VARIANCE ESTIMATION
   coxmcem.object$var <- list()

   coxmcem.object$var$coef <- 
  
 	  coxmeta.info.coef(beta.index,beta.ipd,beta.ad,
		cox=coxmeta.objects$model[[1]],
		surv=coxmeta.objects$model[[2]],
		B.cox,B.surv,w)
   
   coxmcem.object$var$coef <- solve(coxmcem.object$var$coef)

   coxmcem.object$var$vcov <- 
   	 coxmeta.info.vcov(B.cox,B.surv,ipd.groups,meta.groups,D,w)

   coxmcem.object$var$vcov <- 1/coxmcem.object$var$vcov

return(coxmcem.object)
}


coxmeta.Estep <- function(
              cox,
              surv,
              B.cox,
              B.surv,
              weights)
{

   #Patient-level offsets

   u <- apply(B.cox,1,function(b){
      cox$Z%*%b
   })
   
   #RETURNS THE SAMPLED OFFSETS BY COLUMN
   #OBTAIN EACH ROWS WEIGHTED AVERAGE
   u <- u%*%weights

   #STUDY-LEVEL OFFSETS
   v <- apply(B.surv,1,function(b){
      surv$Z%*%b
   })

   v <- v%*%weights
   
   #WEIGHTED MEANS
   cox.means <- t(B.cox)%*%weights
   surv.means <- t(B.surv)%*%weights

   #WEIGHTED SECOND MOMENT
   #RETURNS A VECTOR OF THE STACKED SUBVECTORS

   V.cox <- apply(B.cox,1,function(b){
      b <- cluster.subvector.matrix(b,n.groups=cox$n.groups)
      apply(b,1,function(x){outer(x,x)})
   })

   V.cox <- V.cox%*%weights


   V.surv <- apply(B.surv,1,function(b){
      b <- cluster.subvector.matrix(b,n.groups=surv$n.groups)
      apply(b,1,function(x){outer(x,x)})
   })

   V.surv <- V.surv%*%weights

    return(
          list(
	  cox.offset=u,
          surv.offset=v,
          cox.mean=cox.means,
          surv.mean=surv.means,
          cox.second=V.cox,
          surv.second=V.surv)
          )
}
 


coxmeta.Mstep <-  function(
              beta.index,
              beta.ipd,
              beta.ad,
              cox,
              surv,
              Estep)
{

    #BETA PARAMETERS
    fit <- coxmeta.offset(beta.index,beta.ipd,
                          beta.ad,cox,surv,
                          Estep$cox.offset,Estep$surv.offset)

    #VARIANCE PARAMETERS
    D <- cluster.variance(cox,surv,Estep$cox.second,Estep$surv.second)

return(list(fit=fit,D=D))
}

#MAXIMIZATION OF VARIANCE COMPONENTS

cluster.variance <- function(
                  cox,
                  surv,
                  cox.second,
                  surv.second)
{
   ncox <- length(cox.second)
   nsurv <- length(surv.second)
   D.cox <- cluster.subvector.matrix(cox.second,ncox/cox$n.groups)
   D.surv <- cluster.subvector.matrix(cox.second,nsurv/surv$n.groups)

   D <- cbind(D.cox,D.surv)
   D <- apply(D,1,mean)
   Ddim <- sqrt(length(D))
   
   matrix(D,Ddim,Ddim)
}

#MAXIMIZATION OF EFFECTS

coxmeta.offset <-  function(
               beta.index,
               beta.ipd,
               beta.ad,
               cox,
               surv,
               cox.offset,
               surv.offset,
               error=1/1000)
{


    #KNOWN COVARIANCE-VARIANCE FOR SURVIVAL OUTCOMES
    Sigma <- surv$Sigma

####NEWTON-RAPHSON MAXIMIZATION

    converge <- error*2
    loglik <- 0
    last <- FALSE

    monitor <- list(coef=c(),loglik=c())

    while(converge>error|!last){
    
    if(converge<error) last = TRUE

    #PATIENT LEVEL SCORE/INFO/LIKELIHOOD
    NR.cox <- coxph.with.offset(cox$f,cox$data,beta=beta.ipd,cox.offset)        

    #STUDY LEVEL SCORE/INFO/LIKELIHOOD
    NR.surv <- surv.with.offset(surv$f,surv$data,beta=beta.ad,Sigma,surv.offset)

    #COMBINED
    NR <- combined.score.info(
          beta.index,
	  U.cox=NR.cox$score,
	  I.cox=NR.cox$info,
	  loglik.cox=NR.cox$loglik,
	  U.surv=NR.surv$score,
	  I.surv=NR.surv$info,
	  loglik.surv=NR.surv$loglik)

    #UPDATE
    update <- as.vector(solve(NR$info)%*%NR$score)
 
    beta.ipd <- beta.ipd+update[beta.index$ipd]
    beta.ad <- beta.ad+update[beta.index$ad]
    beta <- beta.construct(beta.ipd,beta.ad,beta.index$shared)

    converge = abs(NR$loglik-loglik)
    loglik = NR$loglik
  
    monitor$coef <- cbind(monitor$coef,beta)
    monitor$loglik <- c(monitor$loglik,loglik)
    }

return(
        list(
	coef=beta,
	var=solve(NR$info),
	loglik=loglik,
	monitor=monitor)
        )
}

coxmeta.initialize <- function(

                   ipd.formula,
                   meta.formula,
                   random.formula,
                   ipd.data,
                   meta.data,
                   ipd.groups,
                   meta.groups,
                   sigma2,
                   study.group.interaction,
                   beta.index)
{

	model.objects <-  coxmeta.data.objects(
	   ipd.formula,
  	   meta.formula,
   	   random.formula,
   	   ipd.data,
   	   meta.data,
   	   ipd.groups,
   	   meta.groups,
   	   sigma2,
   	   study.group.interaction
	   )

        parameters <-  parameter.initialize(
	        random.formula,
                model.objects[[1]],
                model.objects[[2]],
                beta.index)

  return(
         list(
          model=model.objects,
          params=parameters)
          )
}

parameter.initialize <- function(
                     random.formula,
                     cox.data.object,
                     surv.data.object,
                     beta.index)
{

       f.random <- as.character(random.formula)[2]
       phmm.fixed <- sub("(.*\\| ?)(.*)","cluster(\\2",f.random)
       phmm.fixed <- paste(c(cox.data.object$f,phmm.fixed),collapse="+")
       
       phmm.random <- sub("\\((.*)\\| ?(.*)","~\\1",f.random)

       coxme.formula <- paste(
       		     c(cox.data.object$f,as.character(random.formula)[2]),
		     collapse="+")

        coxme.formula <- formula(coxme.formula)
        phmm.formula <- list(fixed=formula(phmm.fixed),
       		    random=formula(phmm.random))
 
 
   if(cox.data.object$n.groups>1){

       cox.params <- coxmcem.mixed.initialize(
       		  coxme.formula,
		  phmm.formula,
		  cox.data.object$data
		  )

       
       D = as.matrix(cox.params$D)
       beta.ipd = cox.params$beta

       cox.loc <- cox.params$mean
       cox.scale <- cox.params$Sigma

         }
     else{
       D = make.diagonal.matrix(runif(ncol(cox.data.object$Z)))
       beta.ipd = rnorm(length(beta.index$ipd))
       }
  
       beta.ad = runif(length(beta.index$ad))
       beta.ad[beta.index$shared.index] <- beta.ipd[beta.index$ad[beta.index$shared.index]]

       beta = beta.construct(beta.ipd,beta.ad,beta.index$shared.index)

       DA <- diag(rep(diag(D),each=surv.data.object$n.groups))

       scale <- solve(solve(DA)+t(surv.data.object$Z)%*%solve(surv.data.object$Sigma)%*%surv.data.object$Z)
       R <- log(-log(surv.data.object$S))-surv.data.object$X%*%beta.ad
       loc <- scale%*%t(surv.data.object$Z)%*%solve(surv.data.object$Sigma)%*%R
       loc <- as.numeric(loc)

       n.frailties = length(loc)
 
       N = surv.data.object$n.groups*(nrow(cox.data.object$data)/cox.data.object$n.groups)

       scale = diag(rep(n.frailties/N,n.frailties))

       if(cox.data.object$n.groups==1){

         cox.loc <- cluster.subvector.matrix(loc,surv.data.object$n.groups)
         cox.loc <- apply(cox.loc,2,mean)
         n.scale <- length(cox.loc)
         cox.scale <- make.diagonal.matrix(rep(scale[1,1],n.scale))
       }

return(list(
	beta.ipd = beta.ipd,
	beta.ad = beta.ad,
	D = D,
	cox.loc = cox.loc,
	cox.scale = cox.scale,
	surv.loc = loc,
	surv.scale = scale)
        )
}

###CREATING DATA OBJECTS

coxmeta.data.objects <- function(

                     ipd.formula,
                     meta.formula,
                     random.formula,
                     ipd.data,
                     meta.data,
                     ipd.groups,
                     meta.groups,
                     sigma2,
                     study.group.interaction)
{

#MAKING DATA OBJECTS

	cox.data.object <- list(
	   f=ipd.formula,
	   Z=frailty.model.matrix(random.formula,ipd.data),
	   data=ipd.data,
	   n.groups=ipd.groups
	   )

#MAKE SURV DATA OBJECT

	surv <- meta.data[,all.vars(meta.formula)[1]]

        surv.data.object <- list(

	   f=meta.formula,
	   data=meta.data,
	   S=surv,
	   X=model.matrix(terms(meta.formula),meta.data),
	   Z=frailty.model.matrix(random.formula,meta.data),
	   Sigma=surv.covariance.variance(surv,sigma2,study.group.interaction),
	   n.groups=meta.groups
)

return(
       list(
        cox.data.object,
        surv.data.object)
        )
}


coxmeta.info.coef <-  function(
                 beta.index,
                 beta.ipd,
                 beta.ad,
                 cox,surv,
                 B.cox,
                 B.surv,
                 weights)
{

   cox.offset <- cox$Z%*%t(B.cox)%*%weights
   surv.offset <- surv$Z%*%t(B.surv)%*%weights

   info.beta <- coxmeta.score.info.coef(
          beta.index,
	  beta.ipd,
	  beta.ad,
	  cox,surv,
	  cox.offset,surv.offset
	  )

   info.beta <- as.vector(info.beta$info)
  
   b.cox.offset <- lapply(1:nrow(B.cox),function(i){cox$Z%*%B.cox[i,]})
   b.surv.offset <- lapply(1:nrow(B.surv),function(i){surv$Z%*%B.surv[i,]})
   

   scores.beta <- mapply(
            coxmeta.score.info.coef,
	    cox.offset=b.cox.offset,
	    surv.offset=b.surv.offset,
	    MoreArgs=list(
            beta.index=beta.index,
            beta.ipd=beta.ipd,
            beta.ad=beta.ad,
            cox=cox,surv=surv,score.only=TRUE))

    U2.beta <- if(is.matrix(scores.beta)) apply(scores.beta,2,function(x){outer(x,x)}) else scores.beta^2


    U2.beta <- U2.beta%*%weights

###LOUIS INFO FOR EFFECTS

    info.beta <- info.beta-U2.beta
     info.beta <- matrix(info.beta,
                  sqrt(length(info.beta)),sqrt(length(info.beta)))

return(info.beta)
}

coxmeta.score.info.coef <- function(
                        beta.index,
                        beta.ipd,
                        beta.ad,
                        cox,
                        surv,
                        cox.offset,
                        surv.offset,
                        score.only=FALSE)
{

    #KNOWN COVARIANCE-VARIANCE FOR SURVIVAL OUTCOMES
    Sigma <- surv$Sigma

    #PATIENT LEVEL SCORE/INFO/LIKELIHOOD
    NR.cox <- coxph.with.offset(cox$f,cox$data,beta=beta.ipd,cox.offset)        

    #STUDY LEVEL SCORE/INFO/LIKELIHOOD
    NR.surv <- surv.with.offset(surv$f,surv$data,beta=beta.ad,Sigma,surv.offset)

    #COMBINED
    NR <- combined.score.info(
          beta.index=beta.index,
	  U.cox=NR.cox$score,
	  I.cox=NR.cox$info,
	  loglik.cox=NR.cox$loglik,
	  U.surv=NR.surv$score,
	  I.surv=NR.surv$info,
	  loglik.surv=NR.surv$loglik)

    if(!score.only)  return(list(info=NR$info,score=NR$score)) else NR$score
}


#####VARIANCE PARAMETERS

coxmeta.info.vcov <- function(
                  B.cox,
                  B.surv,
                  ipd.groups,
                  meta.groups,
                  D,
                  weights)
{

	info.var <- mixed.vcov.info(D,ipd.groups+meta.groups)
	U2.var <- mixed.vcov.score(B.cox,B.surv,D,ipd.groups,meta.groups)
	U2.var <- U2.var%*%weights
	

	info.var <- ifelse(info.var-U2.var>0,info.var-U2.var,info.var)

return(info.var)
}


mixed.vcov.info <- function(
                D,
                n.groups)
{
     
     f <- function(x){n.groups/(2*x^2)}
     theta <- get.diag(D)
     sapply(theta,f)

}

mixed.vcov.score <- function(
                 B.cox,
                 B.surv,
                 D,
                 cox.groups,
                 surv.groups)
{

###RETURNS THE QUADRATIC OF VCOV SCORES FOR DIAGONAL COMPONENTS

       score <- function(
             b.cox,
             b.surv,
             D,
             cox.groups,
             surv.groups)
          {

              n <- dim(D)[1]

	      b.mat <- 
	         combined.subvector.matrix(b.cox,b.surv,cox.groups,surv.groups)

	      scr <- function(i){
	      	    -1/2*(length(b.mat[,i])/D[i,i]-sum(b.mat[,i]^2)/D[i,i]^2)
	      }

	      sapply(1:n,scr)^2
	   }
   
   b.cox.list <- lapply(1:nrow(B.cox),function(i){B.cox[i,]})
   b.surv.list <- lapply(1:nrow(B.surv),function(i){B.surv[i,]})

   mapply(score,b.cox=b.cox.list,b.surv=b.surv.list,
	MoreArgs=list(D=D,cox.groups=cox.groups,surv.groups=surv.groups))
}


coxmcem.mixed.initialize <- function(
		   coxme.formula,
		   phmm.formula,
		   data,
		   init.coef,
		   init.vcov)
{

#TRY COXME FOR INITIAL PARAMETERS AND FRAILTY PROPOSAL
#IF FAILURE GO TO PHMM

#COXME FORMULA
#PHMM FORMULA IS A LIST WITH FIXED AND RANDOM PARTS

	no.error <- tryCatch(coxme(coxme.formula,data),
			error=function(void){FALSE})

	if(!is.logical(no.error)) no.error = TRUE

   	if(no.error){

          fit <- coxme(coxme.formula,data)
	  coef <- fit$coef$fixed
	  vcov <- fit$coef$random[[1]]
	  mean <- as.vector(fit$frail[[1]])
         }
        else{

	  fit <- phmm(phmm.formula$fixed,phmm.formula$random,data)
	  
	  coef <- fit$coef
	  vcov <- fit$Sigma
	  
	  mean <- as.vector(fit$bhat)
         }

	 if(!missing(init.coef)) coef <- init.coef
	 if(!missing(init.vcov)) vcov <- init.vcov

	 n.frailties <- length(mean)		#SCALE BASED ON EFFECTIVE
	 n.eff <- nrow(data)/n.frailties	#OBS PER FRAILTY
	 Sigma <- diag(rep(1/n.eff,n.frailties))

return(
       list(
	beta=coef,
	D=vcov,
	mean=mean,
	Sigma=Sigma)
       )
}


coxmeta.importance.weights <- function(
                           cox,
                           surv,
                           B.cox,
                           B.surv,
                           D,
                           beta.cox,
                           beta.surv,
                           loc.cox,
                           scale.cox,
                           loc.surv,
                           scale.surv,
                           df)
{

   LL <- coxmeta.frailty.loglik(cox,surv,B.cox,B.surv,D,beta.cox,beta.surv)

   LL.proposal.cox <- proposal.dmt(B.cox,loc.cox,scale.cox,df)
   LL.proposal.surv <- proposal.dmt(B.surv,loc.surv,scale.surv,df)
   LL.proposal <- LL.proposal.cox+LL.proposal.surv


   w <- importance.weights(LL,LL.proposal)
   
return(list(weights=w,sd.loglik=sd(LL),loglik=LL%*%w))
}

######LIKELIHOOD COMPONENTS OF DATA FOR GIVEN FRAILTY

coxmeta.frailty.loglik <- function(
                       cox,
                       surv,
                       B.cox,
                       B.surv,
                       D,
                       beta.cox,
                       beta.surv)
{

   loglik.cox <- apply(B.cox,1,function(b){
      cox.mixed.loglik(cox$f,cox$data,cox$Z,b,beta.cox)
   })

   loglik.surv <- apply(B.surv,1,function(b){
      surv.mixed.loglik(surv$S,surv$X,beta.surv,surv$Z,b,surv$Sigma) 
   })

   loglik.penalty.cox <- penalty.mixed.loglik(t(B.cox),D,cox$n.groups)

   loglik.penalty.surv <- penalty.mixed.loglik(t(B.surv),D,surv$n.groups)

   loglik = loglik.cox+loglik.surv+loglik.penalty.cox+loglik.penalty.surv

return(loglik)
}


cox.mixed.loglik <- function(
                 f,
                 data,
                 Z,
                 b,
                 beta)
{

   offset <- Z%*%b
   f <- paste(c(f,"offset(offset)"),collapse="+")
   f <- as.formula(f)
     				#LOGLIK DATA CONDITIONAL ON CLUSTER EFFECTS
   L1 <- coxph(f,data,init=beta,iter=0)$loglik[2]
   return(L1)
}


surv.mixed.loglik <- function(
                  S,
                  X,
                  beta,
                  Z,
                  b,
                  Sigma)
{

   if(is.matrix(X)) u = X%*%beta else u=X*beta

   R = log(-log(S))-u-Z%*%b

   -1/2*t(R)%*%solve(Sigma)%*%R-1/2*log(det(Sigma))   
}

penalty.mixed.loglik <- function(
                     B,
                     D,
                     n.groups)
{

   #B HAS EACH COLUMN AS A SAMPLE OF FRAILTIES
   #FRAILTIES ARE ORDERED BY COVARIATE, I.E. INTERCEPT, TRT, ETC.

   sigma <- rep(get.diag(D),each=n.groups)

   apply(-1/2*B^2/sigma,2,sum)-1/2*sum(log(sigma))
}

#######LIKELIHOOD FOR MULTIVARIATE T PROPOSAL

proposal.dmt <- function(
             B,
             loc,
             scale,
             df,
             log=TRUE)
{

   apply(B,1,function(b){
   	       dmt(b,mean=loc,S=scale,df=df,log=log)
	       })
}


importance.weights <- function(
                   loglik,
                   loglik.proposal)
{

   w <- loglik-loglik.proposal-max(loglik-loglik.proposal)

   exp(w)/sum(exp(w))
}




















