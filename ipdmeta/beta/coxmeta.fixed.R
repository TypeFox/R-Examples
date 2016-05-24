coxmeta.fixed <- 

function(
ipd.formula,
meta.formula,
ipd.data,
meta.data,
sigma2,
study.group.interaction,
error=1/1000
)

{

    beta.index <- beta.indices(ipd.formula,meta.formula)
    beta.ad = runif(max(beta.index$ad))
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

return(list(
	coef=beta,
	var=solve(NR$info),
	loglik=loglik,
	monitor=monitor))
}

####DEPENDENT FUNCTIONS

combined.score.info <- 

function(beta.index,U.cox,I.cox,loglik.cox,U.surv,I.surv,loglik.surv){

   p <- max(beta.index$ad)
   U <- rep(0,p)
   U[beta.index$ipd] <- U.cox
   U[beta.index$ad] <- U[beta.index$ad]+U.surv

   #COMBINING INFO
   I <- matrix(0,p,p)
   I[beta.index$ipd,beta.index$ipd] <- I.cox
   I[beta.index$ad,beta.index$ad] <- I[beta.index$ad,beta.index$ad]+I.surv

   loglik <- loglik.cox+loglik.surv

   return(list(score=U,info=I,loglik=loglik))
}


coxph.with.offset <- function(formula,data,beta,offset){

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

   return(list(score=U,info=I,loglik=fit$loglik[2]))
}


surv.with.offset <- function(formula,data,beta,Sigma,offset){

   X <- model.matrix(terms(formula),data)

   S <- data[,all.vars(terms(formula))[1]]
   
   R <- log(-log(S))-X%*%beta
   
   if(!missing(offset)) R <- R-offset

   U <- as.vector(t(X)%*%solve(Sigma)%*%R)
   I <- as.vector(t(X)%*%solve(Sigma)%*%X)
   I <- matrix(I,length(beta),length(beta))

   loglik <- -1/2*t(R)%*%solve(Sigma)%*%R

   return(list(score=U,info=I,loglik=loglik))
}

