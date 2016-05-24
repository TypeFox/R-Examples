coxmcem <- 

function(

	fixed,
	random,
	data,
	df,
	n.groups,
	max.iter=10,
	min.sample=100,
	mc.step=2,
	est.delta=5/100,
	echo = TRUE,
	init.coef,
	init.vcov
){

   ### MISSING REMOVAL
   vars <- c(all.vars(fixed),all.vars(random))
   missing.cases <- sapply(vars,function(x){any(is.na(data[,x]))})

   if(any(missing.cases)){
     i <- which(missing.cases)
     missing.index <- sapply(vars[i],function(x){which(is.na(data[,x]))})
     missing.index <- unique(as.vector(unlist(missing.index)))
     data <- data[-missing.index,]
     print(paste(length(missing.index),"Missing cases excluded"))
   }

   p.beta <- length(all.vars(fixed))-2
   Z <- frailty.model.matrix(random,data)

#PARAMETER INITIALIZATION

   initialized <- 
      coxmcem.initialize.coxmcem(fixed,random,data,p.beta,init.coef,init.vcov)

   beta <- initialized$beta
   D <- initialized$D

###LOGLIKELIHOOD FOR INITIALIZED PARAMETERS
   LL1 <- coxme.loglik.cox.part(fixed,data,Z,initialized$mean,beta)
      sigma <- rep(get.diag(D),each=n.groups)

   LL2 <- sum(-1/2*initialized$mean^2/sigma)-1/2*sum(log(sigma))

   loglik.max <- LL1+LL2

   proposal.rmt <- function(n){

   	       rmt(n,mean=initialized$mean,S=initialized$Sigma,df=df)
	       }
	      
   proposal.dmt <- function(x){
   	       dmt(x,mean=initialized$mean,S=initialized$Sigma,df=df,log=TRUE)
	       }

 
#CONVERGENCE/SAMPLE SIZE INITIALIZATION
   N = min.sample
   n.iter <- 0
   converged = FALSE

#COXMCEM OUTPUT
   coxmcem.object <- 
      list(
      max.weight=c(),
      mc.samples=c(),
      est.converge=c(),
      loglik=c(),
      sd.loglik=c()
      )
  
    
while(!converged&n.iter<max.iter){

   n.iter <- n.iter+1
#SAMPLING RANDOM EFFECTS

b.tilde <- proposal.rmt(N)
   

#LIKELIHOOD FUNCTIONS: DEPEND ON CURRENT BETA AND D
loglik.cox <- apply(b.tilde,1,function(b){
	L1 <- coxme.loglik.cox.part(fixed,data,Z,b,beta)
	L1
})

loglik.penalty <- penalty.loglik.coxmcem(t(b.tilde),D,n.groups)

LL <- loglik.cox+loglik.penalty

#UPDATE LOCATION OF PROPOSAL WITH MOST LIKELY SET OF FRAILTIES
#IF MORE LIKELY THAN IN PREVIOUS ITERATION

if(loglik.max<max(LL)){

	 loglik.max=max(LL)
	 i.loglik=which(LL==loglik.max)
	 initialized$mean <- b.tilde[i.loglik,]
}


loglik.proposal <- apply(b.tilde,1,proposal.dmt)

#IMPORTANCE SAMPLING WEIGHTS
   
   w <- LL-loglik.proposal-max(LL-loglik.proposal)
   w <- exp(w)/sum(exp(w))

#E STEP
   E <- importance.Estep(t(b.tilde),Z,w,n.groups)

#M STEP
   M <- importance.Mstep(fixed,E,data)

###RELATIVE CHANGE

zeta.now = c(M$cox.fit$coef,get.diag(M$D))
zeta.last = c(beta,get.diag(D))

est.converge <- relative.criteria(zeta.now,zeta.last)

###UPDATE COX OBJECT

coxmcem.object$loglik <- append(coxmcem.object$loglik,sum(LL*w))
coxmcem.object$sd.loglik <- append(coxmcem.object$sd.loglik,sd(LL))
coxmcem.object$max.weight <- append(coxmcem.object$max.weight,max(w))

tries <- length(coxmcem.object$max.weight)
m <- length(coxmcem.object$est.converge)

if(tries>3) cv.last <- sd(coxmcem.object$est.converge[(m-2):m])/mean(coxmcem.object$est.converge[(m-2):m])

coxmcem.object$est.converge <- append(coxmcem.object$est.converge,est.converge)

if(tries>3) cv.now <- sd(coxmcem.object$est.converge[(m-1):(m+1)])/mean(coxmcem.object$est.converge[(m-1):(m+1)])


coxmcem.object$mc.samples <- append(coxmcem.object$mc.samples,N)

###CHECK IF CONVERGENCE IS MET; UPDATE SAMPLE SIZE IF NOT


if(tries>2){

   est.convergence <- max(coxmcem.object$est.converge[(m-1):(m+1)])
  
      if(est.convergence<est.delta){
         converged = TRUE
	}

if(tries>3){

      if(cv.now>cv.last){
	N = ceiling(N+N/mc.step)
       }
 }	
}

#UPDATE ESTIMATES; PRINT CONVERGENCE INFO
   if(tries>3&echo) print(paste(c("E-step n","Max relative change %"),
   round(c(coxmcem.object$mc.samples[length(coxmcem.object$mc.samples)],
	est.converge*100),5)))

   beta <- coef(M$cox.fit)
   D <- M$D
}

#FINAL ESTIMATES

   b0 <- matrix(0,length(b.tilde[1,]),1)
   coxmcem.object$loglik.fixed <- coxme.loglik.cox.part(fixed,data,Z,b0,beta)
   coxmcem.object$iterations <- n.iter
   coxmcem.object$coef <- beta
   coxmcem.object$vcov <- D
   coxmcem.object$cluster <- E$mean
   
   var <- coxme.variance(b.tilde,w,fixed,data,beta,D,Z,n.groups)
   coxmcem.object$var$coef <- solve(var$info.beta)
   coxmcem.object$var$vcov <- 1/var$info.var

return(coxmcem.object)
}

####### DEPENDENT FUNCTIONS


process.formula <- 

function(f){

   #RETURNS FORMULA PART GIVEN A COXME-PARANTHETICAL FORMULA		
   str <- as.character(f)[2]
   str <- strsplit(strsplit(str,"\\|")[[1]][1],"\\(")[[1]][2]
   as.formula(paste(c("~",str),collapse=""))
}


coxmcem.formula <- 

function(fixed,random,type=c("coxme","phmm"))
{
	if(type=="coxme"){

	f <- formula(paste(c(fixed,as.character(random)[2]),collapse="+"))
	return(f)

	}
	else{

	fixed.part <- formula(paste(c(fixed,"cluster(cluster)"),collapse="+"))
	random.part <- process.formula(random)

	return(list(fixed=fixed.part,random=random.part))
	}
}

relative.criteria <- function(now,last){

   denom = abs(last)

   if(any(abs(last)<.001)) denom = denom+.001

   Delta <- abs(now-last)/denom

   max(Delta)
}

coxme.loglik.cox.part <- function(f,data,Z,b,beta){

   offset <- Z%*%b
   f <- paste(c(f,"offset(offset)"),collapse="+")
   f <- as.formula(f)
     				#LOGLIK DATA CONDITIONAL ON CLUSTER EFFECTS
   L1 <- coxph(f,data,init=beta,iter=0)$loglik[2]
   return(L1)
}

###LOGLIKELIHOOD OF NORMAL CONTRIBUTION

penalty.loglik.coxmcem <- function(B,D,n.groups){

   #B HAS EACH COLUMN AS A SAMPLE OF FRAILTIES
   #FRAILTIES ARE ORDERED BY COVARIATE, I.E. INTERCEPT, TRT, ETC.

   sigma <- rep(get.diag(D),each=n.groups)

   apply(-1/2*B^2/sigma,2,sum)-1/2*sum(log(sigma))
}


coxmcem.initialize.coxmcem <- function(fixed,random,data,p.beta,init.coef,init.vcov){

#TRY COXME FOR INITIAL PARAMETERS AND FRAILTY PROPOSAL
#IF FAILURE GO TO PHMM

	f <- coxmcem.formula(fixed,random,"coxme")

	no.error <- tryCatch(coxme(f,data),
			error=function(void){FALSE})

	if(!is.logical(no.error)) no.error = TRUE

   	if(no.error){

          fit <- coxme(f,data)
	  coef <- fit$coef$fixed
	  vcov <- if(is.matrix(fit$coef$random[[1]])) diag(fit$coef$random[[1]]) else fit$coef$random[[1]]
	  vcov <- make.diagonal.matrix(vcov)

	  mean <- as.vector(fit$frail[[1]])
         }

        else{

   	  f <- coxmcem.formula(fixed,random,"phmm")
	  fit <- phmm(f$fixed,f$random,data)
	  
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
			 Sigma=Sigma
		 )
	 )
}


coxme.variance <- function(B,weights,formula,data,beta,D,Z,n.groups){

   p.beta <- length(beta)
   p.D <- ncol(cluster.subvector.matrix(B[1,],n.groups))

###EFFECTS

info.beta <- coxph.with.offset.coxmcem(formula,data,beta,Z%*%t(B)%*%weights)$info

scores.beta <- apply(B,1,function(b){
      coxph.with.offset.coxmcem(formula,data,beta,Z%*%b)$score
})

U2.beta <- if(is.matrix(scores.beta)) apply(scores.beta,2,function(x){outer(x,x)}) else scores.beta^2


U2.beta <- U2.beta%*%weights

###LOUIS INFO FOR EFFECTS

info.beta <- info.beta-U2.beta
info.beta <- matrix(info.beta,p.beta,p.beta)

###VARIANCE

info.var <- my.vcov.info.coxmcem(D,n.groups)
U2.var <- my.vcov.score.coxmcem(B,D,n.groups)

U2.var <- U2.var%*%weights

info.var <- ifelse(info.var-U2.var>0,info.var-U2.var,info.var)

return(list(info.beta=info.beta,info.var=info.var))
}


coxph.with.offset.coxmcem <- function(formula,data,beta,offset){

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


   return(list(score=U,info=I,loglik=fit$loglik[2]))
}


my.vcov.info.coxmcem <- function(D,n.groups){
     
     f <- function(x){n.groups/(2*x^2)}
     theta <- get.diag(D)
     sapply(theta,f)

}

my.vcov.score.coxmcem <- function(B,D,n.groups){

###RETURNS THE QUADRATIC OF VCOV SCORES FOR DIAGONAL COMPONENTS

	   score <- function(b,D,n.groups){

      	      n <- dim(D)[1]
	      b.mat <- cluster.subvector.matrix(b,n.groups)

	      scr <- function(i){
	      	    -1/2*(length(b.mat[,i])/D[i,i]-sum(b.mat[,i]^2)/D[i,i]^2)
	      }

	      sapply(1:n,scr)^2
	   }
   
   b.list <- lapply(1:nrow(B),function(i){B[i,]})

   mapply(score,b=b.list,MoreArgs=list(D=D,n.groups=n.groups))
}

importance.Estep <- function(b,Z,weights,n.groups){

  U <- apply(b,2,function(x){Z%*%x})	#OFFSET FOR EACH SAMPLED RANDOM EFFECT
  offset <- U%*%weights  
 
  means <- b%*%weights
  means <- cluster.subvector.matrix(means,n.groups)	#MEANS BY GROUP
  
  E2 <- apply(b,2,function(x){
  	M <- cluster.subvector.matrix(x,n.groups)
	V <- apply(M,1,function(x){outer(x,x)})
	return(V)
	})
 
  E2 <- E2%*%weights
  E2 <- cluster.subvector.matrix(E2,nrow(E2)/n.groups)
 
return(list(offset=offset,mean=means,second.moment=E2))
}


importance.Mstep <- function(f,importance.object,data){

	offset <- importance.object$offset
	f.offset <- as.formula(paste(c(f,"offset(offset)"),collapse="+"))

	#FIXED EFFECTS
	fit <- coxph(f.offset,data=data)

	#VARIANCE-COVARIANCE
	p <- ncol(importance.object$mean)
	D <- apply(importance.object$second.moment,1,mean)
	D <- matrix(D,p,p)

return(list(cox.fit=fit,D=D))
}


make.square.matrix <- function(x){
   n = length(x)
   p = sqrt(n)
   matrix(x,p,p)		   
}

make.diagonal.matrix <- function(x){
   n = length(x)
   D = matrix(0,n,n)
   diag(D) = x
   D
}

dmt <- 

function (x, mean = rep(0, d), S, df = Inf, log = FALSE) 
{
    if (df == Inf) 
        return(dmnorm(x, mean, S, log = log))
    d <- if (is.matrix(S)) 
        ncol(S)
    else 1
    if (d > 1 & is.vector(x)) 
        x <- matrix(x, 1, d)
    n <- if (d == 1) 
        length(x)
    else nrow(x)
    X <- t(matrix(x, nrow = n, ncol = d)) - mean
    Q <- apply((solve(S) %*% X) * X, 2, sum)
    logDet <- sum(logb(abs(diag(qr(S)$qr))))
    logPDF <- (lgamma((df + d)/2) - 0.5 * (d * logb(pi * df) + 
        logDet) - lgamma(df/2) - 0.5 * (df + d) * logb(1 + Q/df))
    if (log) 
        logPDF
    else exp(logPDF)
}

rmt <- function (n = 1, mean = rep(0, d), S, df = Inf) 
{
    d <- if (is.matrix(S)) 
        ncol(S)
    else 1
    if (df == Inf) 
        x <- 1
    else x <- rchisq(n, df)/df
    z <- rmnorm(n, rep(0, d), S)
    y <- t(mean + t(z/sqrt(x)))
    return(y)
}

dmnorm <- 
function (x, mean = rep(0, d), varcov, log = FALSE) 
{
    d <- if (is.matrix(varcov)) 
        ncol(varcov)
    else 1
    if (d > 1 & is.vector(x)) 
        x <- matrix(x, 1, d)
    n <- if (d == 1) 
        length(x)
    else nrow(x)
    X <- t(matrix(x, nrow = n, ncol = d)) - mean
    Q <- apply((solve(varcov) %*% X) * X, 2, sum)
    logDet <- sum(logb(abs(diag(qr(varcov)[[1]]))))
    logPDF <- as.vector(Q + d * logb(2 * pi) + logDet)/(-2)
    if (log) 
        logPDF
    else exp(logPDF)
}

rmnorm <- 
function (n = 1, mean = rep(0, d), varcov) 
{
    d <- if (is.matrix(varcov)) 
        ncol(varcov)
    else 1
    z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
    y <- t(mean + t(z))
    return(y)
}
