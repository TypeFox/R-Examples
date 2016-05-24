mbscore <- function(x,lambda,subset=NULL,tol=1e-8,...){
	beta <- scores <- scorevar.mod <- scorevar.sand <- p.model <- p.sand <- matrix(NA,ncol(x),ncol(x))
	
	if(is.null(subset)) subset <- matrix(TRUE,ncol(x),ncol(x))

	for(i in 1:ncol(x)){
		mod <- lassoscore(x[,i],x[,-i],lambda=lambda,subset=which(subset[-i,i]),tol=tol,...)
		scores[-i,i] <- mod$scores
		p.model[-i,i] <- mod$p.model
		p.sand[-i,i] <- mod$p.sand
		scorevar.mod[-i,i] <- mod$scorevar.model
		scorevar.sand[-i,i] <- mod$scorevar.sand
		beta[-i,i] <- mod$fit$beta	
	}
	re <- list("scores"=scores,
				"p.model" = p.model,
				"p.sand" = p.sand,
				"scorevar.mod" = scorevar.mod,
				"scorevar.sand" = scorevar.sand,
				"beta" = beta)
	re$lambda <- lambda
	class(re) <- "mbscore"			
	return(re)
}



qqpval <- function(p, cone=TRUE, log=TRUE, add=FALSE,col=1,pch=1,...){
	if(length(col)==1) col <- rep(col,length(p))
	if(length(pch)==1) pch <- rep(pch,length(p))
	stopifnot(length(col) == length(p) & length(pch) == length(p))
	col <- col[0< p & p < 1]
	pch <- pch[0 < p & p <1]
	p <- p[ 0< p & p < 1]
	
	oo <- order(p,decreasing=log)
	p <- p[oo]
	col <- col[oo]
	pch <- pch[oo]
	
	if(log) p <- -log10(p)
	N <- length(p)
	if(log){
		Ep <- sort(-log10(ppoints(N)))
		cone.lower <- -log10(qbeta(0.025, 1:N, N:1))
		cone.upper <- -log10(qbeta(0.975, 1:N, N:1))
		cone.x <- sort(-log10(ppoints(N)))
	} else {
		Ep <- sort(ppoints(N))
		cone.lower <- sort(qbeta(0.025, 1:N, N:1))
		cone.upper <- sort(qbeta(0.975, 1:N, N:1))
		cone.x <- rev(sort(ppoints(N)))
	}
	if(add){
		points(y=p, x=Ep,col=col,pch=pch,...)
	} else{
		if(log){
			plot(y=p,x=Ep,type="n", xlab = "Expected p-value, -log10", 
           ylab = "Observed p-value, -log10",ylim=range(c(p,cone.upper,cone.lower)),...)
		} else {
		  plot(y=p,x=Ep,type="n", xlab = "Expected p-value", 
           ylab = "Observed p-value",...)
		}
		abline(0,1)
		polygon(y=c(rev(cone.lower),cone.upper),x=c(cone.x,rev(cone.x)), 
            col=1,lty=2,density=0)
		points(y=p, x=Ep,col=col,pch=pch,...)
	}
}

glassoscore <- function(x,lambda,subset=NULL,penalize.diagonal=FALSE,tol=1e-8){
  if(!is.matrix(x)) stop("x must be a matrix")
  if(lambda[1] < 0 | length(lambda) > 1) stop("lambda must be a non-negative scalor")
  n <- nrow(x)
  d <- ncol(x)
  x <- scale(x)*sqrt((n-1)/n)
  s <- var(x)*((n-1)/n)

  if(is.null(subset)) subset <- matrix(TRUE,d,d)
  stopifnot(all(dim(subset)==dim(s)) & is.logical(subset))
  
  mod0 <- glasso(s,lambda,thr=tol,penalize.diagonal=penalize.diagonal)
  
  scores <- matrix(0,d,d)
  scorevars.sand <- scorevars.model <- matrix(0,d,d)
  
  lambdamat <- matrix(lambda,d,d)
  diag(lambdamat) <- 0
  
  wh <- which(mod0$wi !=0)
  adj <- cbind(floor((wh-0.1)/d)+1,(wh-1)%%d+1)
  wh <- wh[adj[,1] <= adj[,2]]
  adj <- adj[adj[,1] <= adj[,2],]
  
  U0 <- with(mod0,{
 		w[adj[,1],adj[,1]]*w[adj[,2],adj[,2]] + w[adj[,2],adj[,1]]*w[adj[,1],adj[,2]]
      })
  Ui0 <- solve(U0)   
  V0 <- with(mod0, var(x[,adj[,1]]*x[,adj[,2]]))*((n-1)/n)
  sand0 <- V0%*%Ui0%*%V0
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      if(subset[i,j] | subset[j,i]){
        if(mod0$wi[i,j] ==0){ 
          mod <- mod0
        } else {
          mod <- glasso(s,rho=lambda,zero=c(i,j),thr=tol,start="warm",w.init=mod0$w,
                        wi.init=mod0$wi,penalize.diagonal=penalize.diagonal)
        }
        
        scores[i,j] <- with(mod,sqrt(n)*(w[i,j]-s[i,j]))
        
        wh <- c((i-1)*d+j,which(mod$wi !=0))
        adj <- cbind(floor((wh-0.1)/d)+1,(wh-1)%%d+1)
        wh <- wh[adj[,1] <= adj[,2]]
        adj <- adj[adj[,1] <= adj[,2],]
        
        if(mod0$wi[i,j] ==0){
          Ua <- with(mod,{
            w[adj[1,1],adj[,1]]*w[adj[1,2],adj[,2]] + w[adj[1,2],adj[,1]]*w[adj[1,1],adj[,2]]
          })
          V <- with(mod, var(x[,adj[,1]]*x[,adj[,2]]))*((n-1)/n)

          scorevars.sand[i,j] <- (V[1,1] + t(Ua[-1])%*%sand0%*%Ua[-1]  - 2*t(Ua[-1] )%*%(Ui0)%*%V[-1,1])
          scorevars.model[i,j] <- Ua[1] - t(Ua[-1])%*%Ui0%*%Ua[-1]
        } else {
          U <- with(mod,{
            w[adj[,1],adj[,1]]*w[adj[,2],adj[,2]] + w[adj[,2],adj[,1]]*w[adj[,1],adj[,2]]
          }) 
          V <- with(mod, var(x[,adj[,1]]*x[,adj[,2]]))*((n-1)/n)
          
          Ui <- solve(U[-1,-1])
          scorevars.sand[i,j] <- V[1,1] + t(U[-1,1])%*%Ui%*%(V[-1,-1]%*%Ui%*%U[-1,1] - 2*V[-1,1])
          scorevars.model[i,j] <- U[1,1] - t(U[-1,1])%*%Ui%*%U[-1,1]
        }    
      }
    }
  }
  
  scores <- as.matrix(forceSymmetric(scores))
  scorevars.sand <- as.matrix(forceSymmetric(scorevars.sand))
  scorevars.model <- as.matrix(forceSymmetric(scorevars.model))
  
  mod0$p.sand <- pchisq(scores^2/scorevars.sand,1,lower.tail=FALSE)
  mod0$p.model <- pchisq(scores^2/scorevars.model,1,lower.tail=FALSE)
  
  mod0$scores <- scores
  mod0$scorevar.sand <- scorevars.sand
  mod0$scorevar.model <- scorevars.model  
  mod0$lambda <- lambda
  
  class(mod0) <- "glassoscore"
  return(mod0)
}

print.glassoscore <- function(x,...){
  cat("An object of class `glassoscore' \n")
  cat("lambda = ", x$lambda[1])
}

print.mbscore <- function(x,...){
  cat("An object of class `mbscore' \n")
  cat("lambda = ", x$lambda[1])
}

plot.glassoscore <- function(x,cone=TRUE, log=TRUE, add=FALSE,...){
  qqpval(x$p.sand[upper.tri(x$p.sand)],cone,log,add, main = "QQ-plot of p-values",...)
  qqpval(x$p.model[upper.tri(x$p.sand)],add=TRUE,col=2,...)
  legend("topleft",legend=c("sandwich","model-based"),col=c(1,2),pch=1)
}

plot.mbscore <- function(x,cone=TRUE, log=TRUE, add=FALSE,...){
  qqpval(x$p.sand[upper.tri(x$p.sand)],cone,log,add, main = "QQ-plot of p-values",...)
  qqpval(x$p.model[upper.tri(x$p.sand)],add=TRUE,col=2,...)
  legend("topleft",legend=c("sandwich","model-based"),col=c(1,2),pch=1)
}



lassoscore <- function(y,X, lambda=0, family=c("gaussian","binomial","poisson"), tol = .Machine$double.eps, maxit=1000, 
                       resvar = NULL, verbose=FALSE, subset = NULL){
  family = match.arg(family)
  family.obj <- get(family,mode="function",envir=parent.frame())()
  
  if(!(family %in% c("gaussian","binomial","poisson"))) stop("family not supported!")                   
  X <- scale(X)*sqrt(nrow(X)/(nrow(X)-1))
  if(family=="gaussian") y <- scale(y,scale=FALSE)
  
  ##initial fit
  out0 <- glmnet(y=y,x=X, family=family, lambda=lambda, thresh = tol, maxit=maxit,standardize=FALSE)
  out0$r <- as.vector(family.obj$linkinv(as.vector(out0$a0 + X%*%out0$beta))-y)
  out0$v <- family.obj$var(family.obj$linkinv(as.vector(out0$a0 + X%*%out0$beta)))
  
  
  if(family =="gaussian" & !is.numeric(resvar)){
    resvar <- sum(out0$r^2)/(length(y)-sum(out0$beta !=0))
  } else if(family !="gaussian"){
    resvar <- 1
  }
  
  out0$n <- nrow(X)
  out0$beta <- as.vector(out0$beta)
  
  wh <- as.vector(out0$beta != 0)
  if(is.null(subset)){
    subset <- 1:ncol(X)
  }
  if(verbose){
    cat("\nProgress:\n")
    pb <- txtProgressBar(min = 0, max = length(subset), style = 3)
    pb.i <- 0
  }
  scores <- scorevar.sand.cons <- scorevar.model.cons <- scorevar.sand <- scorevar.model <- numeric(ncol(X))
  scores[-subset] <- NA
  
  for(i in subset){
  	if(out0$beta[i] != 0){
    	out <- glmnet(y=y,x=X[,-i], family=family, lambda=lambda, thresh = tol, maxit=maxit,standardize=FALSE)
    	out$r <- as.vector(family.obj$linkinv(out$a0 +as.vector(X[,-i]%*%out$beta))-y)
	    out$v <- family.obj$var(family.obj$linkinv(out$a0 +as.vector(X[,-i]%*%out$beta)))  		
    	out$n <- nrow(X)
    	out$beta <- as.vector(out$beta)
    	Xs <- X[,-i,drop=FALSE][,as.vector(out$beta !=0),drop=FALSE]
    } else {
    	out <- out0
    	Xs <- X[,as.vector(out$beta !=0),drop=FALSE]
    }
    
    xx <- X[,i,drop=FALSE]
    
    scores[i] <- sum(out$r*X[,i])/sqrt(out$n)
       
    scorevar.model.cons[i] <- with(out, resvar*(crossprod(xx*sqrt(v))/n))
    scorevar.sand.cons[i] <- with(out, var(xx*r)*(n-1)/n)

    if(ncol(Xs) == 0){
      scorevar.sand[i] <- scorevar.sand.cons[i]
      scorevar.model[i] <- scorevar.model.cons[i]	
    } else if(nrow(Xs) > ncol(Xs) & ncol(Xs) >0){ 
      Ui <- with(out, solve(crossprod(Xs*sqrt(v))/n))
      V <- with(out, var(r*Xs)*(n-1))
      Va <- with(out,cov(r*Xs,r*xx)*(n-1))
      Ua <-  with(out, crossprod(Xs,v*xx)/n)
      va <-  with(out, var(xx*r)*(n-1))
      
      scorevar.sand[i] <- (va + t(Ua)%*%(Ui)%*%(V%*%Ui%*%Ua - 2*Va))/with(out,n-sum(beta !=0))
      scorevar.model[i] <- with(out, resvar*(crossprod(xx*sqrt(v))/n - t(Ua)%*%Ui%*%Ua)	)  
    }
    if(verbose){
      pb.i <- pb.i+1
      setTxtProgressBar(pb, pb.i)
    }
  }
  re <- list(
    "fit"=out0,
    "scores" = scores,
    "scorevar.model.cons" = scorevar.model.cons,
    "scorevar.sand.cons" = scorevar.model.cons,
    "scorevar.model" = scorevar.model,
    "scorevar.sand" = scorevar.sand,	
    "p.model.cons" = pchisq(scores^2/scorevar.model.cons,df=1,lower.tail=FALSE),
    "p.sand.cons" = pchisq(scores^2/scorevar.sand.cons,df=1,lower.tail=FALSE),
    "p.model" = pchisq(scores^2/scorevar.model,df=1,lower.tail=FALSE),
    "p.sand" = pchisq(scores^2/scorevar.sand,df=1,lower.tail=FALSE),
    "lambda" = lambda)
  
  class(re) <- "lassoscore"
  if(verbose) close(pb)
  return(re)
}


print.lassoscore <- function(x,...){
	cat("An object of class `lassoscore'\nbased on n =",x$fit$n, " observations on d =",x$fit$dim[1], " features, using lambda = ",x$lambda,fill=TRUE)
}

summary.lassoscore <- function(object,...){
  cat("An object of class `lassoscore'\nbased on n =",object$fit$n, " observations on d =",object$fit$dim[1], " features, with",sum(object$fit$beta !=0),"non-zero coefficients in regression of `y' on `X'",fill=TRUE)
  cat("\nModel-based p-values:\n")
      print(summary(object$p.model))
  cat("\nSandwich p-values:\n")
    print(summary(object$p.sand))
   }

plot.lassoscore <- function(x, type=c("all","model","sand","model.cons","sand.cons"),cone=TRUE, log=TRUE, add=FALSE,...){
	type <- match.arg(type,choices = c("all","model","sand","model.cons","sand.cons"))
	switch(type,
		sand = qqpval(x$p.sand,cone,log,add, main = "QQ-plot of p-values",...),
		model = qqpval(x$p.model,cone,log,add,main = "QQ-plot of p-values",...),
		sand.cons = qqpval(x$p.sand.cons,cone,log,add,main = "QQ-plot of p-values",...),
		model.cons = qqpval(x$p.model.cons,cone,log,add,main = "QQ-plot of p-values",...),
		all = { 
			qqpval(x$p.model,cone,log,main = "QQ-plot of p-values",pch=16,col="#1B9E77",...)
			qqpval(x$p.sand,cone,log,add=TRUE,col="#D95F02",pch=18,...)
			qqpval(x$p.model.cons,cone,log,add=TRUE,col="#7570B3",...)
			qqpval(x$p.sand.cons,cone,log,add=TRUE,col="#E7298A",pch=5,...)
			legend("topleft",col=c("#1B9E77", "#D95F02" ,"#7570B3", "#E7298A"),legend=c("model-based", "sandwich", "conservative, model based","conservative, sandwich"),pch=c(16,18,1,5))}
		)		
}

