
## low.dim is logic, indicating whether this is a low-dimensional problem (p<<n). The function will pick the better implementation accordingly.
library(Matrix)

rnc.linear <- function(X=NULL,Y,A,lambda,gamma=0,cv=NULL,cv.seed=999,low.dim=NULL){
    Adj <- A
	if(!is.null(X)){
	    n <- nrow(X)
	    p <- ncol(X)
	    if(is.null(low.dim)){
		    low.dim <- TRUE
		    if(p>(n/5)) low.dim <- FALSE
	    }
	    Y <- matrix(Y,ncol=1)
	    D <- diag(rowSums(Adj))
	    L <- D - Adj + diag(rep(gamma,n))
	    W <- diag(rep(1,n))
	    if(low.dim){
	    	    tmp.result <- rnc_solver_naive(X,L,Y,lambda,W)
	    }else{
	        tmp.result <- rnc_solver_X(X,L,Y,lambda,W)
	    }
	    alpha <- tmp.result[1:n]
	    beta <- tmp.result[(n+1):(n+p)]
	    cv.MSE <- 0
            cv.sd <- NULL
	    if(!is.null(cv)){
	    	    set.seed(cv.seed)
	        K <- cv
	        eligible.nodes <- 1:n
	        eligible.n <- length(eligible.nodes)
	        if(cv==-1) K <- eligible.n
	        set.seed(500)
	        cv.order <- sample(eligible.n,size=eligible.n)
	        cv.index <- 1:eligible.n
	        cv.index[cv.order] <- 1:K
                cv.MSE.seq <- rep(0,K)
	        for(k in 1:K){
	            valid.index <- eligible.nodes[which(cv.index==k)]
	            current.index <- (1:n)[-valid.index]
	            s.A <- Adj[current.index,current.index]
	            cv.lm.net <- rnc.linear(X=matrix(X[current.index,],ncol=ncol(X)),Y=Y[current.index],A=s.A,lambda=lambda,gamma=gamma,cv=NULL)

	            L11 <- Matrix(L[valid.index,valid.index],sparse=TRUE)
	            L12 <- -1*L[valid.index,current.index]
	            valid.alpha <- solve(L11,L12,sparse=TRUE)%*%cv.lm.net$alpha
	            valid.y <- valid.alpha + X[valid.index,]%*%matrix(cv.lm.net$beta,ncol=1)
	            cv.MSE.seq[k] <- mean((valid.y-Y[valid.index])^2)
	            #print(mean((valid.y-Y[valid.index])^2))
	            #print(cv.MSE)
	        }
	        cv.MSE <- mean(cv.MSE.seq)
                cv.sd <- sd(cv.MSE.seq)/sqrt(K)
	    }
	    return(list(alpha=alpha,beta=beta,cv.MSE=cv.MSE,cv.sd=cv.sd))
    }else{
	    n <- length(Y)
	    Y <- matrix(Y,ncol=1)
	    D <- diag(rowSums(Adj))
	    L <- D - Adj + diag(rep(gamma,n))
	    W <- diag(rep(1,n))
	    tmp.result <- rnc_solver_noX(L,Y,lambda,W)
	    alpha <- tmp.result[1:n]
	    beta <- NULL
	    cv.MSE <- 0
            cv.sd <- NULL

            if(!is.null(cv)){
	        K <- cv
	        eligible.nodes <- 1:n
	        eligible.n <- length(eligible.nodes)
	        if(cv==-1) K <- eligible.n
	        set.seed(500)
	        cv.order <- sample(eligible.n,size=eligible.n)
	        cv.index <- 1:eligible.n
	        cv.index[cv.order] <- 1:K
                cv.MSE.seq <- rep(0,K)
	        for(k in 1:K){
	            valid.index <- eligible.nodes[which(cv.index==k)]
	            current.index <- (1:n)[-valid.index]
	            s.A <- Adj[current.index,current.index]
	            cv.lm.net <- rnc.linear(X=NULL,Y=Y[current.index],A=s.A,lambda=lambda,gamma=gamma,cv=NULL)

	            L11 <- Matrix(L[valid.index,valid.index],sparse=TRUE)
	            L12 <- -1*L[valid.index,current.index]
	            valid.alpha <- solve(L11,L12,sparse=TRUE)%*%cv.lm.net$alpha
	            valid.y <- valid.alpha
	            cv.MSE.seq[k] <- mean((valid.y-Y[valid.index])^2)
	            #print(mean((valid.y-Y[valid.index])^2))
	            #print(cv.MSE)
	        }
	        cv.MSE <- mean(cv.MSE.seq)
                cv.sd <- sd(cv.MSE.seq)/sqrt(K)
	    }
	    return(list(alpha=alpha,beta=beta,cv.MSE=cv.MSE,cv.sd=cv.sd))

    }
}





rnc.logistic <- function(A,lambda,Y,X=NULL,gamma=0.05,max.iter=20,tol=1e-4,init=NULL,cv=NULL,cv.seed=999,verbose=FALSE){
    #print(paste("max iteration number:",max.iter))
    n <- nrow(A)
    if(is.null(X)){
    	    p <- 0
    	}else{
    		p <- ncol(X)
    		}
    if(is.null(init)){
        init <- matrix(0,nrow=n+p,ncol=1)
    }

    D <- diag(rowSums(A))
    L <- D-A+gamma*diag(rep(1,n))
    if(is.null(X)){
    	theta <- rnc_logistic_fit_noX(L=L,Y=Y,lambda=lambda,theta_init = init,tol=tol,max_iter=max.iter,verbose=verbose)
    }else{
        theta <- rnc_logistic_fit(X=X,L=L,Y=Y,lambda=lambda,theta_init = init,tol=tol,max_iter=max.iter,verbose=verbose)
    }

    alpha <- theta[1:n]
    beta <- theta[-(1:n)]
    cv.dev <- NULL
    cv.sd <- NULL
    if(!is.null(cv)){
        K <- cv
        set.seed(cv.seed)
        cv.dev <- 0
        cv.order <- sample(n,size=n)
        cv.index <- 1:n
        cv.index[cv.order] <- 1:K
        cv.dev.seq <- rep(0,K)
        #print("Beging cross-validation!")
        for(k in 1:K){
        	    #print(paste(k,"th cross-validation...."))
            current.index <- which(cv.index!=k)
            valid.index <- which(cv.index==k)
            s.A <- A[current.index,current.index]
            if(is.null(X)){
	            cv.logit.net <- rnc.logistic(X=NULL,A=s.A,Y=matrix(Y[current.index],ncol=1),lambda=lambda,init = matrix(init[-valid.index,1],ncol=1),tol=tol,max.iter=max.iter,verbose=FALSE,cv=NULL)

		        L11 <- Matrix(L[valid.index,valid.index],sparse=TRUE)
		        L12 <- -1*L[valid.index,current.index]
		        valid.alpha <- solve(L11,L12%*%cv.logit.net$alpha,sparse=TRUE)
	            valid.eta <- valid.alpha

            }else{
	            cv.logit.net <- rnc.logistic(X=matrix(X[current.index,],ncol=ncol(X)),A=s.A,Y=matrix(Y[current.index],ncol=1),lambda=lambda,init = matrix(init[-valid.index,1],ncol=1),tol=tol,max.iter=max.iter,verbose=FALSE,cv=NULL)
	            n.valid <- length(valid.index)
		        L11 <- Matrix(L[valid.index,valid.index],sparse=TRUE)
		        L12 <- -1*L[valid.index,current.index]
		        valid.alpha <- solve(L11,L12%*%cv.logit.net$alpha,sparse=TRUE)
	            valid.eta <- valid.alpha + X[valid.index,]%*%matrix(cv.logit.net$beta,ncol=1)
            }
            valid.p <- exp(valid.eta)/(1+exp(valid.eta))
            valid.p[valid.p<1e-6] <- 1e-6
            valid.p[valid.p>(1-1e-6)] <- 1-1e-6
            cv.dev.seq[k] <-  - mean(Y[valid.index]*log(valid.p)) - mean((1-Y[valid.index])*log(1-valid.p))
            #print(cv.dev.seq[k])
    }
        cv.dev <- mean(cv.dev.seq)
        cv.sd <- sd(cv.dev.seq)/sqrt(K)
    }


    return(list(alpha=alpha,beta=beta,cv.dev=cv.dev,cv.sd=cv.sd))
}




rnc.cox <- function(A,lambda,dt,X=NULL,max.iter=500,tol=5e-6,init=NULL,gamma=0.05,cv=NULL,cv.seed=999,verbose=FALSE){
    n <- nrow(A)
    if(!is.null(X)){
        if(class(X)!="matrix") stop("X should be a matrix!")
        if(nrow(X)!=n) stop("Dimensions do not match!")
    }
    if(nrow(dt)!=n) stop("Dimensions do not match!")
    if(class(dt)!="data.frame") stop("dt must be a data.frame!")
    if(gamma<=1e-5) stop("The input gamma has to be larger than 1e-5!")

    if(is.null(X)){
    	    p <- 0
    	}else{
    		p <- ncol(X)
    		}
    if(is.null(init)){

        alpha.null <- rep(0,n)
        beta.null <- rep(0,p)
        theta.init <- matrix(c(alpha.null,beta.null),ncol=1)
    }else{
        theta.init <- matrix(init,ncol=1)
        if(is.null(X)) theta.init <- theta.init[1:n,1]
    }

    D <- diag(rowSums(A))
    L <- ((D - A) + gamma*diag(rep(1,n)))
    Y <- matrix(dt$y,ncol=1)
    delta <- matrix(dt$delta,ncol=1)
    if(!is.null(X)){
    	theta <-rnc_cox_fit(X=X,L=L,Y=Y,delta=delta,lambda=lambda,theta_init=theta.init,tol=tol,max_iter=max.iter,verbose=verbose)
    }else{
    	theta <-rnc_cox_fit_noX(L=L,Y=Y,delta=delta,lambda=lambda,theta_init=theta.init,tol=tol,max_iter=max.iter,verbose=verbose)
    	}
    alpha <- theta[1:n]
    beta <- theta[-(1:n)]
    if(is.null(X)){
    	eta <- matrix(alpha,ncol=1)
    	}else{
    		eta <- matrix(alpha + X%*%beta,ncol=1)
    	}
    lpl <- cox_pll(eta=eta,Y=matrix(dt$y,ncol=1),delta = matrix(dt$delta,ncol=1))
    CV.lpl <- NULL
    cv.sd <- NULL
if(!is.null(cv)){
set.seed(cv.seed)
K <- cv
if(K<n){
dead.index <- which(dt$delta==1)
dead.cv.index <- sample(x=rep(1:K,ceiling(length(dead.index)/K))[1:length(dead.index)],size=length(dead.index))
censor.index <- which(dt$delta==0)
censor.cv.index <- sample(x=rep(1:K,ceiling(length(censor.index)/K))[1:length(censor.index)],size=length(censor.index))

cv.index <- rep(0,n)
cv.index[dead.index] <- dead.cv.index
cv.index[censor.index] <- censor.cv.index
}
if(K==n){
    cv.index <- 1:n
}

CV.lpl <- SCV.lpl <- 0
CV.lpl.seq <- rep(0,K)
for(cv.k in 1:K){
        #print(paste("LOO CV at iteration ",cv.k))
        valid.index <- which(cv.index==cv.k)
        #print(valid.index)
	    current.index <- seq(1,n)[-valid.index]
        s.A <- A[current.index,current.index]
	    s.dt <- dt[current.index,]
	    if(is.null(X)){
	    	    s.X <- NULL
	    }else{
	        s.X <- matrix(X[current.index,],ncol=p)
	        }
        tmp.init <- NULL
        if(!is.null(init)) {
        	    if(is.null(X)){
        	    	tmp.init <- init[c(current.index)]
        	    }else{
        	        tmp.init <- init[c(current.index,(n+1):(n+p))]
        	    }
        	}
        tmp.model <- rnc.cox(A=s.A,lambda=lambda,dt=s.dt,X=s.X,max.iter=max.iter,gamma=gamma,init=tmp.init,tol=tol)
        #print("CV model fitted!")
        L.valid <- Matrix(L[valid.index,valid.index],sparse=TRUE)

        valid.alpha <- -solve(L.valid,L[valid.index,current.index]%*%tmp.model$alpha)
        new.alpha <- rep(0,n)
        new.alpha[current.index] <- tmp.model$alpha
        new.alpha[valid.index] <- as.numeric(valid.alpha)
        if(is.null(X)){
    	    cv.eta <- matrix(new.alpha,ncol=1)
    	}else{
    		cv.eta <- matrix(new.alpha + X%*%tmp.model$beta,ncol=1)
    	}
        l.beta.k <- cox_pll(eta=cv.eta,Y=matrix(dt$y,ncol=1),delta = matrix(dt$delta,ncol=1))
        l.k.beta.k <- tmp.model$lpl
        CV.lpl.seq[cv.k] <- l.beta.k - l.k.beta.k

    }
    CV.lpl <- mean(CV.lpl.seq)
    cv.sd <- sd(CV.lpl.seq)/sqrt(K)

}
    return(list(alpha=alpha,beta=beta,theta=theta,lpl = lpl,lambda=lambda,gamma=gamma,cv=cv,cv.dev=-1*CV.lpl,cv.sd=cv.sd))
}



# # rnc.cox.path <- function(A,lambda.seq,dt,X=NULL,max.iter=500,tol=5e-6,init=NULL,gamma=0.05,cv=NULL,cv.seed=999,verbose=FALSE){
    # model.list <- list()
    # n.lambda <- length(lambda.seq)
    # for(i in 1:n.lambda){
    	# lambda <- lambda.seq[i]
    	# model.list[[i]] <- rnc.cox(A=A,lambda=lambda,dt=dt,X=X,max.iter=max.iter,tol=tol,init=init,gamma=gamma,cv=cv,cv.seed=cv.seed,verbose=verbose)
    # }
    # cv.seq <- unlist(lapply(model.list,function(x) x$cv.dev))

    # return(list(models=model.list,cv.seq=cv.seq))
# }


## function to fit regression model with network cohesion
rncreg <- function(A,lambda,Y=NULL,X=NULL,dt=NULL,gamma=0.05,model=c("linear","logistic","cox"),max.iter=50,tol=1e-4,init=NULL,cv=NULL,cv.seed=999,low.dim=NULL,verbose=FALSE){
	if(length(model)>1) model <- model[1]
	if(!any(c("linear","logistic","cox")==model)) stop("model must be one of linear, logistic and cox.")
    n <- nrow(A)
    if(norm(t(A)-A,"F")>0) warning("A is not symmetric....")

    if(!is.null(X)){
        if(class(X)!="matrix") stop("X should be a matrix!")
        if(nrow(X)!=n) stop("Dimensions do not match!")
    }
    if(!is.null(Y)){
        if(class(Y)!="matrix") stop("Y should be a matrix!")
        if(nrow(Y)!=n) stop("Dimensions do not match!")
    }
    if(!is.null(dt)){
        if(nrow(dt)!=n) stop("Dimensions do not match!")
        if(class(dt)!="data.frame") stop("dt must be a data.frame!")
    }

    if(model=="linear"){
    	    result <- rnc.linear(X=X,Y=Y,A=A,lambda=lambda,gamma=gamma,cv=cv,cv.seed=cv.seed,low.dim=low.dim)
    	    result.obj <- list(alpha=result$alpha,beta=result$beta,A=A,lambda=lambda,X=X,Y=Y,dt=dt,gamma=gamma,cv=cv,cv.loss=result$cv.MSE,cv.sd=result$cv.sd,model=model)
    	    class(result.obj) <- "rncReg"
    	    return(result.obj)
    }

    if(model=="logistic"){
    	    result <- rnc.logistic(X=X,Y=Y,A=A,lambda=lambda,gamma=gamma,cv=cv,cv.seed=cv.seed,max.iter=max.iter,tol=tol,init=init,verbose=verbose)
    	    result.obj <- list(alpha=result$alpha,beta=result$beta,A=A,lambda=lambda,X=X,Y=Y,dt=dt,gamma=gamma,cv=cv,cv.loss=result$cv.dev,cv.sd=result$cv.sd,model=model)
    	    class(result.obj) <- "rncReg"
    	    return(result.obj)
    }

    if(model=="cox"){
    	    result <- rnc.cox(X=X,dt=dt,A=A,lambda=lambda,gamma=gamma,cv=cv,cv.seed=cv.seed,max.iter=max.iter,tol=tol,init=init,verbose=verbose)
    	    #result.obj <- list(alpha=result$alpha,beta=result$beta,A=A,lambda=lambda,X=X,dt=dt,gamma=gamma,cv=cv,cv.dev=result$cv.dev)
    	    #class(result.obj) <- "rncCox"
    	    result.obj <- list(alpha=result$alpha,beta=result$beta,A=A,lambda=lambda,X=X,Y=Y,dt=dt,gamma=gamma,cv=cv,cv.loss=result$cv.dev,cv.sd=result$cv.sd,model=model)
        class(result.obj) <- "rncReg"
    	    return(result.obj)
    }

}


rncregpath <- function(A,lambdaseq,Y=NULL,X=NULL,dt=NULL,gamma=0.05,model=c("linear","logistic","cox"),max.iter=50,tol=1e-4,init=NULL,cv=NULL,cv.seed=999,low.dim=NULL,verbose=FALSE){
	N <- length(lambdaseq)
	models <- list()
	for(k in 1:N){
		models[[k]] <- rncreg(A=A,lambda=lambdaseq[k],Y=Y,X=X,dt=dt,gamma=gamma,model=model,max.iter=max.iter,tol=tol,init=init,cv=cv,cv.seed=cv.seed,low.dim=low.dim,verbose=verbose)
	}
	cv.seq <- unlist(lapply(models,function(x) x$cv.loss))
        cv.sd.seq <- unlist(lapply(models,function(x) x$cv.sd))
        cv.min.index <- which.min(cv.seq)
        thres <- min(cv.seq) + cv.sd.seq[cv.min.index]
        se.lambda <- max(lambdaseq[cv.seq<=thres])
        cv.1sd.index <- which(lambdaseq==se.lambda)
	return(list(models=models,cv.seq=cv.seq,cv.sd=cv.sd.seq,cv.min.index=cv.min.index,cv.1sd.index=cv.1sd.index))
}



predict_rncLinear <- function(obj,full.X=NULL,full.A){
	    Adj <- full.A
	    n <- nrow(Adj)
	    if(!is.null(full.X)){
	        p <- ncol(full.X)
	        }
	    D <- diag(rowSums(Adj))
	    gamma <- obj$gamma
	    L <- D - Adj + diag(rep(gamma,n))
	    W <- diag(rep(1,n))
        	alpha <- obj$alpha
	    beta <- obj$beta
	    lambda <- obj$lambda
	    n.train <- nrow(obj$A)
	    if(n <= n.train) stop("Dimension mismatch!")
	    if(!is.null(full.X)){
	    		    valid.index <- (n.train+1):n
	            current.index <- 1:n.train

	            L11 <- Matrix(L[valid.index,valid.index],sparse=TRUE)
	            L12 <- -1*L[valid.index,current.index]
	            valid.alpha <- solve(L11,L12%*%alpha,sparse=TRUE)
	            valid.y <- valid.alpha + full.X[valid.index,]%*%beta
                return(list(y=as.matrix(valid.y),terms=as.matrix(valid.y),alpha=as.matrix(valid.alpha)))
	    }else{
	    		    valid.index <- (n.train+1):n
	            current.index <- 1:n.train

	            L11 <- Matrix(L[valid.index,valid.index],sparse=TRUE)
	            L12 <- -1*L[valid.index,current.index]
	            valid.alpha <- solve(L11,L12%*%alpha,sparse=TRUE)
	            valid.y <- valid.alpha
                return(list(y=as.matrix(valid.y),terms=as.matrix(valid.y),alpha=as.matrix(valid.alpha),model="linear"))

	    }

}


predict_rncLogistic <- function(obj,full.X=NULL,full.A){
	    Adj <- full.A
	    n <- nrow(Adj)
	    if(!is.null(full.X)){
	        p <- ncol(full.X)
	        }
	    D <- diag(rowSums(Adj))
	    gamma <- obj$gamma
	    L <- D - Adj + diag(rep(gamma,n))
	    W <- diag(rep(1,n))
        	alpha <- obj$alpha
	    beta <- obj$beta
	    lambda <- obj$lambda
	    n.train <- nrow(obj$A)
	    if(n <= n.train) stop("Dimension mismatch!")
	    if(!is.null(full.X)){
	    		    valid.index <- (n.train+1):n
	            current.index <- 1:n.train

	            L11 <- Matrix(L[valid.index,valid.index],sparse=TRUE)
	            L12 <- -1*L[valid.index,current.index]
	            valid.alpha <- solve(L11,L12%*%alpha,sparse=TRUE)
	            valid.terms <- valid.alpha + full.X[valid.index,]%*%beta
	            valid.p <- exp(valid.terms)/(1+exp(valid.terms))
                return(list(p=as.matrix(valid.p),terms=as.matrix(valid.terms),alpha=as.matrix(valid.alpha)))
	    }else{
	    		    valid.index <- (n.train+1):n
	            current.index <- 1:n.train

	            L11 <- Matrix(L[valid.index,valid.index],sparse=TRUE)
	            L12 <- -1*L[valid.index,current.index]
	            valid.alpha <- solve(L11,L12%*%alpha,sparse=TRUE)
	            valid.terms <- valid.alpha
	            valid.p <- exp(valid.terms)/(1+exp(valid.terms))
                return(list(p=as.matrix(valid.p),terms=as.matrix(valid.terms),alpha=as.matrix(valid.alpha),model="logistic"))

	    }

}




predict_rncCox <- function(obj,full.X=NULL,full.A){
	    Adj <- full.A
	    n <- nrow(Adj)
	    if(!is.null(full.X)){
	        p <- ncol(full.X)
	        }
	    D <- diag(rowSums(Adj))
	    gamma <- obj$gamma
	    L <- D - Adj + diag(rep(gamma,n))
	    W <- diag(rep(1,n))
        	alpha <- obj$alpha
	    beta <- obj$beta
	    lambda <- obj$lambda
	    n.train <- nrow(obj$A)
	    if(n <= n.train) stop("Dimension mismatch!")
	    if(!is.null(full.X)){
	    		    valid.index <- (n.train+1):n
	            current.index <- 1:n.train

	            L11 <- Matrix(L[valid.index,valid.index],sparse=TRUE)
	            L12 <- -1*L[valid.index,current.index]
	            valid.alpha <- solve(L11,L12%*%alpha,sparse=TRUE)
	            valid.terms <- valid.alpha + full.X[valid.index,]%*%beta

                return(list(terms=as.matrix(valid.terms),alpha=as.matrix(valid.alpha)))
	    }else{
	    		    valid.index <- (n.train+1):n
	            current.index <- 1:n.train

	            L11 <- Matrix(L[valid.index,valid.index],sparse=TRUE)
	            L12 <- -1*L[valid.index,current.index]
	            valid.alpha <- solve(L11,L12%*%alpha,sparse=TRUE)
	            valid.terms <- valid.alpha
                return(list(terms=as.matrix(valid.terms),alpha=as.matrix(valid.alpha),model="cox"))
	    }

}


predict.rncReg <- function(object,full.X=NULL,full.A,...){
	if(object$model=="linear"){		return(predict_rncLinear(obj=object,full.X=full.X,full.A=full.A))
	}
	if(object$model=="logistic"){		return(predict_rncLogistic(obj=object,full.X=full.X,full.A=full.A))
	}
	if(object$model=="cox"){		return(predict_rncCox(obj=object,full.X=full.X,full.A=full.A))
	}else{
		stop("Invalid model type!")
	}

}




