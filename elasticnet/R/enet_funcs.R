enet<-
function (x, y, lambda = 0, max.steps, normalize = TRUE, intercept = TRUE, 
    trace = FALSE, eps = .Machine$double.eps) 
{
    call <- match.call()
    nm <- dim(x)
    n <- nm[1]
    m <- nm[2]
    im <- seq(m)
    one <- rep(1, n)
    vn <- dimnames(x)[[2]]
    meanx <- drop(one %*% x)/n
    if (intercept == FALSE) {
        meanx <- rep(0, m)
    }
    x <- scale(x, meanx, FALSE)
    normx <- sqrt(drop(one %*% (x^2)))
    if (normalize == FALSE) {
        normx <- rep(1, m)
    }
    if (any(normx < eps * sqrt(n))) 
        stop("Some of the columns of x have zero variance")
    names(normx) <- NULL
    x <- scale(x, FALSE, normx)
    mu <- mean(y)
    if (intercept == FALSE) {
        mu <- 0
    }
    y <- drop(y - mu)
    d1 <- sqrt(lambda)
    d2 <- 1/sqrt(1 + lambda)
    Cvec <- drop(t(y) %*% x) * d2
    ssy <- sum(y^2)
    residuals <- c(y, rep(0, m))
    if (lambda > 0) {
        maxvars <- m
    }
    if (lambda == 0) {
        maxvars <- min(m, n - 1)
    }
    if (missing(max.steps)) {
        max.steps <- 50 * maxvars
    }
    L1norm <- 0
    penalty <- max(abs(Cvec))
    beta <- rep(0, m)
    betactive <- list(NULL)
    first.in <- integer(m)
    active <- NULL
    Actset <- list(NULL)
    df <- 0
    if (lambda != 0) {
        Cp <- ssy
    }
    ignores <- NULL
    actions <- as.list(seq(max.steps))
    drops <- FALSE
    Sign <- NULL
    R <- NULL
    k <- 0
    while ((k < max.steps) & (length(active) < maxvars)) {
        action <- NULL
        k <- k + 1
        inactive <- if (k == 1) 
            im
        else im[-c(active, ignores)]
        C <- Cvec[inactive]
        Cmax <- max(abs(C))
        if (!any(drops)) {
            new <- abs(C) == Cmax
            C <- C[!new]
            new <- inactive[new]
            for (inew in new) {
                R <- updateRR(x[, inew], R, x[, active], lambda)
                if (attr(R, "rank") == length(active)) {
                  nR <- seq(length(active))
                  R <- R[nR, nR, drop = FALSE]
                  attr(R, "rank") <- length(active)
                  ignores <- c(ignores, inew)
                  action <- c(action, -inew)
                  if (trace) 
                    cat("LARS-EN Step", k, ":\t Variable", inew, 
                      "\tcollinear; dropped for good\n")
                }
                else {
                  if (first.in[inew] == 0) 
                    first.in[inew] <- k
                    active <- c(active, inew)
                    Sign <- c(Sign, sign(Cvec[inew]))
                    action <- c(action, inew)
                  if (trace) 
                    cat("LARS-EN Step", k, ":\t Variable", inew, 
                      "\tadded\n")
                }
            }
        }
        else action <- -dropid
        Gi1 <- backsolve(R, backsolvet(R, Sign))
        A <- 1/sqrt(sum(Gi1 * Sign))
        w <- A * Gi1
        u1 <- drop(x[, active, drop = FALSE] %*% w * d2)
        u2 <- rep(0, m)
        u2[active] <- d1 * d2 * w
        u <- c(u1, u2)
        if (lambda > 0) {
            maxvars <- m - length(ignores)
        }
        if (lambda == 0) {
            maxvars <- min(m - length(ignores), n - 1)
        }
        if (length(active) >= maxvars) {
            gamhat <- Cmax/A
        }
        else {
            a <- (drop(u1 %*% x[, -c(active, ignores)]) + d1 * 
                u2[-c(active, ignores)]) * d2
            gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
            gamhat <- min(gam[gam > eps], Cmax/A)
            Cdrop <- c(C - gamhat * a, -C + gamhat * a) - (Cmax - 
                gamhat * A)
        }
        dropid <- NULL
        b1 <- beta[active]
        z1 <- -b1/w
        zmin <- min(z1[z1 > eps], gamhat)
        if (zmin < gamhat) {
            gamhat <- zmin
            drops <- z1 == zmin
        }
        else drops <- FALSE
        beta[active] <- beta[active] + gamhat * w

        betactive[[k]] <- beta[active]
        Actset[[k]] <- active
        residuals <- residuals - (gamhat * u)
        Cvec <- (drop(t(residuals[1:n]) %*% x) + d1 * residuals[-(1:n)]) * d2
        L1norm <- c(L1norm, sum(abs(beta[active]))/d2)
        penalty <- c(penalty, penalty[k] - abs(gamhat * A))
        if (any(drops)) {
            dropid <- seq(drops)[drops]
            for (id in rev(dropid)) {
                if (trace) 
                  cat("LARS-EN Step", k, ":\t Variable", active[id], 
                    "\tdropped\n")
                  R <- downdateR(R, id)
            }
            dropid <- active[drops]
            beta[dropid] <- 0
            active <- active[!drops]
            Sign <- Sign[!drops]
        }
        if (!is.null(vn)) 
            names(action) <- vn[abs(action)]
        actions[[k]] <- action
    }
    allset <- Actset[[1]]
    for (i in 2:k) {
        allset <- union(allset, Actset[[i]])
    }
    allset <- sort(allset)
    max.p <- length(allset)
    beta.pure <- matrix(0, k + 1, max.p)
    for (i in 2:(k + 1)) {
        for (j in 1:length(Actset[[i - 1]])) {
            l <- c(1:max.p)[allset == Actset[[i - 1]][j]]
            beta.pure[i, l] <- betactive[[i - 1]][j]
        }
    }
    beta.pure <- beta.pure/d2
    dimnames(beta.pure) <- list(paste(0:k), vn[allset])
    k <- dim(beta.pure)[1]
    df <- 1:k
    for (i in 1:k) {
        a <- drop(beta.pure[i, ])
        df[i] <- 1 + length(a[a != 0])
    }
    residuals <- y - x[, allset, drop = FALSE] %*% t(beta.pure)
    beta.pure <- scale(beta.pure, FALSE, normx[allset])
    RSS <- apply(residuals^2, 2, sum)
    R2 <- 1 - RSS/RSS[1]
    Cp <- ((n - m - 1) * RSS)/rev(RSS)[1] - n + 2 * df
    object <- list(call = call, actions = actions[seq(k)], allset = allset, 
        beta.pure = beta.pure, vn = vn, mu = mu, normx = normx[allset], 
        meanx = meanx[allset], lambda = lambda, L1norm = L1norm, 
        penalty = penalty * 2/d2, df = df, Cp = Cp, sigma2 = rev(RSS)[1]/(n - 
            m - 1))
    class(object) <- "enet"
    object
}


predict.enet<-
function (object, newx, s, type = c("fit", "coefficients"), mode = c("step", 
    "fraction", "norm", "penalty"), naive = FALSE, ...) 
{
    mode <- match.arg(mode)
    type <- match.arg(type)
    if (missing(newx) & type == "fit") {
        warning("Type=fit with no newx argument; type switched to coefficients")
        type <- "coefficients"
    }
    betas <- object$beta.pure
    if (naive) {
        lambda <- object$lambda
        betas <- betas/(1 + lambda)
    }
    allset <- object$allset
    sbetas <- scale(betas, FALSE, 1/object$normx)
    kp <- dim(betas)
    k <- kp[1]
    p <- kp[2]
    steps <- seq(k)
    if (missing(s)) {
        s <- steps
        mode <- "step"
    }
    sbeta <- switch(mode, step = {
        if (any(s < 0) | any(s > k)) 
            stop("Argument s out of range")
        steps
    }, fraction = {
        if (any(s > 1) | any(s < 0)) 
            stop("Argument s out of range")
        nbeta <- drop(abs(sbetas) %*% rep(1, p))
        nbeta/nbeta[k]
    }, norm = {
        nbeta <- drop(abs(sbetas) %*% rep(1, p))
        if (any(s > nbeta[k]) | any(s < 0)) 
            stop("Argument s out of range")
        nbeta
    }, penalty = {
        pen <- object$penalty
        if (any(s > pen[1]) | any(s < 0)) 
            stop("Argument s out of range")
        pen
    })
    sfrac <- (s - sbeta[1])/(sbeta[k] - sbeta[1])
    sbeta <- (sbeta - sbeta[1])/(sbeta[k] - sbeta[1])
    usbeta <- unique(sbeta)
    useq <- match(usbeta, sbeta)
    sbeta <- sbeta[useq]
    betas <- betas[useq, ]
    coord <- approx(sbeta, seq(sbeta), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    newbetas <- ((sbeta[right] - sfrac) * betas[left, , drop = FALSE] + 
        (sfrac - sbeta[left]) * betas[right, , drop = FALSE])/(sbeta[right] - 
        sbeta[left])
    newbetas[left == right, ] <- betas[left[left == right], ]
    robject <- switch(type, coefficients = list(s = s, fraction = sfrac, 
        mode = mode, coefficients = drop(newbetas)), fit = list(s = s, 
        fraction = sfrac, mode = mode, fit = drop(scale(newx[, 
            allset], object$meanx, FALSE) %*% t(newbetas)) + 
            object$mu))
    robject
}

plot.enet<-
function (x, xvar=c("fraction","penalty","L1norm","step"),use.color = FALSE, ...) 
{
    object <- x
    xvar <- match.arg(xvar)
    vn <- object$vn
    allset <- object$allset
    vn <- vn[allset]
    coef1 <- object$beta
    coef1 <- scale(coef1, FALSE, 1/object$normx)
    cnums <- seq(ncol(coef1))
    m <- nrow(coef1)
    p <- ncol(coef1)
    s1 <- switch(xvar, fraction = {
        s1 <- object$L1norm
        s1 <- s1/max(s1)
    }, penalty = {
        s1 <- object$penalty
    }, L1norm = {
        s1 <- object$L1norm
    }, step = {
        s1 <- 1:m
    })
    xname <- switch(xvar, fraction = "|beta|/max|beta|", penalty = "lambda (Lasso)", L1norm = "L1norm", step = "Step")
            low <- min(coef1) - 0.01 * abs(min(coef1))
            up <- max(coef1) + 0.01 * abs(max(coef1))
            ylimit <- c(low, up)
        plot(s1, s1, xlab = xname, ylab = "Standardized Coefficients",
             ylim = ylimit, type = "n", ...)
        for (i in 1:p) {
            if (use.color) {
                lines(s1, coef1[, i], col = i, lty = 1)
            }
            else {
                lines(s1, coef1[, i], lty = 1)
            }
            if (!is.null(vn)) {
             axis(4, at = coef1[nrow(coef1), ], labels = vn, cex = 0.8, adj = 0)
            }
          }
        abline(h = 0, lty = 3)

    invisible()
}

print.enet<-      
function(x, ...)   
{ 
        cat("\nCall:\n")
        dput(x$call)
        if(x$lambda==0){
           cat("Cp statistics of the Lasso fit", "\n")
           cat("Cp:", format(round(x$Cp, 3)), "\n")
           cat("DF:", format(round(x$df,3)), "\n")
        }
                                                         
        actions <- x$actions                                                              
        jactions <- unlist(actions)
        jsteps <- rep(seq(along = actions), sapply(actions, length))
        actmat <- rbind(jsteps, jactions)
        vn <- names(jactions)
        if(is.null(vn))
                vn <- rep("", length(jactions))                        
        dimnames(actmat) <- list(c("Step", "Var"), vn)
        cat(paste("Sequence of", x$type, "moves:\n"))
        print(actmat[2:1,  ])
        invisible(x)
}

cv.enet<-
function(x, y, K = 10, lambda, s, mode, 
           trace = FALSE, plot.it = TRUE, se = TRUE, ...)
{
  all.folds <- cv.folds(length(y), K)
  residmat <- matrix(0, length(s), K)
  for(i in seq(K)) {
    omit <- all.folds[[i]]
    fit <- enet(x[ - omit,  ], y[ - omit], lambda=lambda, ...)
    fit <- predict(fit, x[omit,  ,drop=FALSE], s=s, mode = mode)$fit
    if(length(omit)==1){fit<-matrix(fit,nrow=1)}
    residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
    if(trace)
      cat("\n CV Fold", i, "\n\n")
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object<-list(s = s, cv = cv, cv.error = cv.error)
  if(plot.it) {
     plot(s, cv, type = "b", xlab=mode, ylim = range(cv, cv + cv.error, cv - cv.error))
     if (se) 
        error.bars(s, cv + cv.error, cv - cv.error, width = 1/length(s))
  }
  invisible(object)
}


updateRR<-
function(xnew, R = NULL, xold, lambda,eps = .Machine$double.eps)
{
	xtx <- (sum(xnew^2)+lambda)/(1+lambda)
	norm.xnew <- sqrt(xtx)
	if(is.null(R)) {
		R <- matrix(norm.xnew, 1, 1)
		attr(R, "rank") <- 1
		return(R)
	}
	Xtx <- drop(t(xnew) %*% xold)/(1+lambda)
	r <- backsolvet(R, Xtx)
	rpp <- norm.xnew^2 - sum(r^2)
	rank <- attr(R, "rank")	
	if(rpp <= eps)
		rpp <- eps
	else {
		rpp <- sqrt(rpp)
		rank <- rank + 1
	}
	R <- cbind(rbind(R, 0), c(r, rpp))
	attr(R, "rank") <- rank
	R
}


spca<-
function(x,K,para,type=c("predictor","Gram"),sparse=c("penalty","varnum"),use.corr=FALSE, 
         lambda=1e-6,max.iter=200,trace=FALSE,eps.conv=1e-3)
{
      call <- match.call()
      type <- match.arg(type)
      sparse <- match.arg(sparse)
      vn <- dimnames(x)[[2]]
      x<-switch(type,
                predictor = {
                n<-dim(x)[1]
                p<-dim(x)[2]
                if (n/p>=100){
    cat("You may wish to restart and use a more efficient way \n")
    cat("let the argument x be the sample covariance/correlation matrix and set type=Gram \n")
    }
                 x<-scale(x,center=TRUE,scale=use.corr)
                      },
                Gram = {x<-rootmatrix(x)}
                )
      
      svdobj<-svd(x)
      v<-svdobj$v
      totalvariance<-sum((svdobj$d)^2)
      alpha<-as.matrix(v[,1:K,drop=FALSE])      
      beta<-alpha      
      for ( i in 1:K) {
         y<-drop(x%*%alpha[,i])
         beta[,i]<-solvebeta(x,y,paras=c(lambda,para[i]),sparse=sparse)
       }
      xtx<-t(x)%*%x
      temp<-beta
      normtemp<-sqrt(apply(temp^2,2,sum))
      normtemp[normtemp==0]<-1
      temp<-t(t(temp)/normtemp)
      k<-0
      diff<-1
      while((k<max.iter) & (diff>eps.conv)){
         k<-k+1
         alpha<-xtx%*%beta
         z<-svd(alpha)
         alpha<-(z$u)%*%t(z$v)
         for ( i in 1:K) {
            y<-drop(x%*%alpha[,i])
            beta[,i]<-solvebeta(x,y,paras=c(lambda,para[i]),sparse=sparse)       
         }
         normbeta<-sqrt(apply(beta^2,2,sum))
         normbeta[normbeta==0]<-1
         beta2<-t(t(beta)/normbeta)
         diff<-max(abs(beta2-temp))
         temp<-beta2
           if(trace){
              if (k%%10==0){
                  cat("iterations",k,fill=TRUE)
                 }
            }
         }
      normbeta<-sqrt(apply(beta^2,2,sum))
      normbeta[normbeta==0]<-1
      beta<-t(t(beta)/normbeta)
      dimnames(beta)<-list(vn,paste("PC",1:K,sep=""))
      u<-x%*%beta
      R<-qr.R(qr(u))
      pev<-diag(R^2)/totalvariance
      obj<-list(call = call, type=type, K=K,loadings=beta,pev=pev,var.all=totalvariance, vn=vn,para=para,lambda=lambda)
      class(obj) <- "spca"
      obj  
}

solvebeta<-function(x, y, paras, max.steps, sparse=c("penalty","varnum"), eps = .Machine$double.eps)
{
       sparse <- match.arg(sparse)
       if (missing(sparse)) sparse <- "penalty"
	nm <- dim(x)
	n <- nm[1]
	m <- nm[2]
	im <- seq(m)
	one <- rep(1, n)
	vn <- dimnames(x)[[2]]
        lambda<-paras[1]
        if(lambda>0){
	   maxvars <- m
        }
        if (lambda==0) {
           maxvars <- min(m,n-1)
           if (m==n){
             maxvars<-m
           }
        }
	d1 <- sqrt(lambda)
	d2 <- 1/sqrt(1 + lambda)
	Cvec <- drop(t(y) %*% x) * d2
        ssy <- sum(y^2)
	residuals <- c(y, rep(0, m))
        if(missing(max.steps)) {max.steps <- 50 * maxvars}
        penalty<-max(abs(Cvec))
        if (sparse=="penalty" && penalty*2/d2<=paras[2])
          { 
           beta<-rep(0,m)
          }
        else {
               beta <- rep(0, m)
               first.in <- integer(m)
               active <- NULL
               ignores <- NULL
               actions <- as.list(seq(max.steps))
               drops <- FALSE
               Sign <- NULL
               R <- NULL
	       k <- 0
	       while((k < max.steps) & (length(active) < maxvars - length(ignores))) 
		{
		 action <- NULL
		 k <- k + 1
        	 inactive <- if(k == 1) im else im[ - c(active, ignores)]
		 C <- Cvec[inactive]
		 Cmax <- max(abs(C))
		 if(!any(drops)) {
			new <- abs(C) == Cmax 
			C <- C[!new]
			new <- inactive[new]       
			for(inew in new) {
                        R <- updateRR(x[, inew], R, x[, active], lambda) 
                               if(attr(R, "rank") == length(active)) {
					nR <- seq(length(active))
					R <- R[nR, nR, drop = FALSE]
					attr(R, "rank") <- length(active)
					ignores <- c(ignores, inew)
					action <- c(action,  - inew)
				}
				else {
					if(first.in[inew] == 0)
						first.in[inew] <- k
					active <- c(active, inew)
					Sign <- c(Sign, sign(Cvec[inew]))
					action <- c(action, inew)
				}
			}
		}
		else action <-  - dropid
                Gi1 <- backsolve(R, backsolvet(R, Sign))
            	A <- 1/sqrt(sum(Gi1 * Sign))
            	w <- A * Gi1
                u1<-drop(x[,active,drop=FALSE]%*%w*d2)  
                u2<-rep(0,m)
                u2[active]<-d1*d2*w
                u<-c(u1,u2)
                if(lambda>0){
	            maxvars <- m-length(ignores)
                           }
                if (lambda==0){
                    maxvars <- min(m-length(ignores),n-1)
                           }
       		if(length(active) == maxvars - length(ignores)) {
			gamhat <- Cmax/A
		}
		else {
			a <- (drop(u1 %*% x[,  - c(active, ignores)]) + d1 *
				u2[ - c(active, ignores)]) * d2
      			gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
				gamhat <- min(gam[gam > eps], Cmax/A)                        
			Cdrop <- c(C - gamhat * a,  - C + gamhat * a) - (Cmax -
				gamhat * A)
		}
		dropid <- NULL
		b1 <- beta[active]
		z1 <-  - b1/w
		zmin <- min(z1[z1 > eps], gamhat)
		if(zmin < gamhat) {
			gamhat <- zmin
			drops <- z1 == zmin
		}
               else drops <- FALSE
               beta2<-beta          
               beta[active] <- beta[active] + gamhat * w
               residuals <- residuals - (gamhat*u)
               Cvec <- (drop(t(residuals[1:n]) %*% x) + d1 * residuals[ - (1:n)]) * d2
               penalty <- c(penalty,penalty[k]-abs(gamhat*A))
               if(sparse=="penalty" && rev(penalty)[1]*2/d2<=paras[2]){
                   s1<-rev(penalty)[1]*2/d2
                   s2<-rev(penalty)[2]*2/d2
                   beta<-(s2-paras[2])/(s2-s1)*beta+(paras[2]-s1)/(s2-s1)*beta2
                   beta<-beta*d2
                   break
               }
		if(any(drops)) {
			dropid <- seq(drops)[drops]
			for(id in rev(dropid)) {
				R <- downdateR(R, id)
      			}
			dropid <- active[drops]
                        beta[dropid] <- 0
			active <- active[!drops]
			Sign <- Sign[!drops]
                      }
               if(sparse=="varnum" && length(active)>=paras[2]){
                 break
               }
	       if(!is.null(vn))
			names(action) <- vn[abs(action)]
                	actions[[k]] <- action 
	}
             }
     
        return (beta)
}

rootmatrix<-function(x){
     x.eigen<-eigen(x)
     d<-x.eigen$values
     d<-(d+abs(d))/2
     v<-x.eigen$vectors
     return (v%*%diag(sqrt(d))%*%t(v))
} 


print.spca<-function(x, ...){
     cat("\nCall:\n")
     dput(x$call)
     K<-x$K
     cat("\n")
     cat(paste(K,"sparse PCs",sep=" "), "\n")
     cat("Pct. of exp. var. :", format(round(x$pev, 3)*100), "\n")
     v<-x$loadings
     k<-1:K
     for (j in 1:K){
        k[j]<-sum(v[,j]!=0)
     }
     cat("Num. of non-zero loadings : ",k,"\n")
     cat("Sparse loadings \n")
     print(round(x$loadings,3))
     type<-x$type
     if (type=="penalty"){
       cat("lambda_1:", format(round(x$para, 3)), "\n")
     }
}


arrayspc<-
function(x,K=1,para,use.corr=FALSE, max.iter=100,trace=FALSE,eps=1e-3)
{     call <- match.call()
      x<-scale(x,center=TRUE,scale=use.corr)
      svdobj<-svd(x)
      v<-svdobj$v
      totalvariance<-sum((svdobj$d)^2)      
      alpha<-as.matrix(v[,1:K,drop=FALSE])      
      beta<-alpha      
      for ( i in 1:K) {
         y<-drop(x%*%alpha[,i])
         beta[,i]<-soft(drop(t(x)%*%y),para=para[i])
         }
      temp<-beta
      normtemp<-sqrt(apply(temp^2,2,sum))
      normtemp[normtemp==0]<-1
      temp<-t(t(temp)/normtemp)
      k<-0
      diff<-1
      while((k<max.iter) & (diff>eps)){
         k<-k+1
         alpha<-x%*%beta
         alpha<-t(x)%*%alpha
         z<-svd(alpha)
         alpha<-(z$u)%*%t(z$v)
         for ( i in 1:K) {
           y<-drop(x%*%alpha[,i])
           beta[,i]<-soft(drop(t(x)%*%y),para=para[i])
         }     
         normbeta<-sqrt(apply(beta^2,2,sum))
         normbeta[normbeta==0]<-1
         beta2<-t(t(beta)/normbeta)
         diff<-max(abs(beta2-temp))
         temp<-beta2
          if(trace){
              if (k%%10==0){
                  cat("iterations",k,fill=TRUE)
                 }
            }
      } 
      normbeta<-sqrt(apply(beta^2,2,sum))
      normbeta[normbeta==0]<-1
      beta<-t(t(beta)/normbeta)
      u<-x%*%beta
      R<-qr.R(qr(u))
      pev<-diag(R^2)/totalvariance
      obj<-list(call = call, K=K,loadings=beta,pev=pev,var.all=totalvariance,para=para)
      class(obj) <- "arrayspc"
      obj 
}

soft<-function(a,para){
  b<-sort(abs(a))
  b<-abs(a)-para
  b<-(b+abs(b))/2
  b<-sign(a)*b
  b
}

print.arrayspc<-function(x, ...){
     cat("\nCall:\n")
     dput(x$call)
     K<-x$K
     cat("\n")
     cat(paste(K,"sparse PCs",sep=" "), "\n")
     cat("Pct. of exp. var. :", format(round(x$pev, 3)*100), "\n")
     v<-x$loadings
     k<-1:K
     for (j in 1:K){
        k[j]<-sum(v[,j]!=0)
     }
     cat("Num. of non-zero loadings : ",k,"\n")
}

require("lars")
