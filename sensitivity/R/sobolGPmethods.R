plot.sobolGP <- function(x, ...){
	sGP <- x
	dev.new(width = 12, height = 5)
	if(sGP$call$tot){
		par(mfrow=c(1,3))
	} else {
		par(mfrow=c(1,2))
	}
	d <- length(sGP$S$mean)
	ylim=c(-0.05,1.05)
	xlim=c(1,d+0.1)
	plot(1:d,sGP$S$mean, xlim = xlim, ylim = ylim, axes = FALSE,
       		xlab = "", ylab = "", pch = 21)
	axis(side = 1, at = 1:d, labels = 1:d)
	axis(side = 2, at = seq(0,1,0.1), labels = seq(0,1,0.1))
	segments(1:d,sGP$S$ci[1,],1:d,sGP$S$ci[2,])
	box()
	if(sGP$call$tot){
		points((1:d)+0.1,sGP$T$mean, pch = 24)	
		legend("topright", legend = c("main effect", "total effect"), pch = c(21,24))
		segments(1:d+0.1,sGP$T$ci[1,],1:d+0.1,sGP$T$ci[2,],lty=2)
	} else {
		legend("topright", legend = c("main effect"), pch = 21)
	}
	abline(0,0)

nboot <- sGP$call$nboot
if(nboot!=1){
	data <- data.frame(rbind(sGP$S$varPG,sGP$S$varMC))
	names(data) <- names(sGP$S)[1:d]
     	barplot(as.matrix(data), angle = c(45,-45), density = 20, col = "black",
             legend = c("GP variance","MC variance"),ylim=c(0,1.2*max(sGP$S$varPG+sGP$S$varMC)))
     	title(main = list("Variance decomposition of the main effects", font = 4))
	
	if(sGP$call$tot){
		data <- data.frame(rbind(sGP$T$varPG,sGP$T$varMC))
		names(data) <- names(sGP$T)[1:d]
     		barplot(as.matrix(data), angle = c(45,-45), density = 20, col = "black",
       		      legend = c("GP variance","MC variance"),ylim=c(0,1.2*max(sGP$T$varPG+sGP$T$varMC)))
     		title(main = list("Variance decomposition of the total effects", font = 4))
	}
} else {
	data <- data.frame(sGP$S$varPG)
	names(data) <- names(sGP$S)[1:d]
     	barplot(as.matrix(data), angle = 45, density = 20, col = "black",
             legend = c("GP variance"),ylim=c(0,1.2*max(sGP$S$varPG)))
     	title(main = list("Variance of the main effects", font = 4))
	
	if(sGP$call$tot){
		data <- data.frame(sGP$T$varPG)
		names(data) <- names(sGP$T)[1:d]
     		barplot(as.matrix(data), angle = 45, density = 20, col = "black",
       		      legend = c("GP variance"),ylim=c(0,1.2*max(sGP$T$varPG)))
     		title(main = list("Variance of the total effects", font = 4))
	}
}
	
}  

print.sobolGP <- function(x, ...){
	sGP <- x
	cat("\nMethod: ", sGP$call$method, "\n", sep = "")
	cat("\nModel runs:", sGP$call$model@n, "\n")
	cat("\nNumber of GP realizations:", sGP$call$nsim, "\n")
	cat("\nKriging type:", sGP$call$type, "\n")

	d <- length(sGP$S$mean)
	df=data.frame(pointEstimate=t(sGP$S$mean), stdErr=t(sqrt(sGP$S$var)), minCI=as.matrix(sGP$S$ci[1,]), maxCI=as.matrix(sGP$S$ci[2,]))
	colnames(df)=c("estimate","std. error","min. c.i.","max. c.i.")
	rownames(df)=names(sGP$S)[1:d]
	cat("\n")
	print(df)
	cat("\n")
	if(sGP$call$tot){
		df=data.frame(pointEstimate=t(sGP$T$mean), stdErr=t(sqrt(sGP$T$var)), minCI=as.matrix(sGP$T$ci[1,]), maxCI=as.matrix(sGP$T$ci[2,]))
		colnames(df)=c("estimate","std. error","min. c.i.","max. c.i.")
		rownames(df)=names(sGP$T)[1:d]
		print(df)
		cat("\n")
	}
}

ask.sobolGP <- function(x,tot=FALSE, ...){
	sGP <- x
	output <- c()
	if(is.null(sGP$S$xnew)){
		stop("The function sobolGP must be run with the argument sequential = TRUE")
	}
	if(!tot){
		output <- sGP$S$xnew
	} else {
		if(is.null(sGP$T$xnew)){
			stop("The function sobolGP must be run with the argument sequential.tot = TRUE")
		}	
		output <- sGP$T$xnew
	}
	d <- dim(sGP$call$candidate)[2]
	return(output)
}

tell.sobolGP <- function(x,y=NULL,xpoint=NULL,newcandidate=NULL, ...){
	sGP <- x
	x <- xpoint
	if(is.null(newcandidate)){
		d <- dim(sGP$call$candidate)[2]
		test = c(sGP$call$candidate[,1] == x[,1])
		for(i in 2:d){
			test <- test*(sGP$call$candidate[,i] == x[,i])
		}
		newcandidate <- sGP$call$candidate[-which(test == 1),]
	}

	x <- as.matrix(x)
	if(dim(x)[2]==1){
		x <- t(x)
	}
	rownames(x) = NULL
	m  <- update(object = sGP$call$model, newX = x, newy = y, newX.alreadyExist = FALSE,cov.reestim = FALSE)

	

	res <- sobolGP(
		model = m,
		type = sGP$call$type,
		MCmethod = sGP$call$method,
		X1 = sGP$call$X1,
		X2 = sGP$call$X2,
		nsim = sGP$call$nsim,
		conf = sGP$call$conf,
		nboot = sGP$call$nboot,
		sequential = sGP$call$sequential,
		candidate = newcandidate,
		sequential.tot = sGP$call$sequential.tot,
		max_iter = sGP$call$max_iter) 

	return(res)

}

predictGP.sobol <- function (
object,
newdata1 ,
newdata2,
type,
bias.correct = FALSE,
checkNames = FALSE, 
prednewdata1 = FALSE, 
prednewdata2 = TRUE) 
{
    nugget.flag <- object@covariance@nugget.flag
    X <- object@X
    y <- object@y
    if (checkNames) {
        newdata1 <- checkNames(X1 = X, X2 = newdata1, X1.name = "the design", 
            X2.name = "newdata")
        newdata2 <- checkNames(X1 = X, X2 = newdata2, X1.name = "the design", 
            X2.name = "newdata")
    } else {
        newdata1 <- as.matrix(newdata1)
        newdata2 <- as.matrix(newdata2)
        d.newdata <- ncol(newdata1)
        if (!identical(d.newdata, object@d)) {
            stop("newdata must have the same numbers of columns than the experimental design")
        }
        if (!identical(colnames(newdata1), colnames(X))) {
            colnames(newdata1) <- colnames(X)
        }
        if (!identical(colnames(newdata2), colnames(X))) {
            colnames(newdata2) <- colnames(X)
        }
    }
    T <- object@T
    z <- object@z
    M <- object@M
    beta <- object@trend.coef

    F.newdata1 <- model.matrix(object@trend.formula, data = data.frame(newdata1))
    F.newdata2 <- model.matrix(object@trend.formula, data = data.frame(newdata2))

    if(prednewdata1){
    	y.predict.trend1 <- F.newdata1 %*% beta
    } 
    if(prednewdata2){
    	y.predict.trend2 <- F.newdata2 %*% beta
    }
    if (requireNamespace("DiceKriging", quietly = TRUE)){ 
      c.newdata1 <- DiceKriging::covMat1Mat2(object@covariance, X1 = X, X2 = newdata1, 
        nugget.flag = object@covariance@nugget.flag)
      c.newdata2 <- DiceKriging::covMat1Mat2(object@covariance, X1 = X, X2 = newdata2, 
        nugget.flag = object@covariance@nugget.flag)
    }
    Tinv.c.newdata1 <- backsolve(t(T), c.newdata1, upper.tri = FALSE)
    Tinv.c.newdata2 <- backsolve(t(T), c.newdata2, upper.tri = FALSE)

    output.list <- list()
    if(prednewdata1){
    	y.predict.complement1 <- t(Tinv.c.newdata1) %*% z
    	y.predict1 <- y.predict.trend1 + y.predict.complement1
	rm(list=c("y.predict.complement1","y.predict.trend1"))
   	y.predict1 <- as.numeric(y.predict1)
   	output.list$mean1 <- y.predict1	
	rm(list=c("y.predict1"))
    }
    if(prednewdata2){
   	y.predict.complement2 <- t(Tinv.c.newdata2) %*% z
   	y.predict2 <- y.predict.trend2 + y.predict.complement2
	rm(list=c("y.predict.complement2","y.predict.trend2"))
   	y.predict2 <- as.numeric(y.predict2)
    	output.list$mean2 <- y.predict2
	rm(list=c("y.predict2"))
    }
    
	  if (requireNamespace("DiceKriging", quietly = TRUE)){
      C.newdata <- DiceKriging::covMat1Mat2(object@covariance, X1 = newdata1, X2 = newdata2, 
      nugget.flag = object@covariance@nugget.flag)
	  }

	rm(list=c("newdata1","newdata2"))

    cond.cov <- C.newdata - crossprod(Tinv.c.newdata1,Tinv.c.newdata2 )
 
    if (type == "UK") {
        T.M <- chol(t(M) %*% M)

        s2.predict.mat1 <- backsolve(t(T.M), t(F.newdata1 - 
            t(Tinv.c.newdata1) %*% M), upper.tri = FALSE)
        s2.predict.mat2 <- backsolve(t(T.M), t(F.newdata2 - 
            t(Tinv.c.newdata2) %*% M), upper.tri = FALSE)

        cond.cov <- cond.cov + crossprod(s2.predict.mat1,s2.predict.mat2)
        if (bias.correct) 
            cond.cov <- cond.cov * object@n/(object@n - object@p)
    }
    output.list$cov <- cond.cov

    return(output.list)
}

simulateGP.sobol <- function(
object,  
nsim = 100, 
newdata = NULL, 
cond = TRUE, 
checkNames = FALSE, 
max_iter = 1000,
type
) {
	
  if (!is.logical(cond)) stop("'cond' must be TRUE/FALSE")
  if ((!is.null(newdata)) && (checkNames)) {
    newdata <- checkNames(X1 = object@X, X2 = newdata, X1.name = "the design", X2.name = "newdata")
  }

  if (is.null(newdata)) {
    newdata <- object@X
    F.newdata <- object@F
  } else {
    newdata <- as.matrix(newdata)
    m <- nrow(newdata)
        
    if (!identical(ncol(newdata), object@d)) 
      stop("newdata must have the same numbers of columns than the experimental design")
    if (!identical(colnames(newdata), colnames(object@X))) {
      colnames(newdata) <- colnames(object@X)
    }

    newdata <- rbind(newdata, object@X)
    row.names(newdata) <- NULL
    F.newdata <- model.matrix(object@trend.formula, data = data.frame(newdata))
  }
  
  
  n <- object@n
  
# non conditional simulations
    sample_test <- sample.int(m,min(100,m))
    iter <- 1
    test <- FALSE
    yc <- matrix(0,nsim,m+n)
	#RMSE_save <- c() #-!!-#

    while((!test&iter<1000)){  
		a <- sample.int(m+n, size = nsim, replace = FALSE, prob = NULL)
		ya <- rnorm(nsim,0,1)
		if (requireNamespace("DiceKriging", quietly = TRUE)){
      Cab <- DiceKriging::covMat1Mat2(object@covariance, X1 = as.matrix(newdata[a,]), X2 = as.matrix(newdata), 
        		nugget.flag = object@covariance@nugget.flag)/object@covariance@sd2
		}
		yc <- yc + Cab* (as.matrix(ya-diag(yc[,a]))%*%t(as.matrix(rep(1,m+n))))
		diag(yc[1:nsim,a]) <- ya

		if(iter%%20==0){
			Cemp <- var(yc[,sample_test]) 
			Ctheo <- DiceKriging::covMatrix(object@covariance, X = as.matrix(newdata[sample_test,]))[[1]]/
					object@covariance@sd2 
			RMSE_C <- sqrt(mean((Cemp-Ctheo)^2))  
			#RMSE_save <- c(RMSE_save,RMSE_C) #-!!-#
			if(RMSE_C < 0.1){test <- TRUE} 
		}

		iter <- iter+1
    }  
	rm(list=c("Cemp","Ctheo","ya" ,"Cab","a" ,"sample_test"))


	yc <- yc*sqrt(object@covariance@sd2*object@n/(object@n-object@p))

  if (cond){	
	if(type=="UK"){
		y.predict.trend <- F.newdata %*% object@trend.coef
		if (requireNamespace("DiceKriging", quietly = TRUE)){
      c.newdata <- DiceKriging::covMat1Mat2(object@covariance, X1 = object@X, X2 = newdata, 
	          nugget.flag = object@covariance@nugget.flag)
		}

		rm(list=c("newdata"))	

	      Tinv.c.newdata <- backsolve(t(object@T), c.newdata, upper.tri = FALSE)

		rm(list=c("c.newdata"))		

		y.predict.complement <- t(Tinv.c.newdata) %*% object@z
		y.predict <- y.predict.trend + y.predict.complement

		rm(list=c("y.predict.complement"))

		y.predict <- as.numeric(y.predict)
	
		ZD <- t(yc[,(m+1):(m+n)])
		if(!identical(object@noise.var,numeric(0))){
			ZD <- ZD + matrix(rnorm(nsim*n,0,1),nrow=n)*matrix(rep(sqrt(object@noise.var),nsim),ncol=nsim)
		}
	
		xtilde <- backsolve(t(object@T), ZD, upper.tri = FALSE)
		Betatilde <- solve(t(object@M)%*%object@M	)%*%t(object@M)%*%xtilde 

		rm(list=c("xtilde"))
	
		ztilde <-  backsolve(t(object@T), ZD -  object@F%*%Betatilde  , 
	            upper.tri = FALSE)
	
		rm(list=c("ZD"))
	
	      ytilde.predict.complement <- t(Tinv.c.newdata) %*% ztilde

		rm(list=c("ztilde","Tinv.c.newdata"))

	      ytilde.predict.trend <- F.newdata %*% Betatilde  

		rm(list=c("F.newdata","Betatilde"))
	
		ytilde.predict <- ytilde.predict.trend + ytilde.predict.complement
	
		rm(list=c("ytilde.predict.trend","ytilde.predict.complement","y.predict.trend"))
	
		yc <- t(matrix(rep(y.predict,nsim),ncol=nsim)) - t(ytilde.predict) +yc
		rm(list=c("y.predict","ytilde.predict"))
	}
	if(type=="SK"){
		y.predict.trend <- F.newdata %*% object@trend.coef
		if (requireNamespace("DiceKriging", quietly = TRUE)){
      c.newdata <- DiceKriging::covMat1Mat2(object@covariance, X1 = object@X, X2 = newdata, 
	          nugget.flag = object@covariance@nugget.flag)
		}
	      Tinv.c.newdata <- backsolve(t(object@T), c.newdata, upper.tri = FALSE)
		y.predict.complement <- t(Tinv.c.newdata) %*% object@z
		y.predict <- y.predict.trend + y.predict.complement
		y.predict <- as.numeric(y.predict)

		ZD <- t(yc[,(m+1):(m+n)])
		if(!identical(object@noise.var,numeric(0))){
			ZD <- ZD + matrix(rnorm(nsim*n,0,1),nrow=n)*matrix(rep(sqrt(object@noise.var),nsim),ncol=nsim)
		}

		ztilde <-  backsolve(t(object@T), ZD  ,  upper.tri = FALSE)
		ytilde.predict <- t(Tinv.c.newdata) %*% ztilde

		rm(list=c("ztilde","c.newdata","Tinv.c.newdata","newdata","y.predict.trend","F.newdata","ZD"))

		yc <- t(matrix(rep(y.predict,nsim),ncol=nsim)) - t(ytilde.predict) +yc
		rm(list=c("y.predict","ytilde.predict"))
	}
  }
  
   return(yc[,-c((m+1):(m+n))])  

}










