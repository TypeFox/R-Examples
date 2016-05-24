bw.pi <- function(x,M=NULL,lower=0,upper=100,np=500,tol=0.1,outM=FALSE){
	x <- conversion.circular(x, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
	attr(x, "class") <- attr(x, "circularp") <- NULL
	if (!is.numeric(x)) stop("argument 'x' must be numeric")
	x.na <- is.na(x)
	if (sum(x.na)>0) warning("Missing values were removed")
	x <- x[!x.na]
	n <- length(x)
	if (n==0) stop("No observations (at least after removing missing values)")
	if (!is.null(M)){  
		if (M!=as.integer(M) | M<1){ 
			warning("The provided value of 'M' is not valid, 'M' must be a positive integer. 'M' was selected by AIC")
			M=NULL
		}
	}
	if (!is.numeric(upper)){ 
		warning("argument 'upper' must be numeric. Default upper boundary was used")
		upper <- 100
	}
	if (!is.numeric(lower)){
		warning("argument 'lower' must be numeric. Default lower boundary was used")
		lower <- 0
	}
	if (lower<0 | lower>=upper){
      	warning("The boundaries must be positive and 'lower' must be smaller that 'upper'. Default boundaries were used")
		upper <- 100
		lower <- 0
	}
	if (!is.numeric(np)) stop("argument 'np' must be numeric")
	if (np<=0){ 
		warning("'np' must be positive. Default value of 'np' was used")
		np <- 500
	}
	t <- seq(0,2*pi,length=np)
	if (!is.numeric(tol)) stop("argument 'tol' must be numeric")
	z <- cbind(cos(x),sin(x))
	if (is.null(M)){ 	
		x <- circular(x)
		##  Mixture of M=2 components
		y2 <- try(movMF(z, 2, start="S"),TRUE)
		if (class(y2)=="try-error"){
			AIC2 <- NaN
			mean2 <- kappa2 <- p2 <- NA
		}else{ 
			norm2<-mean2<-kappa2<-p2<-numeric(2)
			mu2<-matrix(NA,2,2)
			AIC2<-0
			for (i in 1:2){
				norm2[i]<-sqrt(sum(y2$theta[i,]^2))
				mu2[i,]<-y2$theta[i,]/norm2[i]
				mean2[i]<-atan2(mu2[i,2],mu2[i,1])
				kappa2[i]<-y2$theta[i,1]/mu2[i,1]
				p2[i]<-y2$alpha[i]
				AIC2 <- AIC2+p2[i]*dvonmises(x,circular(mean2[i]),kappa2[i])
			}
			AIC2 <- -2 *sum(log(AIC2)) + 10
		}
		##  Mixture of M=3 components
		y3 <- try(movMF(z, 3, start="S"),TRUE)
		if (class(y3)=="try-error"){
			AIC3 <- NaN
			mean3 <- kappa3 <- p3 <- NA
		}else{ 
			norm3<-mean3<-kappa3<-p3<-numeric(3)
			mu3<-matrix(NA,3,2)
			AIC3<-0
			for (i in 1:3){
				norm3[i]<-sqrt(sum(y3$theta[i,]^2))
				mu3[i,]<-y3$theta[i,]/norm3[i]
				mean3[i]<-atan2(mu3[i,2],mu3[i,1])
				kappa3[i]<-y3$theta[i,1]/mu3[i,1]
				p3[i]<-y3$alpha[i]
				AIC3 <- AIC3+p3[i]*dvonmises(x,circular(mean3[i]),kappa3[i])
			}
			AIC3 <- -2 *sum(log(AIC3)) + 16
		}
		##  Mixture of M=4 components
		y4 <- try(movMF(z, 4, start="S"),TRUE)
		if (class(y4)=="try-error"){
			AIC4 <- NaN
			mean4 <- kappa4 <- p4 <- NA
		}else{ 
			norm4<-mean4<-kappa4<-p4<-numeric(4)
			mu4<-matrix(NA,4,2)
			AIC4<-0
			for (i in 1:4){
				norm4[i]<-sqrt(sum(y4$theta[i,]^2))
				mu4[i,]<-y4$theta[i,]/norm4[i]
				mean4[i]<-atan2(mu4[i,2],mu4[i,1])
				kappa4[i]<-y4$theta[i,1]/mu4[i,1]
				p4[i]<-y4$alpha[i]
				AIC4 <- AIC4+p4[i]*dvonmises(x,circular(mean4[i]),kappa4[i])
			}
			AIC4 <- -2 *sum(log(AIC4)) + 22
		}
		##  Mixture of M=5 components
		y5 <- try(movMF(z, 5, start="S"),TRUE)
		if (class(y5)=="try-error"){
			AIC5 <- NaN
			mean5 <- kappa5 <- p5 <- NA
		}else{ 
			norm5<-mean5<-kappa5<-p5<-numeric(5)
			mu5<-matrix(NA,5,2)
			AIC5<-0
			for (i in 1:5){
				norm5[i]<-sqrt(sum(y5$theta[i,]^2))
				mu5[i,]<-y5$theta[i,]/norm5[i]
				mean5[i]<-atan2(mu5[i,2],mu5[i,1])
				kappa5[i]<-y5$theta[i,1]/mu5[i,1]
				p5[i]<-y5$alpha[i]
				AIC5 <- AIC5+p5[i]*dvonmises(x,circular(mean5[i]),kappa5[i])
			}
			AIC5 <- -2 *sum(log(AIC5)) + 28
		}

		AICs <- c(AIC2,AIC3,AIC4,AIC5)
		mean <- list(mean2,mean3,mean4,mean5)
		kappa <- list(kappa2,kappa3,kappa4,kappa5)
		p <- list(p2,p3,p4,p5)

		f2<- function(t,M) { 
			der2f<-numeric(length(t))
			for (i in 1:M){
				der2f<-der2f+p[[M-1]][i]*kappa[[M-1]][i]*(kappa[[M-1]][i]*exp(kappa[[M-1]][i]*cos(t-mean[[M-1]][i]))*(sin(t-mean[[M-1]][i]))^2 - 
				exp(kappa[[M-1]][i]*cos(t-mean[[M-1]][i]))*cos(t-mean[[M-1]][i]))/(2*pi*besselI(kappa[[M-1]][i],0))
			}
			return(der2f^2)
		}

		x<-as.numeric(x)
		if (class(y2)=="try-error" & class(y3)=="try-error" & class(y4)=="try-error" & class(y5)=="try-error"){
			M <- 1
			bw <- bw.rt(x)
			warning("The smoothing parameter was computed by using the rule of thumb")
		}else{
			AMISE <- function(bw) {1/16 * (1- besselI(bw,2,expon.scaled=T)/besselI(bw,0,expon.scaled=T))^2 * integral + besselI(2*bw,0,expon.scaled=TRUE)/(2*n*pi*(besselI(bw,0,expon.scaled=TRUE))^2)}
			M<-which.min(AICs)+1
			if (length(M) == 0){
				M <- 1
				bw <- bw.rt(x)
				warning("The smoothing parameter was computed by using the rule of thumb")
			}else{
				integral <- int.Simp(t,f2(t,M))
				if(integral=="NaN" | integral=="Inf"){
					while(integral=="NaN" | integral=="Inf"){
						AICs[M-1] <- NaN
						M<-which.min(AICs)+1
						if (length(M)==0) {integral<-0
						}else{ integral <- int.Simp(t,f2(t,M))}	
					}
					if (length(M)==0){
						M <- 1
						bw <- bw.rt(x)
						warning("The smoothing parameter was computed by using the rule of thumb")
					}else{
						bw <- optimize(function(h) AMISE(h),interval=c(lower,upper),tol=tol)$minimum
					}
				}else{
					bw <- optimize(function(h) AMISE(h),interval=c(lower,upper),tol=tol)$minimum
				}
			}
		}
	}else{    
		y2 <- try(movMF(z, M, start="S"),TRUE)
		if (M==1){ 
			bw <- bw.rt(x)
			M <- 1
			warning("The smoothing parameter was computed by using the rule of thumb")
		}else{
			if (class(y2)=="try-error"){		
				bw <- bw.rt(x)
				M <- 1
				warning("The smoothing parameter was computed by using the rule of thumb")
			}else{
				norm<-mean<-kappa<-p<-numeric(M)
				mu<-matrix(NA,M,2)
				for (i in 1:M){
					norm[i]<-sqrt(sum(y2$theta[i,]^2))
					mu[i,]<-y2$theta[i,]/norm[i]
					mean[i]<-atan2(mu[i,2],mu[i,1])
					kappa[i]<-y2$theta[i,1]/mu[i,1]
					p[i]<-y2$alpha[i]
				}
				f2.Mfix<- function(t) {
					der2f<-numeric(length(t))
					for (i in 1:M){
						der2f<-der2f+p[i]*kappa[i]*(kappa[i]*exp(kappa[i]*cos(t-mean[i]))*(sin(t-mean[i]))^2 - exp(kappa[i]*cos(t-mean[i]))*cos(t-mean[i]))/(2*pi*besselI(kappa[i],0))
					}
					return(der2f^2)
				}
				f2t <- f2.Mfix(t)
				integral <- int.Simp(t,f2t)
				if (integral=="Inf" | integral=="NaN"){
					bw <- bw.rt(x)
				}else{
					AMISE <- function(bw) {1/16 * (1- besselI(bw,2,expon.scaled=TRUE)/besselI(bw,0,expon.scaled=TRUE))^2 * integral  + besselI(2*bw,0,expon.scaled=TRUE)/(2*n*pi*(besselI(bw,0,expon.scaled=TRUE))^2)}
					bw <- optimize(function(h) AMISE(h),interval=c(lower,upper),tol=tol)$minimum
				}
			}
		}
	}
	if (bw < lower + tol | bw > upper - tol) 
      warning("minimum occurred at one end of the range")
	if (outM) return(c(bw,M))
	else return(bw)
}
