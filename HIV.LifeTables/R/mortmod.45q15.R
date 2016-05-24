mortmod.45q15 <- function(child.mort, adult.mort, prev, region=1, sex=1, opt=TRUE){
	
lt.mx <-
function(nmx, sex="female", age=c(0,1,seq(5,100,5)), nax=NULL){
	if(is.null(nax)){

		nax <- rep(2.5,length(age))

		if (sex=="male") { # TRUE for male, FALSE for female
			if (nmx[1] >= 0.107) {
				nax[1] <- 0.33
				nax[2] <- 1.352
			} else {
				nax[1] <- 0.045+2.684*nmx[1]
				nax[2] <- 1.651-2.816*nmx[1]
				}
			} 
		if (sex=="female") {
			if (nmx[1] >= 0.107) {
				nax[1] <- 0.35
				nax[2] <- 1.361
			} else {
				nax[1] <- 0.053+2.8*nmx[1]
				nax[2] <- 1.522-1.518*nmx[1]
				}
			}
	} else {nax=nax}
    
    n <- c(diff(age), 999)
    nqx <- (n*nmx)/(1+(n-nax)*nmx)
    nqx <- c(nqx[-(length(nqx))], 1)
    for(i in 1:length(nqx)){
    	if(nqx[i] > 1) nqx[i] <- 1
    	}
    nage <- length(age)

    #nqx <- round(nqx, 4)

    npx <- 1 - nqx
    max.age <- min(which(npx==0))
    l0=100000
    lx <- round(cumprod(c(l0, npx)))
    #lx <- (cumprod(c(l0, npx)))
    ndx <- -diff(lx)
    lxpn <- lx[-1]
    nLx <- n * lxpn + ndx * nax
    lx <- lx[1:length(age)]
    nLx[max.age] <- lx[max.age]/nmx[max.age]
    Tx <- rev(cumsum(rev(nLx)))
    ex <- Tx/lx
    lt <- cbind(Age = age, nax = nax, nmx = nmx, nqx = nqx, npx = npx, ndx = ndx, lx = lx, nLx = round(nLx), Tx = round(Tx), 
        ex = round(ex, 2))
     lt <- lt[lt[,6]!=0,]
     e0 <- lt[1,10]
     lt.45q15 <- 1-(lx[14]/lx[5])
     lt.5q0 <- 1-(lx[3]/lx[1])
     return(list(e0=e0, lt.5q0=lt.5q0, lt.45q15=lt.45q15, lt=lt))
        }	
	
if(region==1 & sex==0){ # Africa, male
	intercept <- median(svd.coeffs.xp[africa.nums,1])
	b1.m <- predict.lm(co1.45q15mod.a.m, newdata=data.frame(cmort.a.m=child.mort, amort.a.m=adult.mort))
	b2.m <- predict.lm(co2.45q15mod.a.m, newdata=data.frame(cmort.a.m=child.mort, amort.a.m=adult.mort, prev.a=prev))
	b3.m <- predict.lm(co3.45q15mod.a.m, newdata=data.frame(cmort.a.m=child.mort, amort.a.m=adult.mort, prev.a=prev))
	
	if(opt==FALSE){
	out.mort <- intercept + b1.m*Mx.svd.scores[,1] + b2.m*Mx.svd.scores[,2] + b3.m*Mx.svd.scores[,3]} # if
	
# optimize the intercent when predicting the mortality rates from the weights
	if(opt==TRUE){
		out.mort.func <- function(intercept.alter){
		out.mort <- intercept.alter + b1.m*Mx.svd.scores[,1] + b2.m*Mx.svd.scores[,2] + b3.m*Mx.svd.scores[,3]
		lt.out <- lt.mx(nmx=exp(out.mort[1:22]), sex="male", age=c(0,1,seq(5,100,5)))
		amort.diff <- abs(adult.mort-unname(lt.out$lt.45q15))
	return(amort.diff)	
	} # function to be optimized
	
	intercept.opt <- optimize(f=out.mort.func, interval=c(-2,2))$minimum
	
	out.mort <- intercept.opt + b1.m*Mx.svd.scores[,1] + b2.m*Mx.svd.scores[,2] + b3.m*Mx.svd.scores[,3]} # if

	
## optimize reduction or addition to 1q0 and 4q1 to make it match
if(opt==TRUE){
	out.mort.start <- out.mort
	lt.out <- lt.mx(nmx=exp(out.mort.start[1:22]), sex="male", age=c(0,1,seq(5,100,5)))

	diff.weight.opt <- log(1-child.mort)/(log(1-lt.out$lt[1,4])+log(1-lt.out$lt[2,4]))
	qx.child.new.opt <- 1-(c(lt.out$lt[1,5]^diff.weight.opt, lt.out$lt[2,5]^diff.weight.opt))
	cmort.opt <- 1-((1-qx.child.new.opt[1])*(1-qx.child.new.opt[2]))
	
	axs <- c(1-lt.out$lt[1,2], 4-lt.out$lt[2,2])
	mx1 <- qx.child.new.opt[1]/(1-axs[1]*qx.child.new.opt[1])
	mx2 <- qx.child.new.opt[2]/(4-axs[2]*qx.child.new.opt[2])

	mx.child.new.opt <- c(mx1, mx2)
	out.mort.new.opt <- unname(c(log(mx.child.new.opt), out.mort.start[3:22]))
	if(out.mort.new.opt[3]>out.mort.new.opt[2]){
		# y <- out.mort.new.opt[1:2]
		# x <- c(1,2)
		# child.mod.lm <- lm(y~x)
	   # new.5q5 <- predict.lm(child.mod.lm, newdata=data.frame(x = 3))
	   new.5m5 <- approx(x=c(2,4),y=out.mort.new.opt[c(2,4)],xout=3)$y
	   out.mort.new.opt[3] <- new.5m5
	}
	lt.out.new.opt <- lt.mx(nmx=exp(out.mort.new.opt), sex="male", age=c(0,1,seq(5,100,5)), nax=c(lt.out$lt[,2], rep(2.5, 22-length(lt.out$lt[,2]))))
	out.mort <- c(out.mort.new.opt, out.mort.start[-c(1:22)])
} # if


	return(exp(out.mort[1:22]))
}

if(region==1 & sex==1){ # Africa, female
	intercept <- median(svd.coeffs.xp[africa.nums,1])
	b1.f <- predict.lm(co1.45q15mod.a.f, newdata=data.frame(cmort.a.f=child.mort, amort.a.f=adult.mort))
	b2.f <- predict.lm(co2.45q15mod.a.f, newdata=data.frame(cmort.a.f=child.mort, amort.a.f=adult.mort, prev.a=prev))
	b3.f <- predict.lm(co3.45q15mod.a.f, newdata=data.frame(cmort.a.f=child.mort, amort.a.f=adult.mort, prev.a=prev))
	
	if(opt==FALSE){
	out.mort <- intercept + b1.f*Mx.svd.scores[,1] + b2.f*Mx.svd.scores[,2] + b3.f*Mx.svd.scores[,3]} # if
	
# optimize the intercent when predicting the mortality rates from the weights
	if(opt==TRUE){
		out.mort.func <- function(intercept.alter){
		out.mort <- intercept.alter + b1.f*Mx.svd.scores[,1] + b2.f*Mx.svd.scores[,2] + b3.f*Mx.svd.scores[,3]
		lt.out <- lt.mx(nmx=exp(out.mort[(22+1):(22*2)]), sex="female", age=c(0,1,seq(5,100,5)))
		amort.diff <- abs(adult.mort-unname(lt.out$lt.45q15))
	return(amort.diff)
	} # function to be optimized
	
	intercept.opt <- optimize(f=out.mort.func, interval=c(-2,2))$minimum
	
	out.mort <- intercept.opt + b1.f*Mx.svd.scores[,1] + b2.f*Mx.svd.scores[,2] + b3.f*Mx.svd.scores[,3]} # if
	
## optimize reduction or addition to 1q0 and 4q1 to make it match
if(opt==TRUE){
	out.mort.start <- out.mort
	lt.out <- lt.mx(nmx=exp(out.mort.start[23:44]), sex="female", age=c(0,1,seq(5,100,5)))

	diff.weight.opt <- log(1-child.mort)/(log(1-lt.out$lt[1,4])+log(1-lt.out$lt[2,4]))
	qx.child.new.opt <- 1-(c(lt.out$lt[1,5]^diff.weight.opt, lt.out$lt[2,5]^diff.weight.opt))
	cmort.opt <- 1-((1-qx.child.new.opt[1])*(1-qx.child.new.opt[2]))
	
	axs <- c(1-lt.out$lt[1,2], 4-lt.out$lt[2,2])
	mx1 <- qx.child.new.opt[1]/(1-axs[1]*qx.child.new.opt[1])
	mx2 <- qx.child.new.opt[2]/(4-axs[2]*qx.child.new.opt[2])

	mx.child.new.opt <- c(mx1, mx2)
	out.mort.new.opt <- unname(c(log(mx.child.new.opt), out.mort.start[25:44]))
	if(out.mort.new.opt[3]>out.mort.new.opt[2]){
		# y <- out.mort.new.opt[1:2]
		# x <- c(1,2)
		# child.mod.lm <- lm(y~x)
	   # new.5q5 <- predict.lm(child.mod.lm, newdata=data.frame(x = 3))
	   new.5m5 <- approx(x=c(2,4),y=out.mort.new.opt[c(2,4)],xout=3)$y
	   out.mort.new.opt[3] <- new.5m5
	}
	lt.out.new.opt <- lt.mx(nmx=exp(out.mort.new.opt), sex="female", age=c(0,1,seq(5,100,5)), nax=c(lt.out$lt[,2], rep(2.5, 22-length(lt.out$lt[,2]))))
	out.mort <- c(out.mort.start[1:22], out.mort.new.opt)
} # if


	return(exp(out.mort[(22+1):(22*2)]))
}

if(region==0 & sex==0){ # Non-Africa, male
	intercept <- median(svd.coeffs.xp[la.nums,1])
	b1.m <- predict.lm(co1.45q15mod.na.m, newdata=data.frame(cmort.na.m=child.mort, amort.na.m=adult.mort))
	b2.m <- predict.lm(co2.45q15mod.na.m, newdata=data.frame(cmort.na.m=child.mort, amort.na.m=adult.mort, prev.na=prev))
	b3.m <- predict.lm(co3.45q15mod.na.m, newdata=data.frame(cmort.na.m=child.mort, amort.na.m=adult.mort, prev.na=prev))
	
	if(opt==FALSE){
	out.mort <- intercept + b1.m*Mx.svd.scores[,1] + b2.m*Mx.svd.scores[,2] + b3.m*Mx.svd.scores[,3]} # if
	
# optimize the intercent when predicting the mortality rates from the weights
	if(opt==TRUE){
		out.mort.func <- function(intercept.alter){
		out.mort <- intercept.alter + b1.m*Mx.svd.scores[,1] + b2.m*Mx.svd.scores[,2] + b3.m*Mx.svd.scores[,3]
		lt.out <- lt.mx(nmx=exp(out.mort[1:22]), sex="male", age=c(0,1,seq(5,100,5)))
		amort.diff <- abs(adult.mort-unname(lt.out$lt.45q15))
	return(amort.diff)	
	} # function to be optimized
	
	intercept.opt <- optimize(f=out.mort.func, interval=c(-2,2))$minimum
	
	out.mort <- intercept.opt + b1.m*Mx.svd.scores[,1] + b2.m*Mx.svd.scores[,2] + b3.m*Mx.svd.scores[,3]} # if
	
## optimize reduction or addition to 1q0 and 4q1 to make it match
if(opt==TRUE){
	out.mort.start <- out.mort
	lt.out <- lt.mx(nmx=exp(out.mort.start[1:22]), sex="male", age=c(0,1,seq(5,100,5)))

	diff.weight.opt <- log(1-child.mort)/(log(1-lt.out$lt[1,4])+log(1-lt.out$lt[2,4]))
	qx.child.new.opt <- 1-(c(lt.out$lt[1,5]^diff.weight.opt, lt.out$lt[2,5]^diff.weight.opt))
	cmort.opt <- 1-((1-qx.child.new.opt[1])*(1-qx.child.new.opt[2]))
	
	axs <- c(1-lt.out$lt[1,2], 4-lt.out$lt[2,2])
	mx1 <- qx.child.new.opt[1]/(1-axs[1]*qx.child.new.opt[1])
	mx2 <- qx.child.new.opt[2]/(4-axs[2]*qx.child.new.opt[2])

	mx.child.new.opt <- c(mx1, mx2)
	out.mort.new.opt <- unname(c(log(mx.child.new.opt), out.mort.start[3:22]))
	if(out.mort.new.opt[3]>out.mort.new.opt[2]){
		# y <- out.mort.new.opt[1:2]
		# x <- c(1,2)
		# child.mod.lm <- lm(y~x)
	   # new.5q5 <- predict.lm(child.mod.lm, newdata=data.frame(x = 3))
	   new.5m5 <- approx(x=c(2,4),y=out.mort.new.opt[c(2,4)],xout=3)$y
	   out.mort.new.opt[3] <- new.5m5
	}
	lt.out.new.opt <- lt.mx(nmx=exp(out.mort.new.opt), sex="male", age=c(0,1,seq(5,100,5)), nax=c(lt.out$lt[,2], rep(2.5, 22-length(lt.out$lt[,2]))))
	out.mort <- c(out.mort.new.opt, out.mort.start[-c(1:22)])
} # if

	return(exp(out.mort[1:22]))
}

if(region==0 & sex==1){ # Non-Africa, female
	intercept <- median(svd.coeffs.xp[la.nums,1])
	b1.f <- predict.lm(co1.45q15mod.na.f, newdata=data.frame(cmort.na.f=child.mort, amort.na.f=adult.mort))
	b2.f <- predict.lm(co2.45q15mod.na.f, newdata=data.frame(cmort.na.f=child.mort, amort.na.f=adult.mort, prev.na=prev))
	b3.f <- predict.lm(co3.45q15mod.na.f, newdata=data.frame(cmort.na.f=child.mort, amort.na.f=adult.mort, prev.na=prev))
	
	if(opt==FALSE){
	out.mort <- intercept + b1.f*Mx.svd.scores[,1] + b2.f*Mx.svd.scores[,2] + b3.f*Mx.svd.scores[,3]} # if
	
# optimize the intercent when predicting the mortality rates from the weights
	if(opt==TRUE){
		out.mort.func <- function(intercept.alter){
		out.mort <- intercept.alter + b1.f*Mx.svd.scores[,1] + b2.f*Mx.svd.scores[,2] + b3.f*Mx.svd.scores[,3]
		lt.out <- lt.mx(nmx=exp(out.mort[(22+1):(22*2)]), sex="female", age=c(0,1,seq(5,100,5)))
		amort.diff <- abs(adult.mort-unname(lt.out$lt.45q15))
	return(amort.diff)
	} # function to be optimized
	
	intercept.opt <- optimize(f=out.mort.func, interval=c(-2,2))$minimum
	
	out.mort <- intercept.opt + b1.f*Mx.svd.scores[,1] + b2.f*Mx.svd.scores[,2] + b3.f*Mx.svd.scores[,3]} # if

## optimize reduction or addition to 1q0 and 4q1 to make it match
if(opt==TRUE){
	out.mort.start <- out.mort
	lt.out <- lt.mx(nmx=exp(out.mort.start[23:44]), sex="female", age=c(0,1,seq(5,100,5)))

	diff.weight.opt <- log(1-child.mort)/(log(1-lt.out$lt[1,4])+log(1-lt.out$lt[2,4]))
	qx.child.new.opt <- 1-(c(lt.out$lt[1,5]^diff.weight.opt, lt.out$lt[2,5]^diff.weight.opt))
	cmort.opt <- 1-((1-qx.child.new.opt[1])*(1-qx.child.new.opt[2]))
	
	axs <- c(1-lt.out$lt[1,2], 4-lt.out$lt[2,2])
	mx1 <- qx.child.new.opt[1]/(1-axs[1]*qx.child.new.opt[1])
	mx2 <- qx.child.new.opt[2]/(4-axs[2]*qx.child.new.opt[2])

	mx.child.new.opt <- c(mx1, mx2)
	out.mort.new.opt <- unname(c(log(mx.child.new.opt), out.mort.start[25:44]))
	if(out.mort.new.opt[3]>out.mort.new.opt[2]){
		# y <- out.mort.new.opt[1:2]
		# x <- c(1,2)
		# child.mod.lm <- lm(y~x)
	   # new.5q5 <- predict.lm(child.mod.lm, newdata=data.frame(x = 3))
	   new.5m5 <- approx(x=c(2,4),y=out.mort.new.opt[c(2,4)],xout=3)$y
	   out.mort.new.opt[3] <- new.5m5
	}
	lt.out.new.opt <- lt.mx(nmx=exp(out.mort.new.opt), sex="female", age=c(0,1,seq(5,100,5)), nax=c(lt.out$lt[,2], rep(2.5, 22-length(lt.out$lt[,2]))))
	out.mort <- c(out.mort.start[1:22], out.mort.new.opt)
} # if

 
return(exp(out.mort[(22+1):(22*2)]))
}
}
