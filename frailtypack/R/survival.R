survival <- function(t,ObjFrailty){

	if (ObjFrailty$typeof == 0){

		nz <- ObjFrailty$n.knots
		the <- ObjFrailty$b[1:(nz+2)] * ObjFrailty$b[1:(nz+2)]
		zi <- ObjFrailty$zi

		res <- NULL
		if(class(ObjFrailty) == "jointPenal"){
			nst <- ObjFrailty$n.strat + 1 # deces
			if((ObjFrailty$xR[,1] > t) || ((max(ObjFrailty$xR[,1])+0.00001) < t)) stop(" Time exceeds the range allowed ")
			if(ObjFrailty$n.strat > 1){
				for (i in 2:ObjFrailty$n.strat){
					if((ObjFrailty$xR[,i] > t) || (max(ObjFrailty$xR[,i]+0.00001) < t)) stop(" Time exceeds the range allowed ")
					b <- ObjFrailty$b[((i-1)*(nz+2)+1):(i*(nz+2))]
					the <- cbind(the,b*b)
				}
			}
			if((ObjFrailty$xD > t) || (max(ObjFrailty$xD) < t)) stop(" Time exceeds the range allowed ")
			b <- ObjFrailty$b[((nst-1)*(nz+2)+1):(nst*(nz+2))]
			the <- cbind(the,b*b)
		}else{
			nst <- ObjFrailty$n.strat
			if((ObjFrailty$x[,1] > t) || ((max(ObjFrailty$x[,1])+0.00001) < t)) stop(" Time exceeds the range allowed ")
			if(ObjFrailty$n.strat > 1){
				for (i in 2:ObjFrailty$n.strat){
					if((ObjFrailty$x[,i] > t) || (max(ObjFrailty$x[,i]) < t)) stop(" Time exceeds the range allowed ")
					b <- ObjFrailty$b[((i-1)*(nz+2)+1):(i*(nz+2))]
					the <- cbind(the,b*b)
				}
			}
		}
		out <- .Fortran("survival2",as.double(t),as.double(the),as.integer(nz+2),
		as.double(zi),survival=as.double(rep(0,nst)),as.integer(nst),PACKAGE = "frailtypack")

		res <- c(res,out$survival)
		return(res)
	}

	if (ObjFrailty$typeof == 1){
		res <- NULL
		nst <- ObjFrailty$n.strat
		if(class(ObjFrailty) == "jointPenal"){
			m <- nst*ObjFrailty$nbintervR + ObjFrailty$nbintervDC
			b <- ObjFrailty$b[1:m]
			time <- ObjFrailty$time
			timedc <- ObjFrailty$timedc
			if((ObjFrailty$xR[,1] > t) || (max(ObjFrailty$xSuR[,1]) < t)) stop(" Time exceeds the range allowed ")
			if((ObjFrailty$xD > t) || (max(ObjFrailty$xSuD) < t)) stop(" Time exceeds the range allowed ")
			out <- .Fortran("survivalj_cpm2",as.double(t),as.double(b),as.integer(nst+1),as.integer(ObjFrailty$nbintervR),
			as.integer(ObjFrailty$nbintervDC),as.double(time),as.double(timedc),
			survival=as.double(rep(0,nst+1)),PACKAGE = "frailtypack")

			res <- c(res,out$survival)
		}else{
			m <- ObjFrailty$n.strat*ObjFrailty$nbintervR
			b <- ObjFrailty$b[1:m]
			time <- ObjFrailty$time
			if((ObjFrailty$x[,1] > t) || ((max(ObjFrailty$xSu[,1])+0.00001) < t)) stop(" Time exceeds the range allowed ")
			if((ObjFrailty$n.strat == 2) & ((ObjFrailty$x[,2] > t) || (max(ObjFrailty$xSu[,2]) < t))) stop(" Time exceeds the range allowed ")
			out <- .Fortran("survival_cpm2",as.double(t),as.double(b),as.integer(nst),as.integer(ObjFrailty$nbintervR),
			as.double(time),survival=as.double(rep(0,nst)),PACKAGE = "frailtypack")

			res <- c(res,out$survival)
			
		}
		return(res)
	}


	if (ObjFrailty$typeof == 2){
		if(!t)stop(" Use only for time greater than 0")
		res <- NULL
		sh1 <- ObjFrailty$shape.weib[1]
		sc1 <- ObjFrailty$scale.weib[1]
		res <- c(res,exp(-(t/sc1)^sh1))
		if(class(ObjFrailty) == "jointPenal"){
			if(ObjFrailty$n.strat > 1){
				for (i in 2:ObjFrailty$n.strat){
					if((ObjFrailty$xR[,i] > t) || (max(ObjFrailty$xSuR[,i]) < t)) stop(" Time exceeds the range allowed ")
					sh1 <- ObjFrailty$shape.weib[i]
					sc1 <- ObjFrailty$scale.weib[i]
					res <- c(res,exp(-(t/sc1)^sh1))
				}
			}
		if((ObjFrailty$xD > t) || (max(ObjFrailty$xSuD) < t)) stop(" Time exceeds the range allowed ")
		sh1 <- ObjFrailty$shape.weib[ObjFrailty$n.strat+1]
		sc1 <- ObjFrailty$scale.weib[ObjFrailty$n.strat+1]
		res <- c(res,exp(-(t/sc1)^sh1))
		}else{
			if(ObjFrailty$n.strat > 1){
				for (i in 2:ObjFrailty$n.strat){
					if((ObjFrailty$x[,i] > t) || (max(ObjFrailty$xSu[,i]) < t)) stop(" Time exceeds the range allowed ")
					sh1 <- ObjFrailty$shape.weib[i]
					sc1 <- ObjFrailty$scale.weib[i]
					res <- c(res,exp(-(t/sc1)^sh1))
				}
			}
		}
		return(res)
	}

}
