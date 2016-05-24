"hazard" <- function(t,ObjFrailty){

	if (ObjFrailty$typeof == 0){

		nz <- ObjFrailty$n.knots
		the <- ObjFrailty$b[1:(nz+2)] * ObjFrailty$b[1:(nz+2)]
		zi <- ObjFrailty$zi

		res <- NULL
		if(class(ObjFrailty) == "jointPenal"){
			nst <- ObjFrailty$n.strat + 1 # deces
			if((ObjFrailty$xR[,1] > t) || ((max(ObjFrailty$xR[,1])+0.00001) <= t)) stop(" Time exceeds the range allowed ")
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
			if((ObjFrailty$x[,1] > t) || ((max(ObjFrailty$x[,1])+0.00001) <= t)) stop(" Time exceeds the range allowed ")
			if(ObjFrailty$n.strat > 1){
				for (i in 2:ObjFrailty$n.strat){
					if((ObjFrailty$x[,i] > t) || (max(ObjFrailty$x[,i]) < t)) stop(" Time exceeds the range allowed ")
					b <- ObjFrailty$b[((i-1)*(nz+2)+1):(i*(nz+2))]
					the <- cbind(the,b*b)
				}
			}
		}
		
		out <- .Fortran("risque2",as.double(t),as.double(the),as.integer(nz+2),
		as.double(zi),risque=as.double(rep(0,nst)),as.integer(nst),PACKAGE = "frailtypack")
		
		if(class(ObjFrailty) == "jointPenal"){
			res <- c(res,out$risque)
		}else{
			res <- out$risque
		}
		return(res)
	}

	if (ObjFrailty$typeof == 1){
		res <- NULL
		if(class(ObjFrailty) == "jointPenal"){
			if((ObjFrailty$xR[,1] > t) || (max(ObjFrailty$xR[,1]) < t)) stop(" Time exceeds the range allowed ")
			x1 <- matrix(ObjFrailty$xR[,1],nrow=3,ncol=ObjFrailty$nbintervR)
			x1 <- rbind(x1,rep(t,ObjFrailty$nbintervR))
			ind <- apply(x1,MARGIN=2, FUN=function(x){which((x[1] <= x[4]) & (x[4] < x[3]))})
			res <- c(res,ObjFrailty$lamR[(which(ind==1)*3-1),1,1])
			if(ObjFrailty$n.strat > 1){
				for (i in 2:ObjFrailty$n.strat){
					if((ObjFrailty$xR[,i] > t) || (max(ObjFrailty$xR[,i]) < t)) stop(" Time exceeds the range allowed ")
					x1 <- matrix(ObjFrailty$xR[,i],nrow=3,ncol=ObjFrailty$nbintervR)
					x1 <- rbind(x1,rep(t,ObjFrailty$nbintervR))
					ind <- apply(x1,MARGIN=2, FUN=function(x){which((x[1] <= x[4]) & (x[4] < x[3]))})
					res <- c(res,ObjFrailty$lamR[(which(ind==1)*3-1),1,i])
				}
			}
			if((ObjFrailty$xD > t) || (max(ObjFrailty$xD) < t)) stop(" Time exceeds the range allowed ")
			x2 <- matrix(ObjFrailty$xD,nrow=3,ncol=ObjFrailty$nbintervDC)
			x2 <- rbind(x2,rep(t,ObjFrailty$nbintervDC))
			ind <- apply(x2,MARGIN=2, FUN=function(x){which((x[1] <= x[4]) & (x[4] < x[3]))})
			res <- c(res,ObjFrailty$lamD[(which(ind==1)*3-1),1])
		}else{
# shared additive nested
			if((ObjFrailty$x[,1] > t) || (max(ObjFrailty$x[,1]) < t)) stop(" Time exceeds the range allowed ")
			x1 <- matrix(ObjFrailty$x[,1],nrow=3,ncol=ObjFrailty$nbintervR)
			x1 <- rbind(x1,rep(t,ObjFrailty$nbintervR))
			ind <- apply(x1,MARGIN=2, FUN=function(x){which((x[1] <= x[4]) & (x[4] < x[3]))})
			res <- c(res,ObjFrailty$lam[(which(ind==1)*3-1),1,1])
			if(ObjFrailty$n.strat > 1){
				for (i in 2:ObjFrailty$n.strat){
					if((ObjFrailty$x[,i] > t) || (max(ObjFrailty$x[,i]) < t)) stop(" Time exceeds the range allowed ")
					x1 <- matrix(ObjFrailty$x[,i],nrow=3,ncol=ObjFrailty$nbintervR)
					x1 <- rbind(x1,rep(t,ObjFrailty$nbintervR))
					ind <- apply(x1,MARGIN=2, FUN=function(x){which((x[1] <= x[4]) & (x[4] < x[3]))})
					res <- c(res,ObjFrailty$lam[(which(ind==1)*3-1),1,i])
				}
			}
		}
		return(res)
	}


	if (ObjFrailty$typeof == 2){
		if(!t)stop(" Use only for time greater than 0")
		res <- NULL
		if(class(ObjFrailty) == "jointPenal"){
			nst <- ObjFrailty$n.strat + 1 # deces
			if((ObjFrailty$xR[,1] > t) || (max(ObjFrailty$xR[,1]) < t)) stop(" Time exceeds the range allowed ")
			sc1 <- ObjFrailty$scale.weib[1]
			sh1 <- ObjFrailty$shape.weib[1]
			res <- c(res,(sh1*(t^(sh1-1)))/(sc1^sh1))
			if(ObjFrailty$n.strat > 1){
				for (i in 2:ObjFrailty$n.strat){
					if((ObjFrailty$xR[,i] > t) || (max(ObjFrailty$xR[,i]) < t)) stop(" Time exceeds the range allowed ")
					sc1 <- ObjFrailty$scale.weib[i]
					sh1 <- ObjFrailty$shape.weib[i]
					res <- c(res,(sh1*(t^(sh1-1)))/(sc1^sh1))
				}
			}
			if((ObjFrailty$xD > t) || (max(ObjFrailty$xD) < t)) stop(" Time exceeds the range allowed ")
			sc1 <- ObjFrailty$scale.weib[nst]
			sh1 <- ObjFrailty$shape.weib[nst]
			res <- c(res,(sh1*(t^(sh1-1)))/(sc1^sh1))
		}else{
			if((ObjFrailty$x[,1] > t) || (max(ObjFrailty$x[,1]) < t)) stop(" Time exceeds the range allowed ")
			sc1 <- ObjFrailty$scale.weib[1]
			sh1 <- ObjFrailty$shape.weib[1]
			res <- c(res,(sh1*(t^(sh1-1)))/(sc1^sh1))
			if(ObjFrailty$n.strat > 1){
				for (i in 2:ObjFrailty$n.strat){
					if((ObjFrailty$x[,i] > t) || (max(ObjFrailty$x[,i]) < t)) stop(" Time exceeds the range allowed ")
					sc1 <- ObjFrailty$scale.weib[i]
					sh1 <- ObjFrailty$shape.weib[i]
					res <- c(res,(sh1*(t^(sh1-1)))/(sc1^sh1))
				}
			}
		}
		return(res)
	}
}


