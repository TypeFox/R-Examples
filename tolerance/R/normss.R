norm.ss <- function(x = NULL, alpha = 0.05, P = 0.99, delta = NULL, P.prime = NULL, side = 1, m = 50, spec = c(NA, NA), hyper.par = list(mu.0 = NULL, sig2.0 = NULL, m.0 = NULL, n.0 = NULL),
	method = c("DIR", "FW", "YGZO")){
	method <- match.arg(method)
    if (side != 1 && side != 2) { 
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
	if (is.na(spec[1])) spec.L <- NULL else spec.L <- spec[1]
	if (is.na(spec[2])) spec.U <- NULL else spec.U <- spec[2]
	if (method=="DIR"){
    if ((is.null(hyper.par$mu.0)|is.null(hyper.par$sig2.0))&(is.null(spec.L)&is.null(spec.U))) { 
        stop(paste("Must specify mu.0 and sig2.0 as well as the appropriate spec limit(s)!", 
            "\n"))
    }
	mu.0 <- hyper.par$mu.0
	s.0 <- hyper.par$sig2.0
			f1 <- function(n,mu,sigma,alpha,P=P,side,spec.U) spec.U - (mu+K.factor(n=n,alpha=alpha,P=P,side=side,method="OCT",m=m)*sigma)
			f2 <- function(n,mu,sigma,alpha,P=P,side,spec.L) (mu-K.factor(n=n,alpha=alpha,P=P,side=side,method="OCT",m=m)*sigma) - spec.L
		if(side==1){
		if(is.null(spec.L)){
			n <- try(ceiling(uniroot(f1,mu=mu.0,sigma=s.0,alpha=alpha,P=P,side=side,spec.U=spec.U,interval=c(2,1e100))$root),silent=TRUE)
			f.calc <- f1(2,mu=mu.0,sigma=s.0,alpha=alpha,P=P,side=side,spec.U=spec.U)
		} else{
			n <- try(ceiling(uniroot(f2,mu=mu.0,sigma=s.0,alpha=alpha,P=P,side=side,spec.L=spec.L,interval=c(2,1e100))$root),silent=TRUE)
			f.calc <- f2(2,mu=mu.0,sigma=s.0,alpha=alpha,P=P,side=side,spec.L=spec.L)
		}
			if(class(n)=="try-error" & f.calc<0) n <- Inf
			if(class(n)=="try-error" & f.calc>=0) n <- 2
	} else{
		dL <- abs(mu.0-spec.L)
		dU <- abs(mu.0-spec.U)
		if(dL <= dU){
			n <- try(ceiling(uniroot(f2,mu=mu.0,sigma=s.0,alpha=alpha,P=P,side=1,spec.L=spec.L,interval=c(2,1e10))$root),silent=TRUE)
			f.calc <- f2(2,mu=mu.0,sigma=s.0,alpha=alpha,P=P,side=1,spec.L=spec.L)
		} else{
			n <- try(ceiling(uniroot(f1,mu=mu.0,sigma=s.0,alpha=alpha,P=P,side=1,spec.U=spec.U,interval=c(2,1e10))$root),silent=TRUE)
			f.calc <- f1(2,mu=mu.0,sigma=s.0,alpha=alpha,P=P,side=1,spec.U=spec.U)
		}
		if(class(n)=="try-error"){ 
			if(f.calc<0){
				n <- Inf
			} else n <- 2
		} else{
		TI.0 <- mu.0 + c(-1,1)*K.factor(n=1e100,alpha=alpha,P=P,side=2,method="HE")*s.0
		if((TI.0[1]<=spec.L) | (TI.0[2]>=spec.U)){
			n <- Inf
		} else{
		within.spec <- FALSE
		new.n <- n
		TI.Test <- function(x,L,U) (x[1]>=L)&(x[2]<=U)
		inc <- 1
		while((sum(within.spec)==0)&(n<Inf)){
			new.n <- new.n:(new.n+inc*1000)
			K2 <- K.factor(n=new.n,alpha=alpha,P=P,side=2,method="HE")
			TI <- cbind(mu.0 - K2*s.0, mu.0 + K2*s.0)
			within.spec <- apply(TI,1,TI.Test,L=spec.L,U=spec.U)
			if(sum(within.spec)==0) new.n <- tail(new.n,1)+1 else n <- new.n[min(which(within.spec))]
			inc <- inc+1
			if(inc>500) n <- Inf
		}
		within.spec <- FALSE
		n.old <- n 
		n <- max(1,n-8)
		brk <- TRUE
		while(!within.spec&brk){
			n <- n + 1
			TI <- try(mu.0 + c(-1,1)*K.factor(n=n,alpha=alpha,P=P,side=2,method="OCT",m=m)*s.0,silent=TRUE)
			if(class(TI)=="try-error"){
				n <- n.old
				brk <- FALSE
			} else within.spec <- (TI[1]>=spec.L) & (TI[2]<=spec.U)		
		}
		}
		}
	}
	new.n <- as.numeric(n)
	P.prime <- delta <- ""	
	} else{
	if (method == "YGZO"){
	    if (is.null(x)) { 
	        stop(paste("Data must be provided to use this method!", 
	            "\n"))
	    }
	    if ((is.null(spec.L)&is.null(spec.U))) { 
	        stop(paste("Must specify the appropriate spec limit(s) for this method!", 
	            "\n"))
	    }
	if(is.null(hyper.par$m.0)&is.null(hyper.par$n.0)){
		TI.out <- as.numeric(normtol.int(x=x, alpha=alpha, P=P, side=side, method="EXACT", m=m)[4:5])
	} else TI.out <- as.numeric(bayesnormtol.int(x=x, alpha=alpha, P=P, side=side, hyper.par=hyper.par)[3:4])
	mu.0 <- hyper.par$mu.0
	sig2.0 <- hyper.par$sig2.0
	s.0 <- sqrt(sig2.0)
	m.0 <- hyper.par$m.0
	n.0 <- hyper.par$n.0
	if (is.null(delta)|is.null(P.prime)){ 
			if((side==1)&length(c(spec.L,spec.U))!=1) stop(paste("You must specify a single value for one (and only one) of spec.L or spec.U.","\n"))
			if((side==2)&(is.null(spec.L)|is.null(spec.U))) stop(paste("Values for both spec.L and spec.U must be specified.","\n"))
			if (is.null(P.prime)){
				if(is.null(c(spec.L,spec.U))){ 
					P.prime <- (1+P)/2
				} else{
				if(side==2){
				P.prime <- pnorm(spec.U,mean=mu.0,sd=s.0)-pnorm(spec.L,mean=mu.0,sd=s.0) ### P.prime 2-sided criterion.
				if((P.prime<=P)|(P.prime>=1)) P.prime <- (1+P)/2
				} else{
					P.prime <- ifelse(!is.null(spec.L),pnorm(spec.L,mean=mu.0,sd=s.0,lower.tail=FALSE),pnorm(spec.U,mean=mu.0,sd=s.0)) ### P.prime 1-sided criterion.
					if((P.prime<=P)|(P.prime>=1)) P.prime <- (1+P)/2
				}
				}
			}
		if (is.null(delta)){
			if(side==1){
				if(!is.null(spec.L)){
					cont <- pnorm(TI.out[1],mean=mu.0,sd=s.0,lower.tail=FALSE)
					delta <- abs(cont-P)/P ### delta 1-sided lower criterion.
				} else{
					cont <- pnorm(TI.out[2],mean=mu.0,sd=s.0)
					delta <- abs(cont-P)/P ### delta 1-sided upper criterion.
				}
			} else{
				    if (is.null(spec.L)&is.null(spec.U)) { 
				        stop(paste("Must specify both spec limits!", 
				            "\n"))
				    }
				cont <- diff(pnorm(TI.out,mean=mu.0,sd=s.0))
				delta <- abs(cont-P)/P ### delta 2-sided criterion.
			}
		}
	} 
	}
	if ((method == "FW")&(is.null(delta)|is.null(P.prime))) stop(paste("You must specify delta and P.prime.","\n"))
		if (side == 1){
			norm.1 <- function(n,P,alpha,P.prime,delta) K.factor(n=n,P=P.prime,alpha=1-delta,side=1)-K.factor(n=n,P=P,alpha=alpha,side=1)
			new.n <- floor(uniroot(norm.1,interval=c(2,1e10),alpha=alpha,P=P,P.prime=P.prime,delta=delta)$root)
		} else{
			norm.2 <- function(n,P,alpha,P.prime,delta) K.factor(n=n,P=P,alpha=alpha,side=2,method="HE")-K.factor(n=n,P=P.prime,alpha=1-delta,side=2,method="HE")
			norm.2.EX <- function(n,P,alpha,P.prime,delta) K.factor(n=n,P=P,alpha=alpha,side=2,method="EXACT")-K.factor(n=n,P=P.prime,alpha=1-delta,side=2,method="EXACT")
			n.star <- ceiling(uniroot(norm.2,interval=c(2,1e10),alpha=alpha,P=P,P.prime=P.prime,delta=delta)$root)
			new.n <- n.star+c(-2:2)
			new.n <- new.n[new.n>3] 
			out <- try(cbind(new.n,K.factor(n=new.n,P=P,alpha=alpha,side=2,method="EXACT",m=m),K.factor(n=new.n,P=P.prime,alpha=1-delta,side=2,method="EXACT",m=m)),silent=TRUE)
			if(class(out)=="try-error"){ 
				new.n <- n.star
			} else{
			diff <- (out[,2]-out[,3])
			if(sum(diff<0)==0){
				new.n <- n.star+c(2:6)
				out <- cbind(new.n,K.factor(n=new.n,P=P,alpha=alpha,side=2,method="EXACT",m=m),K.factor(n=new.n,P=P.prime,alpha=1-delta,side=2,method="EXACT",m=m))
				diff2 <- (out[,2]-out[,3])
				if(sum(diff2<0)==0){
					new.n <- tail(new.n,1)
					tst <- 1
					while(tst>0){
						new.n <- new.n+20
						tst <- K.factor(n=new.n,P=P,alpha=alpha,side=2,method="EXACT",m=m)-K.factor(n=new.n,P=P.prime,alpha=1-delta,side=2,method="EXACT",m=m)
					}
					tst <- 1
					new.n <- new.n-20
					while(tst>0){
						new.n <- new.n+1
						tst <- K.factor(n=new.n,P=P,alpha=alpha,side=2,method="EXACT",m=m)-K.factor(n=new.n,P=P.prime,alpha=1-delta,side=2,method="EXACT",m=m)
					}
				} else new.n <- out[min(which(diff2<0)),1]
			} else if(min(diff<0)==1){
				new.n <- n.star+c((-6):(-2))
				out <- cbind(new.n,K.factor(n=new.n,P=P,alpha=alpha,side=2,method="EXACT",m=m),K.factor(n=new.n,P=P.prime,alpha=1-delta,side=2,method="EXACT",m=m))
				new.n <- out[min(which((out[,2]-out[,3])<0)),1]
			} else new.n <- out[min(which(diff<0)),1]	
		}	
		}
	}		
	temp <- data.frame(alpha=alpha, P=P, delta=delta, P.prime=P.prime, n=new.n)
	rownames(temp) <- NULL
	temp
}

