tolerance <- function(rs, test = "upper", method = "ASY", pooled = NULL, symchoice = NULL, hpsim = NULL, h0sim = NULL, reduce = 1, ITER = 1000, exactL2 = TRUE, comment = TRUE){ 
	
	t1 <- Sys.time()
	
	if(class(rs)!="rrs") stop("'rs' must be of class 'rrs'")
	
	if(all(c("upper","lower","double")!=test)) stop("'test' must be one of 'upper', 'lower' or 'double'")
	
	if(reduce<=0) stop("'reduce' must be greater than zero or less than or equal to one")
	if(reduce>1) stop("'reduce' must be greater than zero or less than or equal to one")

	adaptive <- (length(rs$f$hypoH)>1)
	fvec <- as.vector(t(rs$f$Zm))
	gvec <- as.vector(t(rs$g$Zm))

	xr <- range(rs$f$X)
	yr <- range(rs$f$Y)
	gsize <- ceiling(reduce*length(rs$f$X))

	if(gsize < 10) warning("given reduction results in smaller p-value grid than 10x10")

	grx <- sort(rep(seq(xr[1],xr[2],length=gsize),gsize))
	gry <- rep(seq(yr[1],yr[2],length=gsize),gsize)
		
	if(reduce==1){
		corrGridSpec <- 1:(length(rs$f$X)^2)
	} else {
		corrGridSpec <- apply(as.matrix(data.frame(cbind(grx,gry))),1,getNearest,gridx=sort(rep(rs$f$X,length(rs$f$X))),gridy=rep(rs$f$Y,length(rs$f$Y)),WIN=rs$f$WIN,anypoint=T)
	}


	if(method=="ASY"){
		
		if(is.null(pooled)) stop("'pooled' density estimate required for asymptotics")
		
		if(class(pooled)!="bivden") stop("'pooled' must be an object of class 'bivden'")
	
		if(!all(c(length(rs$f$zVec)==length(pooled$zVec),length(rs$g$zVec)==length(pooled$zVec)))){
			stop("'pooled' appears to have been estimated using a different evaluation grid to that used in 'rs' - evaluation grids must be identical")
		}
		
		if(!all(c(identical_windows(rs$f$WIN,pooled$WIN),identical_windows(rs$g$WIN,pooled$WIN)))){
			stop("'pooled$WIN' does not appear to match study region window used in 'rs' - study regions must be identical")
		}
		
		if(length(pooled$hypoH)!=length(rs$f$hypoH)) stop("smoothing approach (fixed or adaptive) of 'pooled' must match approach used in 'rs'")
	
		edgef <- range(as.vector(rs$f$qhz),na.rm=T)[1]!=range(as.vector(rs$f$qhz),na.rm=T)[2]
		edgep <- range(as.vector(pooled$qhz),na.rm=T)[1]!=range(as.vector(pooled$qhz),na.rm=T)[2]
	
		if((edgef+edgep==1)) stop("edge-correction is inconsistent. all densities must be either edge-corrected or not.")
	
		datarange.list <- list(x=seq(xr[1],xr[2],length=gsize),y=seq(yr[1],yr[2],length=gsize))
		
		if(edgep){
			if(adaptive){
				if(comment) cat("\n--Adaptive-bandwidth asymptotics--\n")
				if(comment) cat("calculating integrals K2...\n")
			
				if(comment) cat("--f--\n")
				hypoQuan <- unique(quantile(sqrt(0.5*as.vector(t(rs$f$hypoH))^2)[corrGridSpec],(1:100)/100,na.rm=T))
				corrQuan <- apply(as.matrix(sqrt(0.5*as.vector(t(rs$f$hypoH))^2)[corrGridSpec]),1,idQuan,q=hypoQuan)
				qQuan <- apply(as.matrix(hypoQuan),1,run_ppp,data=rs$f$data,xy=datarange.list,WIN=rs$f$WIN,counts=rs$f$counts)
				fk2 <- rep(-1,gsize*gsize)
				fk2[is.na(sqrt(0.5*as.vector(t(rs$f$hypoH))^2)[corrGridSpec])] <- NA
				for(i in 1:length(hypoQuan)) fk2[which(corrQuan==i)] <- as.vector(qQuan[[i]]$edg$v)[which(corrQuan==i)]
				fk2 <- (1/(4*pi))*fk2 
				
				if(comment) cat("--g--\n")
				hypoQuan <- unique(quantile(sqrt(0.5*as.vector(t(rs$g$hypoH))^2)[corrGridSpec],(1:100)/100,na.rm=T))
				corrQuan <- apply(as.matrix(sqrt(0.5*as.vector(t(rs$g$hypoH))^2)[corrGridSpec]),1,idQuan,q=hypoQuan)
				qQuan <- apply(as.matrix(hypoQuan),1,run_ppp,data=rs$g$data,xy=datarange.list,WIN=rs$f$WIN,counts=rs$g$counts)
				gk2 <- rep(-1,gsize*gsize)
				gk2[is.na(sqrt(0.5*as.vector(t(rs$g$hypoH))^2)[corrGridSpec])] <- NA
				for(i in 1:length(hypoQuan)) gk2[which(corrQuan==i)] <- as.vector(qQuan[[i]]$edg$v)[which(corrQuan==i)]
				gk2 <- (1/(4*pi))*gk2 
				
				if(exactL2){
					xrL <- sort(rep(seq(-4,4,length=10),10))
					yrL <- rep(seq(-4,4,length=10),10)
					grL <- matrix(c(xrL,yrL),100,2)
					Lsq_gr <- apply(grL,1,Lsq,uh=c(0,0,1),WIN=NULL)
					
					
					if(comment) cat("calculating integrals L2...\n--f--\n")
					coords <- as.matrix(data.frame(cbind(grx,gry,as.vector(t(rs$f$hypoH))[corrGridSpec])))
					fL2 <- c()
					for(i in 1:nrow(coords)){
						if(!inside.owin(coords[i,1],coords[i,2],rs$f$WIN)||is.na(coords[i,3])){
							tempLsq <- NA
						} else {
							temp.xrL <- xrL*coords[i,3]+coords[i,1]
							temp.yrL <- yrL*coords[i,3]+coords[i,2]
							tempLsq1 <- Lsq_gr[inside.owin(temp.xrL,temp.yrL,rs$f$WIN)]
							temp.ygap <- temp.yrL[2]-temp.yrL[1]
							tempLsq <- sum((tempLsq1/(coords[i,3]^4))*temp.ygap^2)
						}
						fL2 <- append(fL2,tempLsq)
					}		
					
					S1rzK <- (1/(as.vector(t(rs$f$qhz))[corrGridSpec]^2))*(2*fk2 + 0.25*fL2)
				
					
					
					if(comment) cat("--g--\n\n")
					coords <- as.matrix(data.frame(cbind(grx,gry,as.vector(t(rs$g$hypoH))[corrGridSpec])))
					gL2 <- c()
					for(i in 1:nrow(coords)){
						if(!inside.owin(coords[i,1],coords[i,2],rs$f$WIN)||is.na(coords[i,3])){
							tempLsq <- NA
						} else {
							temp.xrL <- xrL*coords[i,3]+coords[i,1]
							temp.yrL <- yrL*coords[i,3]+coords[i,2]
							tempLsq <- Lsq_gr[inside.owin(temp.xrL,temp.yrL,rs$f$WIN)]
							temp.ygap <- temp.yrL[2]-temp.yrL[1]
							tempLsq <- sum((tempLsq/(coords[i,3]^4))*temp.ygap^2)
						}
						gL2 <- append(gL2,tempLsq)
					}
					S2rzK <- (1/(as.vector(t(rs$g$qhz))[corrGridSpec]^2))*(2*gk2 + 0.25*gL2)
				
					
					
					
				} else {
					if(comment) cat("calculating integrals L2...\n--f--\n")
					S1rzK <- (1/(as.vector(t(rs$f$qhz))[corrGridSpec]^2))*(2*fk2 + .5*fk2)
					if(comment) cat("--g--\n\n")
					S2rzK <- (1/(as.vector(t(rs$g$qhz))[corrGridSpec]^2))*(2*gk2 + .5*gk2)
				}
				
				denominator <- sqrt(((S1rzK*rs$f$gamma^2)/(sum(rs$f$counts)*rs$f$globalH^2))+((S2rzK*rs$g$gamma^2)/(sum(rs$g$counts)*rs$g$globalH^2)))
			} else {
				if(comment) cat("\n--Fixed-bandwidth asymptotics--\n")
				if(comment) cat("calculating integrals K2...")
			
				k2fix <- (1/(4*pi))*as.vector(run_ppp(data=pooled$data,xy=datarange.list,h=(sqrt(0.5*pooled$pilotH^2)),WIN=pooled$WIN,counts=pooled$counts)$edg$v)
				
				if(comment) cat("done.\n\n")
				h <- unique(pooled$h)
				RrzK <- k2fix/(as.vector(t(pooled$qhz))[corrGridSpec]^2)
				denominator <- sqrt(RrzK*(sum(rs$f$counts)^(-1)+sum(rs$g$counts)^(-1)))/(h*sqrt(as.vector(t(pooled$Zm)))[corrGridSpec])
			}
		} else {
			
			if(adaptive){
				if(comment) cat("\n--Adaptive-bandwidth asymptotics--\n")
				if(comment) cat("calculating integrals K and K2...\n")
			
				if(comment) cat("--f--\n")
				hypoQuan <- unique(quantile(as.vector(t(rs$f$hypoH))[corrGridSpec],(1:100)/100,na.rm=T))
				corrQuan <- apply(as.matrix(as.vector(t(rs$f$hypoH))[corrGridSpec]),1,idQuan,q=hypoQuan)
				qQuan <- apply(as.matrix(hypoQuan),1,run_ppp,data=rs$f$data,xy=datarange.list,WIN=rs$f$WIN,counts=rs$f$counts)
				fk <- rep(-1,gsize*gsize)
				fk[is.na(as.vector(t(rs$f$hypoH))[corrGridSpec])] <- NA
				for(i in 1:length(hypoQuan)) fk[which(corrQuan==i)] <- as.vector(qQuan[[i]]$edg$v)[which(corrQuan==i)]
				
				hypoQuan <- unique(quantile(sqrt(0.5*as.vector(t(rs$f$hypoH))^2)[corrGridSpec],(1:100)/100,na.rm=T))
				corrQuan <- apply(as.matrix(sqrt(0.5*as.vector(t(rs$f$hypoH))^2)[corrGridSpec]),1,idQuan,q=hypoQuan)
				qQuan <- apply(as.matrix(hypoQuan),1,run_ppp,data=rs$f$data,xy=datarange.list,WIN=rs$f$WIN,counts=rs$f$counts)
				fk2 <- rep(-1,gsize*gsize)
				fk2[is.na(sqrt(0.5*as.vector(t(rs$f$hypoH))^2)[corrGridSpec])] <- NA
				for(i in 1:length(hypoQuan)) fk2[which(corrQuan==i)] <- as.vector(qQuan[[i]]$edg$v)[which(corrQuan==i)]
				fk2 <- (1/(4*pi))*fk2 
				
				
				if(comment) cat("--g--\n")
				hypoQuan <- unique(quantile(as.vector(t(rs$g$hypoH))[corrGridSpec],(1:100)/100,na.rm=T))
				corrQuan <- apply(as.matrix(as.vector(t(rs$g$hypoH))[corrGridSpec]),1,idQuan,q=hypoQuan)
				qQuan <- apply(as.matrix(hypoQuan),1,run_ppp,data=rs$g$data,xy=datarange.list,WIN=rs$g$WIN,counts=rs$g$counts)
				gk <- rep(-1,gsize*gsize)
				gk[is.na(as.vector(t(rs$g$hypoH))[corrGridSpec])] <- NA
				for(i in 1:length(hypoQuan)) gk[which(corrQuan==i)] <- as.vector(qQuan[[i]]$edg$v)[which(corrQuan==i)]
				
				hypoQuan <- unique(quantile(sqrt(0.5*as.vector(t(rs$g$hypoH))^2)[corrGridSpec],(1:100)/100,na.rm=T))
				corrQuan <- apply(as.matrix(sqrt(0.5*as.vector(t(rs$g$hypoH))^2)[corrGridSpec]),1,idQuan,q=hypoQuan)
				qQuan <- apply(as.matrix(hypoQuan),1,run_ppp,data=rs$g$data,xy=datarange.list,WIN=rs$g$WIN,counts=rs$g$counts)
				gk2 <- rep(-1,gsize*gsize)
				gk2[is.na(sqrt(0.5*as.vector(t(rs$g$hypoH))^2)[corrGridSpec])] <- NA
				for(i in 1:length(hypoQuan)) gk2[which(corrQuan==i)] <- as.vector(qQuan[[i]]$edg$v)[which(corrQuan==i)]
				gk2 <- (1/(4*pi))*gk2 
				
				if(exactL2){
					xrL <- sort(rep(seq(-4,4,length=10),10))
					yrL <- rep(seq(-4,4,length=10),10)
					grL <- matrix(c(xrL,yrL),100,2)
					Lsq_gr <- apply(grL,1,Lsq,uh=c(0,0,1),WIN=NULL)
					
					if(comment) cat("calculating integrals L2...\n--f--\n")
					coords <- as.matrix(data.frame(cbind(grx,gry,as.vector(t(rs$f$hypoH))[corrGridSpec])))
					fL2 <- c()
					for(i in 1:nrow(coords)){
						if(!inside.owin(coords[i,1],coords[i,2],rs$f$WIN)||is.na(coords[i,3])){
							tempLsq <- NA
						} else {
							#h <<- coords[i,3]
							temp.xrL <- xrL*coords[i,3]+coords[i,1]
							temp.yrL <- yrL*coords[i,3]+coords[i,2]
							tempLsq1 <- Lsq_gr[inside.owin(temp.xrL,temp.yrL,rs$f$WIN)]
							temp.ygap <- temp.yrL[2]-temp.yrL[1]
							tempLsq <- sum((tempLsq1/(coords[i,3]^4))*temp.ygap^2)
						}
						fL2 <- append(fL2,tempLsq)
					}		
					S1rzK <- (1/(fk^2))*(2*fk2 + 0.25*fL2)
					
					
					if(comment) cat("--g--\n\n")
					coords <- as.matrix(data.frame(cbind(grx,gry,as.vector(t(rs$g$hypoH))[corrGridSpec])))
					gL2 <- c()
					for(i in 1:nrow(coords)){
						if(!inside.owin(coords[i,1],coords[i,2],rs$f$WIN)||is.na(coords[i,3])){
							tempLsq <- NA
						} else {
							temp.xrL <- xrL*coords[i,3]+coords[i,1]
							temp.yrL <- yrL*coords[i,3]+coords[i,2]
							tempLsq <- Lsq_gr[inside.owin(temp.xrL,temp.yrL,rs$f$WIN)]
							temp.ygap <- temp.yrL[2]-temp.yrL[1]
							tempLsq <- sum((tempLsq/(coords[i,3]^4))*temp.ygap^2)
						}
						gL2 <- append(gL2,tempLsq)
					}
					S2rzK <- (1/(gk^2))*(2*gk2 + 0.25*gL2)
				} else {
					if(comment) cat("calculating integrals L2...\n--f--\n")
					S1rzK <- (1/(fk^2))*(2*fk2 + .5*fk2)
					if(comment) cat("--g--\n\n")
					S2rzK <- (1/(gk^2))*(2*gk2 + .5*gk2)
				}
				
				denominator <- sqrt(((S1rzK*rs$f$gamma^2)/(sum(rs$f$counts)*rs$f$globalH^2))+((S2rzK*rs$g$gamma^2)/(sum(rs$g$counts)*rs$g$globalH^2)))
			} else {
				if(comment) cat("\n--Fixed-bandwidth asymptotics--\n")
				if(comment) cat("calculating integrals K and K2...")
			
				kfix <- as.vector(run_ppp(data=pooled$data,xy=datarange.list,h=pooled$pilotH,WIN=pooled$WIN,counts=pooled$counts)$edg$v)
				k2fix <- (1/(4*pi))*as.vector(run_ppp(data=pooled$data,xy=datarange.list,h=(sqrt(0.5*pooled$pilotH^2)),WIN=pooled$WIN,counts=pooled$counts)$edg$v)
				
				if(comment) cat("done.\n\n")
				
				RrzK <- k2fix/(kfix^2)
				denominator <- sqrt(RrzK*(unique(rs$f$h)^(-2)*sum(rs$f$counts)^(-1)+unique(rs$g$h)^(-2)*sum(rs$g$counts)^(-1)))/(sqrt(as.vector(t(pooled$Zm)))[corrGridSpec])
			}
		}
		
		if(rs$log) {
			numerator <- as.vector(t(rs$rsM))[corrGridSpec]
		} else {
			numerator <- as.vector(t(rs$rsM))[corrGridSpec]-1
		}
		Zstandard <- numerator/denominator
		
		if(test=="upper"){
			P <- pnorm(Zstandard,lower.tail=F)
		} else if (test=="lower"){
			P <- pnorm(Zstandard,lower.tail=T)
		} else {
			P <- 2*pnorm(abs(Zstandard),lower.tail=F)
		}
		if(comment){ 
			t2 <- Sys.time()
			cat("\n")
			print(t2-t1)
		}
		return(list(X=seq(xr[1],xr[2],length=gsize),Y=seq(yr[1],yr[2],length=gsize),Z=matrix(Zstandard,gsize,gsize,byrow=T),P=matrix(P,gsize,gsize,byrow=T)))
	
	} else if(method=="MC"){
		t1 <- Sys.time()
		
		ITER <- round(ITER)
		
		if(!adaptive){
			mcvals <- rsmc.fix(rs,hpsim,ITER,corrGridSpec,comment)
		} else if(!is.null(symchoice)){
			if(is.null(rs$f$pdef)||is.null(rs$g$pdef)) stop("Cannot perform symmetric adaptive simulations -- need identical 'pdef' objects present in rs$f and rs$g objects")
			mcvals <- rsmc.sym(rs,symchoice=symchoice,hpsim,h0sim,ITER,corrGridSpec,comment)
		} else {
			mcvals <- rsmc.asym(rs,hpsim,h0sim,ITER,corrGridSpec,comment)
		}
		
		p <- p.temp <- as.vector(t(mcvals))
		
		if(test=="lower") p <- 1-p
		if(test=="double"){
			p[which(p.temp<=0.5)] <- 2*p.temp[which(p.temp<=0.5)]
			p[which(p.temp>0.5)] <- 2-2*p.temp[which(p.temp>0.5)]
		}
		
		if(comment){ 
			t2 <- Sys.time()
			cat("\n")
			print(t2-t1)
		}	
		return(list(X=seq(xr[1],xr[2],length=gsize),Y=seq(yr[1],yr[2],length=gsize),Z=NA,P=matrix(p,gsize,gsize,byrow=T)))
	} else {
		stop("'method' must be one of 'ASY' or 'MC'")
	}
}

