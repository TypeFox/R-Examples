HIest <- function(G,P,type,method="SANN",iterations=1000,Cscale=NULL,surf=FALSE,startgrid=99,start=c(.5,.5),control=list(fnscale=-1,maxit=iterations)){
	g <- G
	k <- dim(P)[1]
	nonzer0 <- function(x){replace(x,x<.Machine$double.xmin,.Machine$double.xmin)}
	P[,3] <- nonzer0(P[,3])
	P[,4] <- nonzer0(P[,4])
	HIproposal <- function(par,G,P){
		ps <- c(par[1]-par[2]/2,par[2],(1-par[1])-par[2]/2)
		ps <- replace(ps,ps<.Machine$double.xmin,.Machine$double.xmin)
		pn <- as.numeric(rmultinom(1,Cscale,ps)/Cscale)
		c(S=pn[1]+pn[2]/2,H=pn[2])
	}
	if(is.null(Cscale)){Cscale <- 100}

	if(type=="allele.count"){
		N <- nrow(G)
		Nindex <- 1:N
		L.fn <- function(par,G,P){
			G <- G[!is.na(G)]
			S <- par[1]
			H <- par[2]
			if(H>min(2*S,2-2*S)+1e-8) {return(-1000*k*(H+abs(S-1/2)))}
			p12 <- max(H,1e-107)
			p11 <- max(S-p12/2,1e-107)
			p22 <- max(1-p11-p12,1e-107)
			k <- length(G)
			rij <- P[,3]
			sij <- P[,4]
			rik <- 1-rij
			sik <- 1-sij
			Lf <- rep(NA,k)
			for(i in 1:k){
				if(G[i]==2) Lf[i] <- log(p11*rij[i]^2+p12*rij[i]*sij[i]+p22*sij[i]^2)
				if(G[i]==1) Lf[i] <- log(p11*2*rij[i]*rik[i]+p12*(rij[i]*sik[i]+rik[i]*sij[i])+p22*2*sij[i]*sik[i])
				if(G[i]==0) Lf[i] <- log(p11*rik[i]^2+p12*rik[i]*sik[i]+p22*sik[i]^2)
				}
			sum(Lf,na.rm=TRUE)
			}
		}
	if(type=="dominant"){
		N <- nrow(G)
		Nindex <- 1:N
		L.fn <- function(par,G,P){
			S <- par[1]
			H <- par[2]
			if(H>min(2*S,2-2*S)+1e-8) {return(-1000*k*(H+abs(S-1/2)))}
			p12 <- max(H,1e-107)
			p11 <- max(S-p12/2,1e-107)
			p22 <- max(1-p11-p12,1e-107)
			k <- length(G)
			rij <- P[,3]
			sij <- P[,4]
			rik <- 1-rij
			sik <- 1-sij
			Lf <- rep(NA,k)
			for(i in 1:k){
				if(G[i]==1) Lf[i] <- log(p11*rij[i]^2+p12*rij[i]*sij[i]+p22*sij[i]^2+p11*2*rij[i]*rik[i]+p12*(rij[i]*sik[i]+rik[i]*sij[i])+p22*2*sij[i]*sik[i])
				if(G[i]==0) Lf[i] <- log(p11*rik[i]^2+p12*rik[i]*sik[i]+p22*sik[i]^2)
				}
			sum(Lf,na.rm=TRUE)
			}
		}
	if(type=="codominant") {
		N <- nrow(G)/2
		k <- ncol(G)
		Nindex <- rep(1:N,each=2) 
		Loci <- levels(as.factor(P[,1]))
		P.dat <- P
		P <- list()
		for(i in 1:k){
			Locus <- Loci[i]
			temp <- P.dat[P.dat[,1]==Locus,]
			rownames(temp) <- temp[,2]
			temp <- temp[,3:4]
			P[[i]] <- temp
			names(P)[i] <- Locus
			}
		L.fn <- function(par,G,P){
			S <- par[1]
			H <- par[2]
			if(H>min(2*S,2-2*S)+1e-8) {return(-1000*k*(H+abs(S-1/2)))}
			p12 <- max(H,1e-107)
			p11 <- max(S-p12/2,1e-107)
			p22 <- max(1-p11-p12,1e-107)
			k <- ncol(G)
			Lf <- rep(NA,k)
			for(i in 1:k){
				marker <- P[[i]]
				if(is.na(G[1,i])==FALSE){
				if(is.na(G[2,i])==FALSE){
				if(G[1,i]==G[2,i]) { # locus is homozygous
				allele <- G[1,i]
				rij <- marker[rownames(marker)==allele,1]
				sij <- marker[rownames(marker)==allele,2]
				Lf[i] <- log(p11*rij^2+p12*rij*sij+p22*sij^2)
				}
				if(G[1,i]!=G[2,i]){ #locus is heterozygous
				allele1 <- G[1,i]
				allele2 <- G[2,i]
				rij <- marker[rownames(marker)==allele1,1]
				sij <- marker[rownames(marker)==allele1,2]
				rik <- marker[rownames(marker)==allele2,1]
				sik <- marker[rownames(marker)==allele2,2]
				Lf[i] <- log(p11*2*rij*rik+p12*(rij*sik+rik*sij)+p22*2*sij*sik)
				}}}}
			sum(Lf,na.rm=TRUE)
			}
		}
	g.out <- matrix(nrow=N,ncol=3)
	cat("please wait: estimating S and H for",N,"individuals\n")
	for(n in 1:N){
	if(n/10 == round(n/10)) {cat("\n",n)}else{cat(" *")}
	G <- g[Nindex==n,]
	if(method=="L-BFGS-B"){
		S <- start[1]
		H <- start[2]
		est <- optim(par=c(S,H),fn=L.fn,G=G,P=P,method=method,lower=c(0,0),upper=c(1,1),control=control)
		output <- c(S=est$par[1],H=est$par[2],logL=est$value)
		par <- c(S=est$par[1],H=est$par[2])
		}
	if(surf | method=="surf"){
		S <- H <- seq(from=0,to=1,length.out= startgrid)
		surface <- matrix(NA,ncol= startgrid,nrow= startgrid)
		SH <- numeric()
		for(j in 1: startgrid){
			for(i in 1: startgrid){
				surface[i,j] <- L.fn(c(S[i],H[j]),G,P)
				SH <- rbind(SH,c(S[i],H[j]))
				}}
		maxL <- which.max(surface)
		par <- SH[maxL,]
		output <- c(S=par[1],H=par[2],logL=max(surface))
		}else{par <- start}
	if(method=="SANN"){
		est <- optim(par=par,fn=L.fn,gr=HIproposal,G=G,P=P,method="SANN",control=control)
		output <- c(S=est$par[1],H=est$par[2],logL=est$value)
		par <- c(S=est$par[1],H=est$par[2])
	}
	if(method=="mcmc"){
		chain = array(dim = c(iterations+1,2))
		LLik <- L.fn(par,G,P)
		chain[1,] = par
		for(i in 1:iterations){
			proposal <- HIproposal(chain[i,],G,P)
			newLL <- L.fn(proposal,G,P)
			probab <- exp(newLL-LLik[i])
			if(runif(1) < probab){
				chain[i+1,] <- proposal
				LLik[i+1] <- newLL
			}else{
				chain[i+1,] <- chain[i,]
				LLik[i+1] <- LLik[i]
			}
			}
		chain <- cbind(chain,LLik)
		output <- chain[which.max(chain[,3]),]
		}
	g.out[n,] <- output
	}
	cat("\n")
	colnames(g.out) <- c("S","H","logLik")
	as.data.frame(g.out)
	}
