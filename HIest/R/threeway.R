threeway <- function(G,P,type="codominant",surf=TRUE,method="SANN",iterations=500,start=rep(1/6,6),Cscale=NULL,props=NULL,control=list(fnscale=-1,maxit=iterations)){
	# method can be "mcmc", "surf", or "SANN"
	g <- G
	k <- dim(P)[1]
	if(is.null(Cscale)){Cscale <- 100}
	GR <- function(par,G,P){
		par <- replace(par,par<.Machine$double.xmin,.Machine$double.xmin)
		as.numeric(rmultinom(1,Cscale,par)/Cscale)
		}
	if(surf){
		if(is.null(props)){
			X <- 0:9
			X <- X/9
			props <- expand.grid(X,X,X,X,X,X)
			props <- as.matrix(props[rowSums(props)==1,])
			props <- rbind(props,start)
			}
		}

	if(type=="dominant"){
		N <- nrow(G)
		Nindex <- 1:N
		L.fn <- function(par,G,P){
			G <- G[!is.na(G)]
			k <- length(G)
			p11 <- par[1]
			p22 <- par[2]
			p33 <- par[3]
			p12 <- par[4]
			p13 <- par[5]
			p23 <- par[6]
			rij <- P[,3]
			sij <- P[,4]
			qij <- P[,5]
			rik <- 1-rij
			sik <- 1-sij
			qik <- 1-qik
			Lf <- rep(NA,k)
			for(i in 1:k){
				if(G[i]==1) Lf[i] <- log(p11*rij[i]^2+p12*rij[i]*sij[i]+p22*sij[i]^2+p11*2*rij[i]*rik[i]+p12*(rij[i]*sik[i]+rik[i]*sij[i])+p22*2*sij[i]*sik[i]+p33*qij[i]^2+p13*(rij[i]*qik[i]+rik[i]*qij[i])+p23*(sij[i]*qik[i]+sik[i]*qij[i]))
				if(G[i]==0) Lf[i] <- log(p11*rik[i]^2+p12*rik[i]*sik[i]+p22*sik[i]^2+p33*qik[i]^2+p13*rik[i]*qik[i]+p23*sik[i]*qik[i])
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
			temp <- temp[,3:5]
			P[[i]] <- temp
			names(P)[i] <- Locus
			}
		L.fn <- function(par,G,P){
			p11 <- par[1]
			p22 <- par[2]
			p33 <- par[3]
			p12 <- par[4]
			p13 <- par[5]
			p23 <- par[6]
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
				qij <- marker[rownames(marker)==allele,3]
				Lf[i] <- log(p11*rij^2+p22*sij^2+p33*qij^2+p12*rij*sij+p13*rij*qij+p23*sij*qij)
				}
				if(G[1,i]!=G[2,i]){ #locus is heterozygous
				allele1 <- G[1,i]
				allele2 <- G[2,i]
				rij <- marker[rownames(marker)==allele1,1]
				sij <- marker[rownames(marker)==allele1,2]
				qij <- marker[rownames(marker)==allele1,3]
				rik <- marker[rownames(marker)==allele2,1]
				sik <- marker[rownames(marker)==allele2,2]
				qik <- marker[rownames(marker)==allele2,3]
				Lf[i] <- log(p11*2*rij*rik+p22*2*sij*sik+p33*2*qij*qik+p12*(rij*sik+rik*sij)+p13*(rij*qik+rik*qij)+p23*(sij*qik+sik*qij))
				}}}}
			sum(Lf,na.rm=TRUE)
			}
		}
	g.out <- matrix(nrow=N,ncol=10)
	cat("please wait: estimating S and H for",N,"individuals\n")
	for(n in 1:N){
	if(n/10 == round(n/10)) {cat("\n",n)}else{cat(" *")}
	G <- g[Nindex==n,]
	if(surf | method=="surf"){
		pL <- numeric()
		for(i in 1:nrow(props)){
				pL[i] <- L.fn(props[i,],G,P)
				}
		maxL <- which.max(pL)
		par <- props[maxL,]
		output <- c(par,maxL=pL[maxL])
		}else{par <- start}

	if(method=="mcmc"){
		chain = array(dim = c(iterations+1,6))
		LLik <- L.fn(par,G,P)
		chain[1,] = par
		for(i in 1:iterations){
			proposal <- GR(chain[i,],G,P)# as.numeric(rmultinom(1,Cscale,chain[i,])/Cscale)
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
		output <- chain[which.max(chain[,7]),]
		}
	
	if(method=="SANN"){
		est <- optim(par,L.fn,gr=GR,G=G,P=P,method="SANN",control=control)
		output <- c(est$par,est$val) 
	}
	
	S1 <- output[1]+output[4]/2 + output[5]/2
	S2 <- output[2]+output[4]/2+output[6]/2
	S3 <- output[3]+output[5]/2+output[6]/2
	g.out[n,] <- c(output,S1=S1,S2=S2,S3=S3)
	}
	cat("\n")		
	colnames(g.out) <- c("p11","p22","p33","p12","p13","p23","logLik","S1","S2","S3")
	g.out <- g.out[,c(1:6,8:10,7)]
	as.data.frame(g.out)
	}
