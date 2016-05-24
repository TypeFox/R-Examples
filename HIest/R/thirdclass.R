thirdclass <- function(G,P,type="codominant"){
		g <- G
	k <- dim(P)[1]
	# nonzer0 <- function(x){replace(x,x<1e-107,1e-107)}
	# P[,3] <- nonzer0(P[,3])
	# P[,4] <- nonzer0(P[,4])
	# P[,5] <- nonzer0(P[,5])

	if(type=="dominant"){
		N <- nrow(G)
		Nindex <- 1:N
		L.fn <- function(par,G,P){
			G <- G[!is.na(G)]
			k <- length(G)
			# S <- par[1]
			# H <- par[2]
			# if(H>min(2*S,2-2*S)+1e-8) {return(-1000*k*(H+abs(S-1/2)))}
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
			# S <- par[1]
			# H <- par[2]
			# if(H>min(2*S,2-2*S)+1e-8) {return(-1000*k*(H+abs(S-1/2)))}
			# p12 <- max(H,1e-107)
			# p11 <- max(S-p12/2,1e-107)
			# p22 <- max(1-p11-p12,1e-107)
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
	g.out <- matrix(nrow=N,ncol=15)
		for(i in 1:N){
		if(i/10 == round(i/10)) {cat("\n",i)}else{cat(" *")}
		g <- G[Nindex==i,]
		c100000 <- L.fn(c(1,0,0,0,0,0),g,P) # P1
		c010000 <- L.fn(c(0,1,0,0,0,0),g,P) # P2
		c001000 <- L.fn(c(0,0,1,0,0,0),g,P) # P3

		c000100 <- L.fn(c(0,0,0,1,0,0),g,P) # F1 1x2
		c110200 <- L.fn(c(.25,.25,0,.5,0,0),g,P) # F2 1x2
		c100100 <- L.fn(c(.5,0,0,.5,0,0),g,P) # BC 1xF12
		c010100 <- L.fn(c(0,.5,0,.5,0,0),g,P) # BC 2xF12

		c000010 <- L.fn(c(0,0,0,0,1,0),g,P) # F1 1x3
		c101020 <- L.fn(c(.25,0,.25,0,.5,0),g,P) # F2 1x3
		c100010 <- L.fn(c(.5,0,0,0,.5,0),g,P) # BC 1xF13
		c001010 <- L.fn(c(0,0,.5,0,.5,0),g,P) # BC 3xF13

		c000001 <- L.fn(c(0,0,0,0,0,1),g,P) # F1 2x3
		c011002 <- L.fn(c(0,.25,.25,0,0,.5),g,P) # F2 2x3
		c010001 <- L.fn(c(0,.5,0,0,0,.5),g,P) # BC 2xF23
		c001001 <- L.fn(c(0,0,.5,0,0,.5),g,P) # BC 3xF23

		g.out[i,] <- c(c100000, c010000, c001000, c000100, c110200, c100100,c010100,c000010,c101020 ,c100010,c001010,c000001,c011002,c010001,c001001)
		}
	Class <- apply(g.out,1,which.max)
	fun.fn <- function(x){x[rank(x,ties.method="first")==15]-x[rank(x,ties.method="first")==14]}
	LLD <- apply(g.out,1,fun.fn)
	
	colnames(g.out) <- c("c100000", "c010000", "c001000", "c000100", "c110200", "c100100","c010100","c000010","c101020" ,"c100010","c001010","c000001","c011002","c010001","c001001")
	data.frame(g.out,Best=c("c100000", "c010000", "c001000", "c000100", "c110200", "c100100","c010100","c000010","c101020" ,"c100010","c001010","c000001","c011002","c010001","c001001")[Class],LLD=unlist(LLD))

}
