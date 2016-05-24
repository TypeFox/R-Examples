HIC3 <-
function(G,P){
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
		g.out <- matrix(nrow=N,ncol=10)
		for(n in 1:N){
			g <- G[Nindex==n,]
			GP <- matrix(0,ncol=k,nrow=6)
			for(i in 1:k){
				marker <- P[[i]]
				if(is.na(g[1,i])==FALSE){
				if(is.na(g[2,i])==FALSE){
				if(g[1,i]==g[2,i]) { # locus is homozygous
				allele <- g[1,i]
				rij <- marker[rownames(marker)==allele,1]
				sij <- marker[rownames(marker)==allele,2]
				qij <- marker[rownames(marker)==allele,3]
				GP[1:3,i] <- c(rij,sij,qij)
				}
				if(g[1,i]!=g[2,i]){ #locus is heterozygous
				allele1 <- g[1,i]
				allele2 <- g[2,i]
				rij <- marker[rownames(marker)==allele1,1]
				sij <- marker[rownames(marker)==allele1,2]
				qij <- marker[rownames(marker)==allele1,3]
				rik <- marker[rownames(marker)==allele2,1]
				sik <- marker[rownames(marker)==allele2,2]
				qik <- marker[rownames(marker)==allele2,3]
				GP[4:6,i] <- c(rij*sik+rik*sij,rij*qik+rik*qij,sij*qik+sik*qij)
				}}}}
			par <- rowMeans(GP)
			g.out[n,1:6] <- par
			g.out[n,7] <- par[1]+par[4]/2+par[5]/2
			g.out[n,8] <- par[2]+par[4]/2+par[6]/2
			g.out[n,9] <- par[3]+par[5]/2+par[6]/2
			g.out[n,10] <- L.fn(par,g,P)
			}
		colnames(g.out) <- c("p11","p22","p33","p12","p13","p23","S1","S2","S3","logLik")
	as.data.frame(g.out)
}
