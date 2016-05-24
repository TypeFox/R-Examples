HILL <- function(par=c(S,H),G,P,type){
	# this function calculates the joint likelihood (-log(L)) for S and H given a genotype as matrix G 
	# if type == "codominant", G must be a two-row matrix with one column for each locus, as in STRUCTURE
	# if type == "dominant", G is a vector of 0,1 for absence,presence of the dominant allele
	# if type == "allele.count", G must be a vector of genotypes coded as 0,1,2 for the number of "j" alleles. That is, genotype 2 is homozygous for allele j, genotype 1 is heterozygous, and genotype 0 has no j alleles.
	# S and H are ancestry and heterozygosity indices (fractions on (0,1))
	# P contains rij and sij are estimates of the parental allele frequencies
	k <- dim(P)[1]
	S <- par[1]
	H <- par[2]
	if(H>min(2*S,2-2*S)+1e-8) {return(-1000*k*(H+abs(S-1/2)))}
	
	nonzer0 <- function(x){replace(x,x<1e-107,1e-107)}
	P[,3] <- nonzer0(P[,3])
	P[,4] <- nonzer0(P[,4])
	p12 <- max(H,1e-107)
	p11 <- max(S-p12/2,1e-107)
	p22 <- max(1-p11-p12,1e-107)

		LCD <- function(G,P,p11,p12,p22){ # assumes P is a data.frame or matrix with one row per locus, with P1 allele frequency in the third column, P2 allele frequency in the fourth column
		if(is.matrix(G)) G<- colSums(G)
		G <- G[!is.na(G)]
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
	LDD <- function(G,P,p11,p12,p22){
		if(is.matrix(G)) G<- colSums(G)
		#G <- G[!is.na(G)]
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
	LLck <- function(G,P,p11,p12,p22){
	# this function calculates the joint likelihood (-log(L)) for S and H given a matrix G of diploid genotypes in the format of STRUCTURE
	# S and H are ancestry and heterozygosity indices (fractions on (0,1))
	# G has two rows (for the two chromosomes of a diploid)
	# the columns of G are markers
	# P is a list of parental allele frequencies in a similar format to the estimated allele frequencies in STRUCTURE output
	# alleles should be represented by integers (arbitrary, but could be the number of repeats in a microsatellite allele)
	k <- ncol(G)
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
			}
	}}
		}
	sum(Lf,na.rm=TRUE)
	}
	if(type=="codominant"){
#		cat("data entered as CODOMINANT\nthe log-likelihood is\n")
		LL <- LLck(G,P,p11,p12,p22)
		}
	if(type=="dominant"){
#		cat("data entered as DOMINANT\nthe log-likelihood is\n")
		LL <- LDD(G,P,p11,p12,p22)
		}
	if(type=="allele.count"){
#		cat("data entered as ALLELE COUNTS\nthe log-likelihood is\n")
		LL <- LCD(G,P,p11,p12,p22)
		}
	LL
	}
