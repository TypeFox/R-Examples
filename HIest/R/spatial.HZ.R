spatial.HZ <- function(minX ,minY, maxX, maxY, XY, Genotypes, 
						beta=0,sel=0, mid=0,  
						h=0, DM = matrix(0,ncol=3,nrow=3),
						sigmad, sigmac, sigmam, R, M, gens, immigrants=FALSE,plotgrowth=FALSE,m=0.10){
	# minX, etc define the space
	# XY is a matrix of initial X and Y coordinates of individuals
	# Genotypes is a matrix of allele counts for individuals (rows) by loci (columns)
	# beta, sel, and mid describe an environmental gradient affecting the first locus
	# h is heterozygote disadvantage affecting locus 2
	# DM describes the deleterious effects of interactions between loci 3 and 4
	# sigmad is the standard deviation of mother-offspring distance
	# sigmac is the sd of the competition function
	# sigmam is the standard deviation of the mother-father distance
	# R and M are the growth and capacity parameters affecting population regulation
	# gens: how many generations to run
	# immigrants: True if parentals should enter the pop whenever a hybrid leaves
	
	rbnorm <- function(n,mean=0,sd=1,min,max){ # random normal with reflecting boundaries
		X <- rnorm(n,mean,sd)
		X <- replace(X,X>max,max)
		X <- replace(X,X<min,min)
		X
		}
		
	if(plotgrowth){N<-numeric()}
	for(j in 1:gens){
		if(plotgrowth){N[j] <- dim(XY)[1]}
		EL <- Genotypes[,1]
		HL <- Genotypes[,2]
		DL <- Genotypes[,3:4]

		Dist <- as.matrix(dist(XY))

		# competition density: for each individual, sum across all others
		Dens <- dnorm(Dist,sd=sigmac)
		Dens <- rowSums(Dens)
	
		# mating density
		Mdens <- dnorm(Dist,sd=sigmam)
		diag(Mdens) <- 0
		
		# # environmental fitness component
		# Wn          <- West[1]+(East[1]-West[1])/(1+exp(beta[1]*(XY[,1]-mid)))
		# Wn[EL==0.5] <- West[1]+(East[1]-West[1])/(1+exp(beta[2]*(XY[EL==0.5,1]-mid)))
		# Wn[EL==0]   <- West[1]+(East[1]-West[1])/(1+exp(beta[3]*(XY[EL==0,1]-mid)))
		Grad <- exp(beta*(XY[,1]-mid))/(1+exp(beta*(XY[,1]-mid)))
		Wn <- exp(-sel*(EL-Grad)^2)
		
		# underdominance component
		Wh <- 1-h*(HL==0.5)
		
		# Dobzansky-Muller component
		Wdm <- 1-(DL[,1]==0)*((DL[,2]==0)*DM[1,1]+(DL[,2]==0.5)*DM[1,2]+(DL[,2]==1)*DM[1,3])-(DL[,1]==0.5)*((DL[,2]==0)*DM[2,1]+(DL[,2]==0.5)*DM[2,2]+(DL[,2]==1)*DM[2,3])-(DL[,1]==1)*((DL[,2]==0)*DM[3,1]+(DL[,2]==0.5)*DM[3,2]+(DL[,2]==1)*DM[3,3])
		
		# reproduction
		lambda <- R*(Wn*Wh*Wdm)/(1+Dens/M)
		lambda <- replace(lambda,lambda<0,0)
		nkids <- rpois(dim(XY)[1],lambda)
		moms <- (1:(dim(XY)[1]))[nkids>0]
		
		X <- Y <- numeric()
		newG <- matrix(NA,nrow=1,ncol=dim(Genotypes)[2])
		for(i in moms){
			X <- c(X,rbnorm(nkids[i],XY[i,1],sigmad,minX,maxX))
			Y <- c(Y,rbnorm(nkids[i],XY[i,2],sigmad,minY,maxY))
			NG <- matrix(NA,nrow=nkids[i],ncol=dim(Genotypes)[2])
			dad <- sample(dim(Genotypes)[1],1,prob=Mdens[,i])
			for(l in 1:(dim(Genotypes)[2])){
				XX <- rbinom(nkids[i],1,Genotypes[i,l])
				PD <- Genotypes[dad,l]#sum(lambda*Mdens[,i]*Genotypes[,l]/sum(lambda*Mdens[,i]))
				#if(PD>1){PD <- 1}
				YY <- rbinom(nkids[i],1,PD)
				NG[,l] <- (XX+YY)/2
			}
			newG <- rbind(newG,NG)
		}
		newG <- newG[-1,]
		if(immigrants==TRUE){
			newG[X>=maxX*(1-m/2),] <- rep(1,dim(Genotypes)[2])
			newG[X<=minX*(1-m/2),] <- rep(0,dim(Genotypes)[2])
		}
		XY <- cbind(X,Y)
		mothers <- Genotypes[moms,]
		Genotypes <- newG
	}
	if(plotgrowth){quartz();plot(1:gens,N)}
	list(XY=XY,Genotypes=Genotypes,mothers=mothers)
}
