kmeansvar <-
function(X.quanti=NULL,X.quali=NULL,init,iter.max=150,nstart=1,matsim=FALSE)
{
#init:  Either the number of clusters or a vector with group memberships
#iter.max: The maximum number of iterations allowed.
#nstart: If init is a number, how many random sets should be chosen?

	cl <- match.call()
	rec <- recod(X.quanti,X.quali)

	n <- rec$n
	p <- rec$p		#total number of variables
	p1 <- rec$p1	#number of numerical variables
	X <- rec$X 		#data matrix (mixed if necessary)
	Xn <- rec$Xn	#data matrix without any replication
	Z <- rec$Z		#data matrix for input in SVD
	indexj <- rec$indexj #variables indices in columns of Z
	G <- rec$G		#dummy variables if X.quanti non null
	Gcod <- rec$Gcod	
	X.quanti<-rec$X.quanti
	X.quali<-rec$X.quali
	
	sqcc <- function(X1,X2) {
		n <- nrow(X1)
		r <- ncol(X1)
		s <- ncol(X2)
		m <- which.min(c(n,r,s))
		if (m==1) {
			A1 <- X1%*%t(X1)/n
			A2 <- X2%*%t(X2)/n
			A <- A1%*%A2
			e <- eigen(A)
			sim <- Re(e$values[2]) 
		} else {
			V12 <- t(X1)%*%X2/n
			V21 <- t(X2)%*%X1/n
			if (m==2) V<-V12%*%V21
			if (m==3) V<-V21%*%V12
			e <- eigen(V)
			sim <- Re(e$values[2]) 
			}
		return(sim)
	}	

	eta2 <- function(x, gpe) {
		moyennes <- tapply(x, gpe, mean)
		effectifs <- tapply(x, gpe, length)
		varinter <- (sum(effectifs * (moyennes - mean(x))^2))
		vartot <- (var(x) * (length(x) - 1))
		res <- varinter/vartot
		return(res)
	}


	partinit <- function(centers)
	{
		A <- matrix(,p,k)
		#only quantitative data
		if (p1==p) {
			A <- (t(Z)%*%Z[,centers]/n)^2
		 	part <- apply(A,1,which.max)
		}
		#only qualitative data
		if (p1==0) {
			for (g in 1:k) {
				X1 <- Gcod[,which(indexj==centers[g])]
				for (j in 1:p) {
					X2 <- Gcod[,which(indexj==j)]
					A[j,g] <- sqcc(X1,X2)
				}
			}
 			part <- apply(A,1,which.max)
		}
		#a mixture of both
		if (p1 %in% 1:(p-1)) {
			nc1 <- length(centers[which(centers<=p1)])
			nc2 <- k-nc1
			if (nc1 != 0) {
				A[1:p1,1:nc1] <- t(Z[,1:p1])%*%Z[,centers[1:nc1]]/n
				for (g in 1:nc1) {
					for (j in (p1+1):p) A[j,g] <- eta2(X.quanti[,centers[g]],X.quali[,(j-p1)])
				}		
				}
			if (nc2 != 0) {
				for (g in (nc1+1):k) {
					for (j in 1:p1) A[j,g] <- eta2(X.quanti[,j], X.quali[,centers[g]-p1])
					X1 <- Gcod[,which(indexj==centers[g])-p1]
					for (j in (p1+1):p)  {
						X2 <- Gcod[,which(indexj==j)-p1]
						A[j,g] <- sqcc(X1,X2)
					}
					}	
			}
			part <- apply(A,1,which.max)
	}
	return(part)
	}

	if (missing(init)) 
        stop("'init' must be a number or a vector")
	if (length(init) == 1) {
		k <- init 
    centers <- sort(sample.int(p, k))
		part <- as.factor(partinit(centers))
		} else {
			if (!is.integer(init))
				stop("init must be a vector of integer")
			if (nstart!=1)
				stop("nstart must be equal to one")
        		part <- as.factor(init)
       		k <- length(levels(part))
	  		if (p < k) 
            		stop("more cluster centers than variables")
	  		if (length(which(init>k))> 0)
				stop("clusters must be numbered from 1 to k")
	  		if (p != length(part)) 
        			stop("the length of init must be equal to the number of variables")
            }

	do_one <- function(part)
	{
		A <- matrix(,p,k)
		iter <- 0
		diff <- TRUE
		while ((iter< iter.max) && diff)
		{
			iter <- iter+1

			#Representation step
			indexk<-NULL
			for (i in 1:length(indexj))
				indexk[i] <- part[indexj[i]]
			latent.var <- matrix(,n,k)
			sv<-rep(NA,k)
			for (g in 1:k) {
				Zclass <- Z[,which(indexk==g)]
				latent <- clusterscore(Zclass)
				latent.var[,g] <- latent$f
				sv[g] <- latent$sv }
		#Affectation step
			if (p1>0) {
					scorestand <- sweep(latent.var,2,STATS=sv,FUN="/")
					A[1:p1,]<-	(t(Z[,1:p1])%*%scorestand/n)^2
				}
			if (p1!=p)	{
				for (g in 1:k){
					sl.qual <- function(col) {
						eta2(latent.var[,g],col)
						}
					A[(p1+1):p,g]<-apply(X.quali,2,sl.qual)
					}	
			}
			part2 <- apply(A,1,which.max)

			diff <- !all(part==part2)
			part <- part2
		}
		wss <- sv^2 #within cluster sum of squares
		return(list(latent.var=latent.var,part=part,wss=wss,iter=iter))
	}

	res<-do_one(part)
	if (nstart >= 2) {
		best <- sum(res$wss)
        	for (i in 2:nstart) {
            	centers <- sort(sample.int(p, k))
			part<-as.factor(partinit(centers))
            	res2 <- do_one(part)
            	if ((z <- sum(res2$wss)) > best) {
               		res <- res2
                		best <- z}
        	}
    	}

	part <- res$part
	iter <- res$iter
	names(part) <- colnames(X)
	res <- descript(part,rec,cl,matsim,iter=iter)
	class(res) <- "clustvar"
	return(res)
}