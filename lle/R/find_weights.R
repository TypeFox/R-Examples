find_weights <-
function(nns,X,m,reg=2,ss=FALSE,p=0.5,id=FALSE,v=0.99)
{
	#get dimensions
	N <- dim(X)[1]
	n <- dim(X)[2]
	
	#matrix of weights
	wgts <- 0*matrix(0,N,N)
	#subset selection vector
	s <- c()
	#intrinsic dim
	intr_dim <- c()

	for (i in (1:N)){
		#number of neighbours
		k <- sum(nns[i,])
		
		#no neighbours (find_nn_k(k=0) or eps-neighbourhood)
		if( k==0 ) next
		
		# calculate the differences  between xi and its neighbours
		Z <- matrix(c(t(X)) - c(t(X[i,])), nrow=nrow(X), byrow = TRUE)
		Z <- matrix( Z[nns[i,],], ncol=n, nrow=k )		
		
		#gram-matrix
		G <- Z%*%t(Z)
		
		#regularisation
		delta <- 0.1
		#calculate eigenvalues of G
		e <- eigen(G, symmetric=TRUE, only.values=TRUE)$values
		#skip if all EV are null
		if( all( e==0 ) ) next
		
		#choose regularisation method
		#see documentation
		if( reg==1 ){ 
			r <- delta*sum(head( e, n-m ))/(n-m)
		} else if( reg==2 ){
			r <- delta^2/k*sum(diag(G))
		} else r <- 3*10^-3
		
		#calculate intrinsic dimension
		if( id ){	
			tmp <- 1
			while( sum(e[1:tmp])/sum(e) <= v ) tmp <- tmp + 1
			intr_dim <- c( intr_dim, tmp )
		}
		
		#basic value for subset selection
		s <- c( s, sum(head( e, n-m ))/(n-m) )
		
		#use regularisation if more neighbourse than dimensions!
		if( k>n ) alpha <- r else alpha <- 0
			
		#regularisation
		G <- G + alpha*diag(1,k)
		
		#calculate weights
		#using pseudoinverse ginv(A): works better for bad conditioned systems
		if( k >= 2) wgts[i,nns[i,]] <- t(ginv(G)%*%rep(1,k)) else wgts[i] <- G
		wgts[i,] <- wgts[i,]/sum(wgts[i,])
	}
	
	#subset selection
	#see documentation
	if( ss ){
		s <- s/sum(s)
		Fs <- ecdf(s)
		choise <- sample( 1:N, round(p*N), replace=FALSE, prob=Fs( seq(0,max(s),length=N) ) )
		X <- X[choise,]
	} else choise <- 0
	
	#print intrinsic dimension
	if( id ) cat("intrinsic dim: mean=",mean(intr_dim),", mode=",names(sort(-table(intr_dim)))[1],"\n",sep="")
	#print(intr_dim)
	
	return(list(X=X, wgts=wgts, choise=choise,id=intr_dim))
}

