
#################################################################
# smirt model
# noncompensatory and compensatory IRT models
smirt <- function( dat , Qmatrix , irtmodel="noncomp" ,
	est.b=NULL , est.a= NULL , 
	est.c=NULL , est.d=NULL , est.mu.i =NULL , b.init=NULL , a.init=NULL ,
	c.init=NULL , d.init=NULL , mu.i.init=NULL , Sigma.init =NULL , 
	theta.k=seq(-6,6,len=20) , 	theta.kDES = NULL , qmcnodes= 0 ,  
	mu.fixed=NULL , variance.fixed=NULL , 
	est.corr=FALSE , max.increment=1 , increment.factor=1 , 
	numdiff.parm=.0001 , maxdevchange=.10 ,
	globconv=.001 , maxiter=1000 , msteps=4 , mstepconv=.001){
	#..........................................................
	s1 <- Sys.time()
	dat0 <- dat <- as.matrix(dat)
	N <- nrow(dat)
	# check inits for item parameters
	res1 <- .smirt.check.inits( a.init , b.init , irtmodel , Qmatrix )	
	a.init <- res1$a.init
	b.init <- res1$b.init
	
	#***
	# process data
	# Q matrix and theta design matrix
	Q <- Qmatrix
	if ( is.null( theta.kDES) ){
		D <- ncol(Qmatrix) 
		TD <- D
			} else { 
			D <- ncol(theta.k) 
			TD <- ncol(theta.kDES )
				}
	if ( is.null(colnames(Qmatrix) ) ){
		colnames(Qmatrix) <- paste0( "Dim" , 1:D )
			}	
    if ( is.null(theta.kDES) ){			
		nodes <- theta.k
		theta.k <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, D) , 
				ncol = D ) ) ) )	
				}
	theta.k <- as.matrix( theta.k )
    Q <- as.matrix(Q)
    Qmatrix <- as.matrix(Qmatrix)	
	QMC <- FALSE
	if (qmcnodes>0){ 
		QMC <- TRUE
		theta.k <- qmc.nodes( snodes=qmcnodes , ndim=D )
#		if (! is.null(Sigma.init) ){
#			theta.k <- theta.k %*% chol(Sigma.init)
#						}
		colnames(theta.k) <- paste0("Dim",1:D)
		theta.k0 <- theta.k 
			}		
	TP <- nrow(theta.k)
	# theta design matrix
	if ( is.null(theta.kDES) ){ theta.kDES <- theta.k }
	theta.kDES <- as.matrix( theta.kDES )	
	# K <- max( dat , na.rm=T )
	K <- 1
	dat2 <- dat
	I <- ncol(dat)
	dat2.resp <- 1 - is.na(dat)
	dat2[ dat2.resp==0 ] <- 0
	# indicator matrix
	dat2.ind0 <- dat2.resp * 1*(dat2==0)
	dat2.ind1 <- dat2.resp * 1*(dat2==1)
	dat2.ind <- as.matrix( cbind( dat2.ind0 , dat2.ind1 ) )
	#*** indicators for calculating probabilities
#	index_ind <- matrix( 0 , N , I )
#	index_indL <- rep(0,I)
#	for (ii in 1:I){
		#ii <- 1
#		vii <- which( dat2.resp[,ii] == 1 )
#		index_indL[ii] <- length(vii)
#		index_ind[ 1:(index_indL[ii]) , ii] <- vii - 1
#				}
	#****
	# init b parameters
	if ( is.null(b.init)){ 
		if (irtmodel!="comp"){
			b <- Qmatrix * matrix( - stats::qlogis( colMeans(dat , na.rm =TRUE ) ) ,I , D )
				} else {	# irtmodel="comp"s
			b <-  - stats::qlogis( colMeans(dat , na.rm =TRUE ) )
						}
					} else { 
				b <- b.init 
					}
	if ( ! is.null(est.b) & is.null(b.init) ){ b <- 0*b } 
	e1 <- matrix( 1:I , I , ncol(Qmatrix) )
	if ( is.null(est.b)){ 	
		if (irtmodel!="comp"){ est.b <- Qmatrix*e1  }
		if (irtmodel=="comp"){ est.b <- 1:I  }		
				}			
	#****
	# init a parameters
	if ( is.null(a.init) ){ 
		a <- matrix( 1 , I , ncol(Qmatrix) )
		a[ Qmatrix == 0 ] <- 0
			} else { a <- a.init }
#	max.a.increment <- a*max.increment
	max.a.increment <- matrix( .4 , nrow=1 , ncol=ncol(a) )
	if ( ! is.null( est.a)){	
		if ( est.a=="2PL"){ 
			est.a <- matrix(1:I , I , ncol(Qmatrix) ) * ( Qmatrix != 0 )		
			est.corr <- TRUE
				}	
			}
	if ( is.null(est.a)){ 
		est.a <- 0 * Qmatrix*e1 
				}
	if ( ! is.null(est.a) & is.null(a.init) ){ a <- 1*a } 				
	#***
	# init c and d parameters
	if ( is.null(c.init)){  c <- rep(0,I)  } else { c <- c.init }
	if ( is.null(d.init)){ 	d <- rep(1,I) } else { d <- d.init }
	# init mu parameters for partially compensatory model
	if ( is.null(mu.i.init)){ 	
			mu.i <- rep(.5,I) 
			mu.i[ ( rowSums(Qmatrix>0) == 1 ) ] <- 0
				} else {
			mu.i <- mu.i.init 
			}
	if ( irtmodel=="partcomp"){
		if ( is.null(est.mu.i) ){
			est.mu.i <- ( 1:I ) * ( rowSums(Qmatrix>0) > 1 )
								}
							}
	
	# calculate pi.k
	mu <- rep(0,D)
	if ( is.null(mu.fixed) ){ 	mu.fixed <- cbind( 1:D , 0 ) }
	if ( is.null(Sigma.init) ){
			Sigma <- .4+.6*diag(1,D) } else { Sigma <- Sigma.init }
	if ( ! is.null(variance.fixed ) ){
			Sigma[ variance.fixed[,1:2,drop=FALSE] ] <- variance.fixed[,3]
			Sigma[ variance.fixed[,c(2,1),drop=FALSE] ] <- variance.fixed[,3]		
							}	
						
	pi.k <- mvtnorm::dmvnorm( theta.k , mean=mu , sigma=Sigma )
	pi.k <- pi.k / sum( pi.k)
	# init counts
	n.ik <- array( 0 , dim=c(TP,I , K+1) )
	probs <- array(0 , dim=c(I,K+1,TP))	
	se.mu.i <- se.d <- se.c <- se.b <- se.a <- NULL
	# inits
	iter <- 0
	dev0 <- dev <- 0
	conv <- devchange <- 1000
	disp <- "...........................................................\n"	
	#****************************************************
	# start EM algorithm
    while( ( ( maxdevchange < devchange ) | (globconv < conv) ) &
			( iter < maxiter )
						){
# zz0 <- Sys.time()						
		cat(disp)	
		cat("Iteration" , iter+1 , "   " , paste( Sys.time() ) , "\n" )	
		
		# previous values
		dev0 <- dev
		b0 <- b
		a0 <- a
		c0 <- c
		d0 <- d
		mu0 <- mu
		mu.i0 <- mu.i
		Sigma0 <- Sigma	
		# calculate probabilities
		probres <- .smirt.calcprob( a , b, Q=Qmatrix , thetak=theta.kDES , c , d ,
					mu.i , 	irtmodel )
		probs <- problong2probarray( probres , I , TP )
# cat("smirt calcprob") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1		
		# calculate likelihood and posterior
		if (QMC){ pi.k <- 1+0*pi.k }
		probsM <- matrix( aperm( probs , c(2,1,3) ) , nrow=I*2 , ncol=TP )
		res1 <- calcpost( dat2 , dat2.resp , probsM , dat2.ind , pi.k , K )
# cat("smirt calcpost") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1		
#		res2 <- calcpost2( dat2 , dat2.resp , probsM , dat2.ind , pi.k , K ,
#				index_ind , index_indL)
#cat("smirt calcpost2") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1		
		f.yi.qk <- res1$fyiqk
		f.qk.yi <- res1$f.qk.yi		
		pi.k <- res1$pi.k
		n.ik <- array( res1$n.ik , dim=dim(n.ik) )
		N.ik <- res1$N.ik
		# estimate b parameters
		res2 <- .smirt.est.b(   b , a , c , d , mu.i , Qmatrix , est.b , 
				theta.k=theta.kDES , 
				n.ik , I , K , TP , D ,  numdiff.parm=numdiff.parm , 
				max.increment=max.increment,
				msteps ,  mstepconv , irtmodel , increment.factor )		
		b <- res2$b
		se.b <- res2$se.b
#		ll <- res2$ll
# cat("smirt est.b") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1		

		# estimate a parameters
		if ( sum (est.a) > 0 ){
			res2 <- .smirt.est.a(   b , a , c , d , mu.i , Qmatrix , est.a , 
			    theta.k=theta.kDES , 
				n.ik , I , K , TP , D=TD ,  numdiff.parm=numdiff.parm , 
				max.a.increment=max.a.increment,
				msteps ,  mstepconv , irtmodel , increment.factor )		
			a <- res2$a
			se.a <- res2$se.a
							}	

						
		# estimate c parameters
		if ( sum (est.c) > 0 ){
			res2 <- .smirt.est.c(   b , a , c , d , mu.i , Qmatrix , est.c , 
				theta.k=theta.kDES , 
				n.ik , I , K , TP , D ,  numdiff.parm=numdiff.parm , 
				max.increment=max.increment,
				msteps ,  mstepconv, irtmodel , increment.factor)		
			c <- res2$c
			se.c <- res2$se.c
							}								

		# estimate d parameters
		if ( sum (est.d) > 0 ){
			res2 <- .smirt.est.d(   b , a , c , d , mu.i , Qmatrix , est.d , 
				theta.k=theta.kDES , 
				n.ik , I , K , TP , D ,  numdiff.parm=numdiff.parm , 
				max.increment=max.increment,
				msteps ,  mstepconv , irtmodel , increment.factor )		
			d <- res2$d
			se.d <- res2$se.d
							}								
							
		# estimate mu.i parameters
		if ( sum(est.mu.i) > 0 ){
			res2 <- .smirt.est.mu.i.partcomp(   b , a , c , d , mu.i , Qmatrix , est.mu.i , theta.k , 
				n.ik , I , K , TP , D , numdiff.parm , max.increment=max.increment,
				msteps ,  mstepconv  , increment.factor)
			mu.i <- res2$mu.i
			se.mu.i <- res2$se.mu.i
							}							
#		 print(a1-a0) ; a0 <- a1 ;									
		flush.console()		

		#***
		# decrease increments in every iteration
		if( increment.factor > 1){	
			max.increment <-  max.increment / increment.factor
			max.a.increment <-  max.a.increment / increment.factor
					}
		
		# Covariance estimation
		m1 <- .smirt.est.covariance( f.qk.yi , Sigma , theta.k , N ,
			mu.fixed , variance.fixed , D , est.corr , irtmodel )
		mu <- m1$mu
		Sigma <- m1$Sigma
		pi.k <- m1$pi.k
# cat("smirt est.covariance") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1	

		#**
		# calculate deviance
		pi.k <- matrix( pi.k , ncol=1 )	
#		ll <- sum( log( rowSums( f.yi.qk * outer( rep(1,nrow(f.yi.qk)) , pi.k ) ) ) )
		ll <- sum( log( rowSums( f.yi.qk * outer( rep(1,nrow(f.yi.qk)) , pi.k[,1] ) ) ) )
		if (QMC){
			ll <- sum( log( rowMeans( f.yi.qk ) ) )
			# beta[ rep(1,snodes) , ] +  t ( t(chol(variance)) %*% t(theta0.samp) )			
			theta.k <- matrix( mu , nrow=TP , ncol=D , byrow=T ) + t ( t(chol(Sigma)) %*% t(theta.k0) )
			theta.kDES <- theta.k
				} 				
		dev <- -2*ll
		# convergence criteria
		conv <- max( abs(a-a0) , abs(b-b0)  , abs(c-c0) , abs(d-d0) , abs(mu.i - mu.i0) )
		iter <- iter+1
		devchange <- abs( ( dev - dev0 ) / dev0  )
		#****
		# print progress			
		cat( paste( "   Deviance = "  , round( dev , 4 ) , 
			if (iter > 1 ){ " | Deviance change = " } else {""} ,
			if( iter>1){round( - dev + dev0 , 6 )} else { ""}	,"\n",sep="") )
		cat( paste( "    Maximum b parameter change = " , 
				paste( round(max(abs(b-b0)) ,6) , collapse=" " ) , "\n" , sep=""))
		cat( paste( "    Maximum a parameter change = " , 
				paste( round(max(abs(a-a0)) ,6) , collapse=" " ) , "\n" , sep=""))
		cat( paste( "    Maximum c parameter change = " , 
				paste( round(max(abs(c-c0)) ,6) , collapse=" " ) , "\n" , sep=""))
		cat( paste( "    Maximum d parameter change = " , 
				paste( round(max(abs(d-d0)) ,6) , collapse=" " ) , "\n" , sep=""))				
		if (irtmodel=="partcomp"){
			cat( paste( "    Maximum mu.i parameter change = " , 
				paste( round(max(abs(mu.i-mu.i0)) ,6) , collapse=" " ) , "\n" , sep=""))				
					}
				
		#--- distribution parameters
		cat( " Means: " , round( mu , 3 ) , "\n")				
		cat( " Standard deviations: " , round( sqrt(diag(Sigma)) , 3 ) , "\n")
		c1 <- stats::cov2cor(Sigma)
		cat( " Correlations " , round( c1[lower.tri(c1)] , 3 ) , "\n")		
		# flush.console()			
				}

			
	# *********
	# arrange OUTPUT
	#---
	# Information criteria
	ic <- list( "deviance" = dev , "n" = nrow(dat2) )
	ic$np.item <- 0
	ic$np.item.b <- ic$np.item.a <- ic$np.item.c <- ic$np.item.d <- 0
	if (irtmodel=="noncomp"){
	for (dd in 1:D){ 
		ic$np.item.b <- ic$np.item.b + length( setdiff( unique( est.b[,dd] ) , 0 )  )
				}	} else {
		ic$np.item.b <- ic$np.item.b + length( setdiff( unique( est.b ) , 0 )  )
						}
	for (dd in 1:TD){ 				
		ic$np.item.a <- ic$np.item.a + length( setdiff( unique( est.a[,dd] ) , 0 ) )		
			}
	ic$np.item.c <- length( setdiff( unique( est.c ) , 0 ) )			
	ic$np.item.d <- length( setdiff( unique( est.d ) , 0 ) )	
	ic$np.item.mu.i <- length( setdiff( unique( est.mu.i ) , 0 ) )	
	# sum over all item parameters
	ic$np.item <- ic$np.item.b + ic$np.item.a + ic$np.item.c + ic$np.item.d + ic$np.item.mu.i
	# covariance matrix
	ic$np.cov <- 0
	
	ic$np.cov.covM <- D*(D+1)/2	
	if ( est.corr ){ ic$np.cov.covM <- ic$np.cov.covM - D }
	if ( ! is.null(variance.fixed) ){ 
			ic$np.cov.covM <- ic$np.cov.covM - nrow(variance.fixed) 
									}
	ic$np.cov.mu <- D
	if ( ! is.null(mu.fixed) ){ ic$np.cov.mu <- ic$np.cov.mu - nrow(mu.fixed) }
	ic$np.cov <- ic$np.cov.covM + ic$np.cov.mu
	ic$np <- ic$np.item + ic$np.cov
    # AIC
    ic$AIC <- dev + 2*ic$np
    # BIC
    ic$BIC <- dev + ( log(ic$n) )*ic$np
    # CAIC (conistent AIC)
    ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
	# corrected AIC
    ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )		
	#---
	# person parameters
	pers <- .smirt.person.parameters( data=dat2 , D=D , theta.k=theta.k ,
		p.xi.aj=f.yi.qk , p.aj.xi=f.qk.yi , weights=rep(1,N) )	
	person <- pers$person
	EAP.rel <- pers$EAP.rel 
	#---
	# item
	item <- data.frame("item" = colnames(dat))
	item$N <- colSums( dat2.resp )
	item$p <- colSums( dat2 , na.rm=TRUE) / item$N
	if ( irtmodel != "comp"){
		for (dd in 1:( ncol(Qmatrix)) ){
			b[ Qmatrix[,dd] == 0 , dd ] <- NA
			se.b[ est.b[,dd] == 0 , dd ] <- NA		
			item[ , paste0("b.",colnames(Qmatrix)[dd]) ] <- b[,dd]
						}
				} else {	# end irtmodel == "noncomp"
		item$b <- b
		se.b[ est.b == 0 ] <- NA
				} # end irtmodel != "noncomp"
	for (dd in 1:( ncol(Qmatrix)) ){
		a[ Qmatrix[,dd] == 0 , dd ] <- NA
		if ( sum( est.a) != 0 ){
			se.a[ est.a[,dd] == 0 , dd ] <- NA		
						}
		item[ , paste0("a.",colnames(Qmatrix)[dd]) ] <- a[,dd]
					}
	item$c <- c
	item$d <- d
	if ( irtmodel == "partcomp" ){
	    mu.i[ rowSums( Qmatrix > 0 ) == 1 ] <- NA
		item$mu.i <- mu.i
				}
	
	
	obji <- item
	for (vv in seq(2,ncol(obji) )){
		obji[,vv] <- round( obji[,vv],3 ) }
	cat("*********************************\n")
	cat("Item Parameters\n")
    print( obji )		
	
	cat("*********************************\n")
	cat("EAP Reliability = " , round(EAP.rel,3) , "\n")		
	
	s2 <- Sys.time()
	
    res <- list("deviance" = dev , "ic"=ic , "item"=item , 
		"person" = person , "EAP.rel"=EAP.rel , 
		"mean.trait"=mu , "sd.trait" = sqrt( diag(Sigma ) ) ,
		"Sigma"=Sigma , "cor.trait"= cov2cor(Sigma ) , 
		"b"=b , "se.b" = se.b , "a"=a , "se.a"=se.a ,
		"c"=c , "se.c" = se.c , "d"=d , "se.d"=se.d ,	
		"mu.i"=mu.i , 	"se.mu.i" = se.mu.i , 
		"f.yi.qk"=f.yi.qk , "f.qk.yi"=f.qk.yi , "probs"=probs ,
		"n.ik"=n.ik ,  "iter"=iter , "dat2"=dat2 , "dat2.resp"=dat2.resp , 
		"dat" = dat0 , 
		"I"=I , "D"=D , "K"=K ,  "G"= 1 , "theta.k" = theta.k , "pi.k"=pi.k , 
		"irtmodel"=irtmodel , "Qmatrix" = Qmatrix , "s1"=s1 , "s2"=s2 
				)
	class(res) <- "smirt"
	return(res)
		}
############################################################