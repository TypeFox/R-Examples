


######################################################################
# calculation of pairwise marginal likelihood
.ll.rasch.pml3.probit.est.b <- function( b , a , sigma , Q ,eps.corr, itempairs  , IP , eps=10^(-14) ,
		h , desb00 , desb01 , desb10  , desb11 , b.items ){
	cor.Sigma <- NULL
	# unidimensional case
	t1 <- a^2*sigma^2	
	# define different b parameters
	xib0 <-   - b  / sqrt( 1 + t1 )	
	xib1 <-   - (b+h)  / sqrt( 1 + t1 )	
	xib2 <-   - (b-h)  / sqrt( 1 + t1 )	
# a00 <- Sys.time()	
	# 00
	res00 <- .pml3.est.b.aux( xib1=xib0 , xib2=xib0 , itempairs , a , sigma , 
				eps.corr , cor.Sigma , eps )	
	ll00 <- res00$ll

	# 10 / 01		
	# 10: b + h ; b
	# 01: b ; b + h
	# 11: b+h ; b+h
	ll10 <- .pml3.est.b.aux( xib1=xib1 , xib2=xib0 , itempairs , a , sigma , 
				eps.corr , cor.Sigma , eps )$ll	
	ll01 <- .pml3.est.b.aux( xib1=xib0 , xib2=xib1 , itempairs , a , sigma , 
				eps.corr , cor.Sigma , eps )$ll	
	ll11 <- .pml3.est.b.aux( xib1=xib1 , xib2=xib1 , itempairs , a , sigma , 
				eps.corr , cor.Sigma , eps )$ll				
	# 20 / 02
	# 20: b-h;b
	# 02: b;b-h
	# 22: b-h;b-h
	ll20 <- .pml3.est.b.aux( xib1=xib2 , xib2=xib0 , itempairs , a , sigma , 
				eps.corr , cor.Sigma , eps )$ll	
	ll02 <- .pml3.est.b.aux( xib1=xib0 , xib2=xib2 , itempairs , a , sigma , 
				eps.corr , cor.Sigma , eps )$ll	
	ll22 <- .pml3.est.b.aux( xib1=xib2 , xib2=xib2 , itempairs , a , sigma , 
				eps.corr , cor.Sigma , eps )$ll								
# cat("calc all ll") ; a11 <- Sys.time(); print(a11-a00) ; a00 <- a11		
# 	ll.temp <- as.matrix( cbind( ll00 , ll10 , ll01 , ll20 , ll02 ) )
	ll0 <- t(ll00) %*% desb00
	ll1 <- t(ll01) %*% desb01 +  t(ll10) %*% desb10 + t(ll11) %*% desb11
	ll2 <- t(ll02) %*% desb01 +  t(ll20) %*% desb10 + t(ll22) %*% desb11
	# calculate increment			
	incr <- nr.numdiff( ll0=as.vector(ll0) , ll1=as.vector(ll1) , 
			ll2=as.vector(ll2) , h=h , eps = 10^(-10) )	
# cat("calc incr") ; a11 <- Sys.time(); print(a11-a00) ; a00 <- a11	
	change <- rep( 0 , length(b) )
	change <- incr[ b.items ]		
	change[ is.na(change ) ] <- 0
	change <- ifelse( abs(change) > 1 , 1*sign(change) , change )
	b <- b + change
	ll <- sum( ll00 )
# cat("add and mult") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		        						
    res <-list( "ll" = ll , "itempairs" = res00$itempairs , "b" = b, "sigma" = sigma # ,
#                "ind1" = ind1 , "ind2" = ind2 
				)
    return(res)
        }

######################################################################


#*********	
.pml3.est.b.aux <- function( xib1 , xib2 , itempairs , a , sigma , 
			eps.corr , cor.Sigma , eps ){
# a0 <- Sys.time()			
    xi1 <- xib1[ itempairs[,"item1"] ] 
    xi2 <- xib2[ itempairs[,"item2"] ] 
	a1 <- a[ itempairs[,"item1"] ] 	
	a2 <- a[ itempairs[,"item2"] ] 
	t1 <- a1*a2*sigma^2	
    cor.Sigma <- ( t1 + eps.corr ) / ( 1 + t1 )
	cor.Sigma[ cor.Sigma > 1 ] <- .99
# cat("cor sigma") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		        	
	pxi1 <- stats::pnorm( xib1 )
	pxi2 <- stats::pnorm( xib2 )	
# cat("pnorm") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		        		
    itempairs$p1.item1 <- pxi1[ itempairs$item1 ]
    itempairs$p1.item2 <- pxi2[ itempairs$item2 ]
	itempairs$p11 <- pbivnorm2( x = xi1 , y = xi2 , rho = cor.Sigma )
# cat("pbivnorm") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		        			
    itempairs$p10 <- itempairs$p1.item1 - itempairs$p11
    itempairs$p01 <- itempairs$p1.item2 - itempairs$p11
    itempairs$p00 <- 1 - itempairs$p11 - itempairs$p01 - itempairs$p10
    ind1 <- which( colnames(itempairs) %in% c( "p11" , "p10" , "p01" , "p00" ) )
# cat("ip2") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		        				
#    itempairs[ itempairs[ , ind1 ] < eps , ind1 ] <- eps
	for (ii in ind1){
		itempairs[ , ii ] <- ifelse( itempairs[,ii] < eps , eps , itempairs[,ii] )
				}
# cat("ifelse") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		        								
    ind2 <- which( colnames(itempairs) %in% c( "f11" , "f10" , "f01" , "f00" ) )
#    ll <-  sum( log( itempairs[,ind1] ) * itempairs[,ind2]  )
    ll <-  rowSums( log( itempairs[,ind1] ) * itempairs[,ind2]  )
# cat("ll") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		
	res <- list("ll" = ll , "itempairs" = itempairs )
# cat("rest") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		        									
# stop("here")
	return(res)
		}
#********************************






######################################################################
# calculation of pairwise marginal likelihood
.ll.rasch.pml3.probit.est.a <- function( b , a , sigma , Q ,eps.corr, itempairs  , IP , eps=10^(-14) ,
		h , desa00 , desa01 , desa10  , desa11 , a.items ){
	cor.Sigma <- NULL
    a0 <- a
	a1 <- a+h
	a2 <- a-h
		
	# 00
	ll00 <- .pml3.est.a.aux( b , itempairs , a01=a0 , a02=a0 , sigma , 
			eps.corr , cor.Sigma , eps )								
	# 10 / 01 / 11			
	ll10 <- .pml3.est.a.aux( b , itempairs , a01=a1 , a02=a0 , sigma , 
			eps.corr , cor.Sigma , eps )	
	ll01 <- .pml3.est.a.aux( b , itempairs , a01=a0 , a02=a1 , sigma , 
			eps.corr , cor.Sigma , eps )		
	ll11 <- .pml3.est.a.aux( b , itempairs , a01=a1 , a02=a1 , sigma , 
			eps.corr , cor.Sigma , eps )		
	# 20 / 02 / 22
	ll20 <- .pml3.est.a.aux( b , itempairs , a01=a2 , a02=a0 , sigma , 
			eps.corr , cor.Sigma , eps )	
	ll02 <- .pml3.est.a.aux( b , itempairs , a01=a0 , a02=a2 , sigma , 
			eps.corr , cor.Sigma , eps )		
	ll22 <- .pml3.est.a.aux( b , itempairs , a01=a2 , a02=a2 , sigma , 
			eps.corr , cor.Sigma , eps )		
	ll0 <- t(ll00) %*% desa00
	ll1 <- t(ll01) %*% desa01 +  t(ll10) %*% desa10 + t(ll11) %*% desa11
	ll2 <- t(ll02) %*% desa01 +  t(ll20) %*% desa10 + t(ll22) %*% desa11
	# calculate increment			
	incr <- nr.numdiff( ll0=as.vector(ll0) , ll1=as.vector(ll1) , ll2=as.vector(ll2) , h=h , eps = 10^(-10) )	
	change <- rep( 0 , length(a) )
	change <- incr[ a.items ]		
	change[ is.na(change ) ] <- 0
	change <- ifelse( abs(change) > .3 , .3*sign(change) , change )	
	a <- a + change
	ll <- sum( ll00 )
    res <-list( "ll" = ll , "itempairs" = itempairs , "a" = a, "sigma" = sigma # ,
#                "ind1" = ind1 , "ind2" = ind2 
				)
    return(res)
        }

######################################################################


#*********	
.pml3.est.a.aux <- function( b , itempairs , a01 , a02 , sigma , 
			eps.corr , cor.Sigma , eps ){
	
	# unidimensional case
	t11 <- a01^2*sigma^2	
	t12 <- a02^2*sigma^2	
	
	# define different b parameters
	xib1 <-  - b  / sqrt( 1 + t11 )				
	xib2 <-  - b  / sqrt( 1 + t12 )				
	
    xi1 <- xib1[ itempairs[,"item1"] ] 
    xi2 <- xib2[ itempairs[,"item2"] ] 

	a1 <- a01[ itempairs[,"item1"] ] 	
	a2 <- a02[ itempairs[,"item2"] ] 
	
	t1 <- a1*a2*sigma^2	
    cor.Sigma <- ( t1 + eps.corr ) / ( 1 + t1 )
	cor.Sigma[ cor.Sigma > 1 ] <- .99
    
	pxi1 <- stats::pnorm( xib1 )
	pxi2 <- stats::pnorm( xib2 )	
	
    itempairs$p1.item1 <- pxi1[ itempairs$item1 ]
    itempairs$p1.item2 <- pxi2[ itempairs$item2 ]
	itempairs$p11 <- pbivnorm2( x = xi1 , y = xi2 , rho = cor.Sigma )
    itempairs$p10 <- itempairs$p1.item1 - itempairs$p11
    itempairs$p01 <- itempairs$p1.item2 - itempairs$p11
    itempairs$p00 <- 1 - itempairs$p11 - itempairs$p01 - itempairs$p10
    ind1 <- which( colnames(itempairs) %in% c( "p11" , "p10" , "p01" , "p00" ) )
    itempairs[ itempairs[ , ind1 ] < eps , ind1 ] <- eps
    ind2 <- which( colnames(itempairs) %in% c( "f11" , "f10" , "f01" , "f00" ) )
#    ll <-  sum( log( itempairs[,ind1] ) * itempairs[,ind2]  )
    ll <-  rowSums( log( itempairs[,ind1] ) * itempairs[,ind2]  )
	return(ll)
		}
#********************************







######################################################################
# calculation of pairwise marginal likelihood
.ll.rasch.pml3.probit.est.corr <- function( b , a , sigma , Q ,eps.corr, itempairs  , 
		IP , eps=10^(-14) , h , deseps00 , eps.items ){
	cor.Sigma <- NULL		
	# 00
	ll00 <- .pml3.est.eps.aux( b , itempairs , a , sigma , 
			eps.corr , cor.Sigma , eps )								
	# 11
	ll11 <- .pml3.est.eps.aux( b , itempairs , a , sigma , 
			eps.corr+h , cor.Sigma , eps )								
	# 20 / 02 / 22
	ll22 <- .pml3.est.eps.aux( b , itempairs , a , sigma , 
			eps.corr-h , cor.Sigma , eps )								
	ll0 <- t(ll00) %*% deseps00
	ll1 <- t(ll11) %*% deseps00
	ll2 <- t(ll22) %*% deseps00
	# calculate increment			
	incr <- nr.numdiff( ll0=as.vector(ll0) , ll1=as.vector(ll1) , ll2=as.vector(ll2) , h=h , eps = 10^(-10) )	
	change <- rep( 0 , IP )
	change <- incr[ eps.items ]		
	change[ is.na(change ) ] <- 0
	change <- ifelse( abs(change) > .3 , .3*sign(change) , change )	
	eps.corr <- eps.corr + change
	ll <- sum( ll00 )
    res <-list( "ll" = ll , "itempairs" = itempairs , "eps.corr" = eps.corr, "sigma" = sigma # ,
#                "ind1" = ind1 , "ind2" = ind2 
				)
    return(res)
        }

######################################################################


#*********	
.pml3.est.eps.aux <- function( b , itempairs , a , sigma , 
			eps.corr , cor.Sigma , eps ){
	
	# unidimensional case
	t11 <- a^2*sigma^2	
	t12 <- a^2*sigma^2	
	
	# define different b parameters
	xib1 <-  - b  / sqrt( 1 + t11 )				
	xib2 <-  - b  / sqrt( 1 + t12 )				
	
    xi1 <- xib1[ itempairs[,"item1"] ] 
    xi2 <- xib2[ itempairs[,"item2"] ] 

	a1 <- a[ itempairs[,"item1"] ] 	
	a2 <- a[ itempairs[,"item2"] ] 
	
	t1 <- a1*a2*sigma^2	
    cor.Sigma <- ( t1 + eps.corr ) / ( 1 + t1 )
	cor.Sigma[ cor.Sigma > 1 ] <- .99
    
	pxi1 <- stats::pnorm( xib1 )
	pxi2 <- stats::pnorm( xib2 )	
	
    itempairs$p1.item1 <- pxi1[ itempairs$item1 ]
    itempairs$p1.item2 <- pxi2[ itempairs$item2 ]
	itempairs$p11 <- pbivnorm2( x = xi1 , y = xi2 , rho = cor.Sigma )
    itempairs$p10 <- itempairs$p1.item1 - itempairs$p11
    itempairs$p01 <- itempairs$p1.item2 - itempairs$p11
    itempairs$p00 <- 1 - itempairs$p11 - itempairs$p01 - itempairs$p10
    ind1 <- which( colnames(itempairs) %in% c( "p11" , "p10" , "p01" , "p00" ) )
    itempairs[ itempairs[ , ind1 ] < eps , ind1 ] <- eps
    ind2 <- which( colnames(itempairs) %in% c( "f11" , "f10" , "f01" , "f00" ) )
#    ll <-  sum( log( itempairs[,ind1] ) * itempairs[,ind2]  )
    ll <-  rowSums( log( itempairs[,ind1] ) * itempairs[,ind2]  )
	return(ll)
		}
#********************************