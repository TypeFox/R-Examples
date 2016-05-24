


####################################################
# probabilistic Guttman analysis
prob.guttman <- function( dat , pid = NULL , guess.equal=FALSE ,
		slip.equal=FALSE , itemlevel = NULL , conv1 = .001 , 
		glob.conv = .001 , mmliter = 500 ){
	# data prep
	dat <- as.matrix(dat)
	I <- ncol(dat)
	dat2.resp <- 1 - is.na( dat )
	dat2 <- dat
	dat2[ is.na(dat2) ] <- 0
	if ( is.null( pid) ){ pid <- seq(1,nrow(dat)) }
	# p values
	p <- colMeans( dat , na.rm=TRUE )
	# items and p value levels
	itemtable <- data.frame( "index" = seq(1,ncol(dat)) , 
					"item" = colnames(dat) , "p" = p )
	if ( is.null( itemlevel )){ 
		itemtable <- itemtable[ order( itemtable$p , decreasing=TRUE ) , ]
		itemtable$level <- match( itemtable$p , unique(itemtable$p ) )
			} else { itemtable$level <- itemlevel }
	itemtable <- itemtable[ order( itemtable$index ) , ]
	theta.k <- c( 0 , sort(unique(itemtable$level ) ))
	# design matrices for items
	itemdes <- matrix(0 , nrow= length(theta.k) , ncol=ncol(dat))
	colnames(itemdes) <- colnames(dat)
	rownames(itemdes) <- paste( "theta_level_" , theta.k , sep="")
	items <- itemtable$item
	for (ii in items){ 
		itemdes[,ii] <- 1 * ( theta.k >= itemtable[ itemtable$item == ii , "level" ] )
					}
	# init slipping and guessing parameters
	guess <- slip <- rep(.2,I)

	# prior distribution pi(k)
	LK <- length(theta.k)
	pi.k <- rep( 1 /LK , LK )
	iter <- 0
	dev0 <- dev.change <- 1000
	par.change <- 1000
    cat("---------------------------------------------------------------------------------------------------------- \n")
    cat("Probabilistic Guttman Model \n")
    cat("---------------------------------------------------------------------------------------------------------- \n")
	while ( ( dev.change > glob.conv | par.change > conv1 ) & ( iter < mmliter ) ){
		#*****
		# E step Guttman scaling
		res <- .e.step.guttman( dat2 , dat2.resp , theta.k , pi.k , I , guess , slip , itemdes )
		
		# M step Guttman scaling
		pi.k <- ( res$n.k / sum(res$n.k ) )[,1]
		
		guess.old <- guess
		slip.old <- slip
		# update guessing parameters
		
		c1 <- colSums( t( (res$r.jk)[,,1]  ) * ( 1 - itemdes ) )
		c2 <- colSums( t( (res$n.jk)[,,1]  ) * ( 1 - itemdes ) )
		if ( ! guess.equal ){	guess <- c1/c2 } else {
			guess <-  rep( sum(c1) / sum(c2) , I )
						}
		# update slipping parameters
		c1 <- colSums( t( (res$r.jk)[,,1]  ) * (  itemdes ) )
		c2 <- colSums( t( (res$n.jk)[,,1]  ) * ( itemdes ) )
#		slip <- 1- c1/c2
		if ( ! slip.equal ){	slip <- 1 - c1/c2 } else {
			slip <-  rep( 1 - sum(c1) / sum(c2) , I )
						}
		par.change <- max( abs( c( guess - guess.old , slip - slip.old) ) )
		dev <- -2*res$ll
		dev.change <- abs( dev - dev0 )
		cat( paste( "Iteration " , iter + 1 ,  " |  Deviance = "  , round( dev , 4 ) , 
									if (iter > 0 ){ " | Deviance change = " } else {""} ,
									if( iter>0){round( - dev + dev0 , 6 )} else { ""}  ," |",sep=""))                       
		cat( paste( " Maximum parameter change = " , 
									round( max(abs( par.change )) , 6 ) ,  " \n"   )  )  
		utils::flush.console()
		dev0 <- dev                                    
		iter <- iter + 1
			}

	# MAP
	map <- theta.k[ whichrowMaxs( res$f.qk.yi )$arg ]
	# MLE
	mle <- theta.k[ whichrowMaxs( res$f.yi.qk )$arg ]

	# person data frame
	person <- data.frame( "pid" = pid ,  
				"score" = rowSums( dat , na.rm=T) , "max" = rowSums(dat2.resp) , 
				"MLE" = mle , "MAP" = map )
	# Information criteria
	if (guess.equal){ np <- 1 } else {
		np <- length(guess) }
	if (slip.equal){ np <- np + 1 } else { np <- np + length(slip) }
	np <- np + length(theta.k) - 1
	ic <- data.frame( "np" = np ,"n" = nrow(dat),
		"deviance" = dev )
    # AIC
    ic$AIC <- dev + 2*ic$np
    # BIC
    ic$BIC <- dev + ( log(ic$n) )*ic$np
    # CAIC (conistent AIC)
    ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
	# corrected AIC
    ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )		
	# item data frame
	item <- itemtable          
	item$guess <- guess
	item$slip <- slip            
	item$discrim <- 1-slip-guess
	#**********
	# trait distribution
	trait <- data.frame( "level" = theta.k , "prob" = pi.k )
	trait$freq.MLE <- sapply( theta.k , FUN = function(tt){ mean( person$MLE == tt ) } )
	trait$freq.MAP <- sapply( theta.k , FUN = function(tt){ mean( person$MAP == tt ) } )
	rownames(trait) <- rownames(itemdes)
	# collect results
	res <- list( "person" = person , "item" =item , "theta.k" = theta.k ,
		    "trait" = trait , "pi.k"= pi.k , "f.yi.qk"=res$f.yi.qk ,
			"f.qk.yi" = res$f.qk.yi , "probs" = res$probs , 
			"ic" = ic , "G"= 1 , "deviance" = dev , "iter" = iter ,
			"itemdesign" =itemdes )
	class(res) <- "prob.guttman"
	return(res)
		}
###########################################################################




#*******************************************************
# Summary for prob.guttman object                         *
summary.prob.guttman <- function( object , ... ){
    # object      ... object from rasch.mml                #
	
    cat("---------------------------------------------------------------------------------------------------------- \n")
		d1 <- utils::packageDescription("sirt")
		cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
#		cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
#		cat("Computation time:" , print(object$s2 - object$s1), "\n\n")
    cat("Probabilistic Guttman Model \n")
    cat("---------------------------------------------------------------------------------------------------------- \n")
	cat( "Number of iterations =" , object$iter , "\n" )
    cat( "Deviance = " , round( object$deviance , 2 ) , "\n" )
    cat( "Number of persons = " , object$ic$n , "\n" )    
    cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
    cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , 
				round( object$ic$AIC - object$ic$deviance ,2 ) , 
			"   | AIC = -2*LL + 2*p  \n" )    
    cat( "AICc = " , round( object$ic$AICc , 2 ) ," | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )
		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
    cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , 
			"   | BIC = -2*LL + log(n)*p  \n" )  
    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   

	#cat( "\n")
    cat("---------------------------------------------------------------------------------------------------------- \n")
    cat("Trait Distribution \n")
	obji <- object$trait
	roundvars <- c("prob" , "freq.MLE" , "freq.MAP" )
    for (vv in roundvars ){ 
		obji[,vv] <- round( obji[,vv] , 3 ) 
				}
#	rownames(obji) <- NULL
    print( obji )                

		
	cat( "\n")
    cat("---------------------------------------------------------------------------------------------------------- \n")
    cat("Item Parameter \n")
	obji <- object$item
	roundvars <- c("p","guess","slip","discrim" )
    for (vv in roundvars ){ 
		obji[,vv] <- round( obji[,vv] , 3 ) 
				}
	rownames(obji) <- NULL
    print( obji )                
                }
#*******************************************************





#*************************************************************************************
.e.step.guttman <- function( dat2 , dat2.resp , theta.k , pi.k , I , guess , slip , itemdes ){
    #...................................                    
    # probabilities of correct item at theta_k
    guessM <- matrix( guess , byrow=T , nrow=length(theta.k) , ncol=I )
    slipM <- matrix( slip , byrow=T , nrow=length(theta.k) , ncol=I )   
    pjk <- itemdes * ( 1 - slipM ) +  (  1 - itemdes ) * guessM   
	TP <- dim(pjk)[1]
#    pjk <- t(pjk)
        #***
        pjkt <- t(pjk)
		pjkL <- array( NA , dim=c( I , 2 , TP  ) )
		pjkL[,1,] <- 1 - pjkt
		pjkL[,2,] <- pjkt	
		probsM <- matrix( aperm( pjkL , c(2,1,3) ) , nrow=I*2 , ncol=TP )
		f.yi.qk <- mml_calc_like( dat2=dat2 , dat2resp = dat2.resp , 
					probs = probsM )$fyiqk
	#******
    f.qk.yi <- 0 * f.yi.qk
    pi.k <- matrix( pi.k , ncol=1 )
    f.qk.yi <- f.yi.qk * outer( rep( 1 , nrow(dat2) ) , pi.k[,1] )        
    f.qk.yi <- f.qk.yi / rowSums( f.qk.yi )     

    # expected counts at theta.k
    G <- 1
    n.k <- matrix( 0 , nrow=length(theta.k) , ncol=G )
    r.jk <- n.jk <- array( 0 , dim=c( ncol(dat2) , length(theta.k) , G) )
    ll <- rep(0,G)
    group <- rep(1,nrow(dat2))
    dat1 <- matrix( 1 , nrow=nrow(dat2) , ncol=2 )
    for (gg in 1:G){ 
		ind.gg <- which( group == gg )
	    res <- mml_raschtype_counts( dat2=dat2[ind.gg,] , dat2resp=dat2.resp[ind.gg,] , 
					dat1=dat1[ind.gg,2] , fqkyi=f.qk.yi[ind.gg,] ,
					pik=pi.k[,gg] , fyiqk=f.yi.qk[ind.gg,]  )
		n.k[,gg] <- res$nk
		n.jk[,,gg] <- res$njk
        r.jk[,,gg] <- res$rjk
		ll[gg] <- res$ll
        }       
    res <- list( "n.k" = n.k , "n.jk" = n.jk , "r.jk" = r.jk , 
            "f.qk.yi" = f.qk.yi , "pjk" = pjk  ,
            "f.yi.qk" = f.yi.qk , "ll" = sum(ll) , "probs" = pjkL )
    return(res)
    }
#*************************************************************************************

