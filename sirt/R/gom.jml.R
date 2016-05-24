
########################################
# GOM JML
gom.jml <- function( dat , K=2 , seed=NULL , globconv = .001 , 
        maxdevchange = .001 , maxiter = 600 , min.lambda = .001 ,
		min.g =.001 ){ 
	s1 <- Sys.time()		
    dev <- 0
    e1 <- environment()
	N <- dat.resp <- p.item <- score <- g.mean <- g1 <- gcut.distr <- NULL	
    # data processing
    datproc <- res <- .gom.proc(dat, seed=seed,K)
    .sirt.attach.environment( res , envir=e1 )
    iter <- 0 ; conv <- 10
	devchange <- -100
        disp <- "...........................................................\n" 
   
    #****
    # BEGIN Algorithm
    while( ( ( maxdevchange < devchange ) | (globconv < conv) ) &
                ( iter < maxiter )
                            ){
        cat(disp)   
        cat("Iteration" , iter+1 , "   " , paste( Sys.time() ) , "\n" )                         
        g0 <- g
        lambda0 <- lambda
        dev0 <- dev
        # update g
        g <- .gom.jml.est.g( lambda , g , N , I , K , dat , dat.resp , min.g )
        # update lambda
        lambda <- .gom.jml.est.lambda( lambda , g , N , I , K ,  dat , dat.resp ,
			min.lambda )
        # compute deviance
        dev <- .gom.deviance( lambda, g , K , dat , dat.resp )
        conv <- max( max( abs( g0-g) ) , max( abs( lambda-lambda0) ) )
        devchange <- abs( ( dev - dev0 ) / dev0 )
        cat( paste( "   Deviance = "  , round( dev , 4 ) , 
            if (iter > 1 ){ " | Deviance change = " } else {""} ,
            if( iter>1){round( - dev + dev0 , 6 )} else { ""}   ,"\n",sep="") )
        cat( paste( "    Maximum membership parameter change = " , 
                    paste( round(max(abs(g-g0)) ,6) , collapse=" " ) , "\n" , sep=""))
        cat( paste( "    Maximum probability change = " , 
                    paste( round(max(abs(lambda-lambda0)) ,6) , collapse=" " ) , "\n" , sep=""))
        flush.console()
        iter <- iter+1
            }
    colnames(g) <- colnames(lambda) <- paste0("Class",1:K)
	#*****
	# information criteria
	ic <- list( "deviance" = dev , "n" = nrow(dat) )
	ic$K <- K
	ic$I <- I
	ic$np.item <- K*I
	ic$np.person <- (K-1)*N
	ic$np <- ic$np.item + ic$np.person
    # AIC
    ic$AICi <- dev + 2*ic$np.item
    ic$AICip <- dev + 2*ic$np
	ic$BICi <- dev + ( log(ic$n) )*ic$np.item 
    # BIC
    ic$BIC <- dev + ( log(ic$n) )*ic$np
    # CAIC (conistent AIC)
    ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
	# corrected AIC
    ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )		
	#***
	lambda <- data.frame("item"=colnames(dat) , "p" = p.item , lambda )	
	#***
	# discretizing membership functions
		
    res <- .gom.discret.membership(g,K,score)				
	.sirt.attach.environment( res , envir=e1 )
	
	s2 <- Sys.time()
    res <- list( "lambda"=lambda , "g"=g , "g.mean" = g.mean , 
			"gcut"=g1 , "gcut.distr"=gcut.distr , "K"=K ,
			"deviance"=dev, "ic"=ic , "N"=N , "score"=score, "iter"=iter , 
			"datproc"=datproc , 
			"s1"=s1 , "s2"=s2)
	class(res) <- "gom.jml"
    return(res)
    }
#####################################################
#####################################################
summary.gom.jml <- function( object , ... ){
	cat("-----------------------------------------------------------------\n")
		d1 <- packageDescription("sirt")
		cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
		cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
		cat("Computation time:" , print(object$s2 - object$s1), "\n\n")
    cat("Grade of Membership Model (Joint Maximum Likelihood Estimation) \n")
	cat("   Function 'gdm.jml'\n")
	cat("-----------------------------------------------------------------\n")
	cat( "Number of iterations =" , object$iter , "\n" )
    cat( "Deviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    cat( "Number of persons = " , object$ic$n , "\n" )    
    cat( "Number of items   = " , object$ic$I , "\n" )    		
    cat( "Number of classes = " , object$ic$K , "\n" )    	
    cat( "Number of estimated parameters  = " , object$ic$np , "\n" ) 
    cat( "    Item parameters (ni)        = " , object$ic$np.item , "\n" ) 
    cat( "    Person parameters (np)      = " , object$ic$np.person , "\n" ) 	
    cat( "AICi  = " , round( object$ic$AICi , 2 ) , " | penalty =" , round( object$ic$AICi - object$ic$deviance ,2 ) , 
			"   | AICi = -2*LL + 2*(ni)  \n" )    
    cat( "AICip = " , round( object$ic$AICip , 2 ) , " | penalty =" , round( object$ic$AICip - object$ic$deviance ,2 ) , 
			"   | AICip = -2*LL + 2*(ni+np)  \n" )    			
#    cat( "AICc = " , round( object$ic$AICc , 2 ) ," | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )
#		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
    cat( "BICi  = " , round( object$ic$BICi , 2 ) , " | penalty =" , round( object$ic$BICi - object$ic$deviance ,2 ) , 
			"   | BICi = -2*LL + log(n)*ni  \n" )  
#    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
#		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   
#	cat( "Trait Distribution\n" , 
#			  "Mean=" , 0 , " SD=" , round( object$sigma , 3) ) 
#	cat( "\nEAP Reliability = ") 
#	cat(round( object$EAP.rel,3 ) )
#	cat( "\n")
	cat("-----------------------------------------------------------------\n")
	cat("Item Parameters \n")
	obji <- object$lambda
	rownames(obji) <- NULL
	rvars <- seq( 2 , ncol(obji ) )
	for (vv in rvars ){ obji[,vv] <- round( obji[,vv] , 3 ) }
	print(obji)
	cat("-----------------------------------------------------------------\n")
	cat("Membership Scores\n\n")
    cat( "Class Proportions:\n" )
	print( round( object$g.mean,3)  )
    cat( "\nDistribution membership scores\n" )	
	obji <- object$gcut.distr	
#	obji <- object$rater
	rvars <- seq( 1 , ncol(obji ) )
	for (vv in rvars ){ obji[,vv] <- round( obji[,vv] , 2 ) }
	print(obji)	
                }
#*******************************************************
