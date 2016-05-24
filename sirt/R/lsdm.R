
#.......................................................................................................................#
# LSDM - Least Squares Distance Method 
# LSDM -- Least Squares Distance Method of Cognitive Validation                           
# Reference: Dimitrov, D. (2007) Applied Psychological Measurement, 31, 367-387.         
lsdm <- function( data , Qmatrix , theta = qnorm(seq(.0005,.9995,len=100)) , quant.list = c( .5 , .65 , .8 ) ,
                    b = NULL , a = rep( 1 , nrow(Qmatrix) ) ,
                    c = rep(0, nrow(Qmatrix) )  ){
        #...................................................#
        # Input:                                            #
        # data ... ( I x L ) matrix of ICCs                 #
        # Qmatrix ... (I x K ) matrix of Q matrix entries   #
        # theta ... L  vector with trait discretization     #
        # quant.list ... vector of quantiles for attribute  #
        #                response curves                    #
        # b     ... item difficulty                         #
        # a     ... item discrimination                     #
        # c     ... guessing parameter                      #
		# theta ... grid of theta values					#
		#####################################################
		# generate sequence for display
        display.separate <- paste( rep("." , each=80 ) , collapse="" )
        # display progress
        cat( display.separate , "\n" )
        cat( "LSDM -- Least Squares Distance Method of Cognitive Validation \n")
        cat("Reference: Dimitrov, D. (2007) Applied Psychological Measurement, 31, 367-387.\n")
        cat( display.separate , "\n" ) 
        if (! is.null(b) ){ 
            eins <- rep(1, length(theta) )
            data <- outer(c,eins) + ( 1 - outer(c,eins) )* 
                        plogis(  outer( a , eins ) * ( outer( rep(1,nrow(Qmatrix)) , theta ) - 
						outer( b , eins ) ) )   
                    }
		Qmatrix <- as.matrix(Qmatrix)
		# print Q matrix
		cat("\nQmatrix\n\n")
		cmax <- apply( Qmatrix , 2 , max )
		Qmatrix <- Qmatrix / outer( rep(1,nrow(Qmatrix)) , cmax )
		print(Qmatrix) ; cat("\n")
		d1 <- det( t(Qmatrix)%*% Qmatrix )
		# warning for singular Q matrices
		if (abs(d1) < 10^(-8)){ stop("You inputted a singular Q matrix. LSDM cannot be computed.\n") }	
        est.icc <- T
        I <- nrow(data)
        L <- ncol(data)
        K <- ncol(Qmatrix)
        if ( is.null( rownames(Qmatrix) ) ){ rownames(Qmatrix) <- rownames(data) }
        # log probability functions
		data1 <- data
        logdata <- log( data1 + .001)
        # estimate item parameter and item quantiles
        cat("Estimation of Item Parameters \n") ; flush.console()
        icc.pars <- est.logist.quant( probcurves = data , theta = theta , 
					quantiles = quant.list , est.icc = est.icc)
        cat( display.separate , "\n" )
		#******************************************
        # Estimate attribute response curves
        cat("Estimation of Attribute Parameters \n") ; flush.console()
        ui <-  - diag( K )
		# fit function for every L (for every theta point)
        log.arc0 <- sapply( 1:L , FUN = function(tt){
            # unrestricted linear model
            mod1.tt <- stats::lm( logdata[, tt ] ~ 0 + as.matrix(Qmatrix ) )
			# include disjunctive version here
			# P = A1 * A2
			# log(P) = log(A1) + log(A2) for every tt
			# including weights leads to
			# log(P) = w1 * log(A1) + w2 * log(A2)
            # restricted linear model
#            mod2.tt <- ic.infer::orlm( mod1.tt , index = 1:K , ui )
			mod2.tt <- ic.infer::orlm.lm( mod1.tt , index = 1:K , ui )			
			mod2.tt$b.restr
                } )	
		#*******************************************
		# estimate "ordinary" LLTM
#        lltm.res1 <- lm( as.numeric(icc.pars$b.1PL) ~ as.matrix(Qmatrix ) )
        lltm.res1 <- stats::lm( as.numeric(icc.pars$b.1PL) ~ 0 + as.matrix(Qmatrix ) )
        slltm.res1 <- summary(lltm.res1)
        cat( display.separate , "\n" )
        # exponentiate attribute response curve
        arc0 <- exp(log.arc0)
        # calculate Rasch data predicted by LLTM
        data.lltm <- outer( lltm.res1$fitted , theta , FUN = function(x1,x2){ plogis( x2 - x1  )  } )
        rownames(arc0) <- colnames(Qmatrix)
        # estimate attribute parameter and attribute quantiles
        arc0.pars <- est.logist.quant( probcurves = arc0 , theta = theta , 
				quantiles = quant.list , est.icc = est.icc )
 #       arc0.pars$eta.LLTM <- coef(lltm.res1)[-1]
        arc0.pars$eta.LLTM <- coef(lltm.res1)
#        arc0.pars$se.LLTM <- slltm.res1[[4]][ -1,2]
#        arc0.pars$pval.LLTM <- slltm.res1[[4]][ -1,4]        
        arc0.pars$se.LLTM <- slltm.res1[[4]][ ,2]
        arc0.pars$pval.LLTM <- slltm.res1[[4]][ ,4]        
        # estimate weight parameters in LSDM
#        logpistar <- rep( 0 , I )
		W <- matrix( NA , nrow=I , ncol=K )
        for (ii in 1:I){
#            index.ii <- which( Qmatrix[ii, ] == 1 )
			index.ii <- which( Qmatrix[ii,] > 0 )
            L.ii <- length(index.ii)
            x.ii <- t( log.arc0[ index.ii , ] )
			x.ii <- matrix( x.ii , ncol = length(index.ii) )
#			m1 <- outer( rep( 1 , nrow(x.ii)) , Qmatrix[ii,index.ii] )
#			x.ii <- x.ii * m1		
#            if (nrow(x.ii) == 1 ){ x.ii <- x.ii[1,] } else { x.ii <- rowSums( x.ii ) }
            y.ii <- as.numeric(logdata[ii,])
            mod1.ii <- stats::lm( y.ii ~ 0 + x.ii )	
#            ui <- diag( c( - 1 ,  -1 ) )
#            ci <-  c( 0 , -1 )
#            mod2.ii <- orlm( mod1.ii , index = 1:2 , ui , ci)
#            cm2.ii <- coef(mod2.ii)
#            logpistar[ii] <- cm2.ii[1]
            W[ii,index.ii] <- coef( mod1.ii )
					}
#        pistar <- exp( logpistar )
        #***************************
        # evaluate goodness of fit
        data0.fitted <- exp( as.matrix( Qmatrix ) %*% log.arc0 )
        rownames(data.lltm) <- rownames(data0.fitted) <- rownames(data)
        # MAD for original model (Dimitrov)
        mad0 <- rowMeans( abs( data - data0.fitted ) )
        md0 <- rowMeans( ( data - data0.fitted ) )        
        mm0 <- mean(mad0)
        mad.lltm <- rowMeans( abs( data - data.lltm ) )
        md.lltm <- rowMeans( ( data - data.lltm ) )
        mm.lltm <- mean(mad.lltm)
        # Model Fit LSDM
        cat(paste( "Model Fit LSDM   -  Mean MAD:" , formatC( round( mm0 , 3 ),digits=3 , width=6) , 
                        "    Median MAD:" , formatC( round( median(mad0) , 3 ),digits=3 , width=6)     , "\n") )
        cat(paste( "Model Fit LLTM   -  Mean MAD:" , formatC( round( mm.lltm , 3 ),digits=3, width=6) , 
                    "    Median MAD:" , formatC( round( median(mad.lltm) , 3 ),digits=3 , width=6) ,
                    "   R^2=" , format( round( slltm.res1$r.squared , 3 ),digits=3) ,   "\n") )
        item <- data.frame( "N.skills" = rowSums( Qmatrix ) , 
                            "mad.lsdm" = mad0 , 
                            "mad.lltm" = mad.lltm , 
                            "md.lsdm" = md0 , 
                            "md.lltm" = md.lltm , 
#                            "pistar" = pistar ,
#                            "W" = W , 
                            icc.pars ) 
	    colnames(W) <- colnames(Qmatrix)
	    rownames(W) <- rownames(Qmatrix)		
        res <- list(   "mean.mad.lsdm0" = mm0 , 
                        "mean.mad.lltm" = mm.lltm , 
                        "attr.curves" = arc0 ,
                        "attr.pars" = arc0.pars  ,
                        "data.fitted" = data0.fitted , 
                        "theta" = theta ,
                        "item" = item ,
                        "data" = data , "Qmatrix" = Qmatrix ,
                        "lltm" = lltm.res1 , "W"=W )
        class(res) <- "lsdm"
        return( res )
        }

#.............................
# summary for LSDM function   
summary.lsdm <- function( object , ... ){
		lsdmobj <- object
        # generate sequence for display
        display.separate <- paste( rep("." , each=80 ) , collapse="" )
        # display progress
        cat( display.separate , "\n" )
        cat( "LSDM -- Least Squares Distance Method of Cognitive Validation \n")
        cat("Reference: Dimitrov, D. (2007) Applied Psychological Measurement, 31, 367-387.\n")
        cat( display.separate , "\n" ) ; flush.console()
        cat("\nModel Fit\n\n")
        cat(paste( "Model Fit LSDM   -  Mean MAD:" , formatC( round( lsdmobj$mean.mad.lsdm0 , 3 ),digits=3 , width=6) , 
                        "    Median MAD:" , formatC( round( median(lsdmobj$item$mad.lsdm) , 3 ),digits=3 , width=6)     , "\n") )
        cat(paste( "Model Fit LLTM   -  Mean MAD:" , formatC( round( lsdmobj$mean.mad.lltm , 3 ),digits=3, width=6) , 
                    "    Median MAD:" , formatC( round( median(lsdmobj$item$mad.lltm) , 3 ),digits=3 , width=6) ,        
                    "   R^2=" , format( round( summary(lsdmobj$lltm)$r.squared , 3 ),digits=3) ,   "\n") )
        cat( display.separate , "\n" )
        cat("\nAttribute Parameters\n\n")
        dfr.a <- data.frame( "N.Items" = colSums(lsdmobj$Qmatrix) , round( lsdmobj$attr.pars , 3 )  )
        elim.a <- union( grep("Q" , colnames(dfr.a) ) , grep( "sigma" , colnames(dfr.a) ) )
        print( dfr.a[ , -elim.a ] )
        cat( display.separate , "\n" )
        cat("\nItem Parameters\n\n")
        dfr.i <- round( lsdmobj$item , 3 )

        displ.i <- sort( union(  which( colnames(dfr.i)   %in% c(  "a.2PL" , "b.2PL" , "N.Items" , "b.1PL" ) )  ,
                            grep( "mad" , colnames(dfr.i)   ) )                        
                                        )
        print( round( dfr.i[ ,  displ.i  ] , 3 ) )
        cat( display.separate , "\n" )
        cat("\nDiscrimination Parameters\n\n")
        dfr.i <- round( lsdmobj$W , 3 )
        print( dfr.i )
		
		
        }


