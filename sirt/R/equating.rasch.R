
#***************************************************************************************************
# Equating (Linking) in the Rasch model
equating.rasch <- function( x , y , theta = seq( -4 , 4 , len=100) ,
		alpha1 = 0 , alpha2 = 0 ){
    # INPUT:
    # x ... data frame: 1st column Item labels group 1, 2nd column: item difficulties group 1
    # y ... data frame: 1st column Item labels group 2, 2nd column: item difficulties group 2
    # theta ... theta values where the test characteristic curves are evaluated
	# alpha ... parameters for gernealized logistic distribution
    #****************************
    # Data preparation 
    x[,1] <- gsub( " " , "" , paste( x[,1] ) )
    y[,1] <- gsub( " " , "" , paste( y[,1] ) )
    b.xy <- data.frame( merge( x , y , by.x = 1 , by.y = 1 ) )
    colnames(b.xy) <- c("item" , "Itempar.Gr1" , "Itempar.Gr2" )
    b.xy <- stats::na.omit( b.xy )
    #    b.xy <- b.xy[ rowSums( is.na( b.xy[,2:3] ) ) < 2 , ]
    # mean-mean method 
        B.mm <- mean(b.xy[,3]) - mean(b.xy[,2])
	g1 <- .prob.raschtype.genlogis( theta = theta , b = b.xy[,2]  , alpha1 = 0 , alpha2 = 0 )
    # Haebara function
        ha <- function(B){
            sum( ( 1 / ( 1+  exp( outer( theta , b.xy[,2] , "-" )  ) ) - 1/ ( 1 + exp( outer( theta , b.xy[,3] - B , "-" ) ) ) )^2 )
            }
		if ( abs( alpha1) + abs(alpha2) > 0 ){ 
			ha <- function(B){ 
				sum( ( .prob.raschtype.genlogis( theta = theta , b = b.xy[,2]  , alpha1 = alpha1 , alpha2 = alpha2 ) -
						.prob.raschtype.genlogis( theta = theta , b = b.xy[,3] - B , alpha1 = alpha1 , alpha2 = alpha2 ) )^2 )
											}
								}
        B.ha <- stats::optimize(  ha , c(-4,4) )$minimum
    # Stocking and Lord Approach
        sl <- function(B){
            sum( ( rowSums( 1 / ( 1+  exp( outer( theta , b.xy[,2] , "-" )  ) )  - 
							1/ ( 1 + exp( outer( theta , b.xy[,3] - B , "-" ) ) ) ) 
											)^2 )
							}
		if ( abs( alpha1) + abs(alpha2) > 0 ){ 
			sl <- function(B){ 
            sum( ( rowSums( .prob.raschtype.genlogis( theta = theta , b = b.xy[,2]  , alpha1 = alpha1 , alpha2 = alpha2 )   - 
							.prob.raschtype.genlogis( theta = theta , b = b.xy[,3] - B , alpha1 = alpha1 , alpha2 = alpha2 ) 
											) )^2 )
											}
								}						
        B.sl <- stats::optimize(  sl , c(-4,4) )$minimum
    # all parameter estimates    
    B.est <- data.frame( B.mm , B.ha , B.sl )
    colnames(B.est) <- c("Mean.Mean" , "Haebara" , "Stocking.Lord")
    # Transformation of item parameters (according to Stocking-Lord)
    b.xy$TransfItempar.Gr1 <- b.xy[,2] + B.est[1,"Stocking.Lord"]
    x[,2] <- x[,2] + B.est[1,"Stocking.Lord"]
    # transformed parameters
    transf.par <- merge( x , y , 1 , 1 , all=T )
    colnames(transf.par) <- c("item" , "TransfItempar.Gr1" , "Itempar.Gr2"  )
    transf.par <- transf.par[ order( paste(transf.par$item ) ) , ]
    # calculate variance and linking error
    des <- data.frame( "N.Items" = nrow(b.xy) , "SD" = stats::sd( b.xy$TransfItempar.Gr1 - b.xy$Itempar.Gr2 ) )
    des$Var <- des$SD^2
    des$linkerror  <- as.vector( sqrt( des["SD"]^2 / des["N.Items"] ) )[1,1]
    # OUTPUT:
    # B.est ... estimated shift parameter according to the three methods    
    # anchor ... original and transformed item parameters of anchor items   
    # transf.par   ... transformed item paramters
    return( list( "B.est" = B.est , "descriptives" = des , 
                 "anchor" = b.xy[ , c(1,2,4,3)] , "transf.par" = transf.par 
                 ) )
        }
#***************************************************************************************************


##NS export(equating.rasch.jackknife)
#-----------------------------------------------------------------------------------------------
# Computation of linking error using Jackknife
equating.rasch.jackknife <- function( pars.data , display = TRUE , se.linkerror = FALSE ,
			alpha1=0 , alpha2 = 0 ){ 
        #.......................................................
		# using generalized logistic model
		#--------------------------------------------------------
        # first column: jackknife unit
        # second column: parameter first scale
        # third column: parameter second scale
		# fourth column: item ?????
        pars.data <- as.data.frame( stats::na.omit( pars.data ) )
        itemunits <- unique( pars.data[,1] )
        N.units <- length( itemunits )
        N.items <- nrow( pars.data )
        pars.data[,4] <- paste("I" , 1:N.items,sep="")
        # display
        if (display == TRUE){ cat( paste( "Jackknife Equating Procedure (Stocking-Lord)\n" , 
						N.items , " Items in ", N.units , " Units\n" , sep="") ) }
        # equating without jackknife
        mod1 <- equating.rasch( pars.data[ , c( 4 , 2) ] , pars.data[ , c(4, 3) ] )
        res1 <- data.frame( "unit" = itemunits , "shift" = 0 , "SD" = 0 , "linkerror" = 0)


        for (nn in 1:N.units){
                    #            nn <- 1
            pars.data1 <- pars.data[ pars.data[,1] != itemunits[nn] , ]
            mod.nn <- equating.rasch( pars.data1[ , c( 4 , 2) ] , pars.data1[ , c(4, 3) ] )
            res1[ nn , "shift" ] <- mod.nn$B.est$Stocking.Lord
            res1[ nn , "SD" ] <- mod.nn$descriptives$SD

            # Jackknife of the linking error
            if (se.linkerror == TRUE){ 
                    itemunits.nn <- itemunits[ - nn ]
                    l1 <- NULL
                        for (ii in itemunits.nn){
                                #       ii <- itemunits.nn[1]
                                pars.data1.ii <- pars.data1[ paste(pars.data1[,1]) != ii , ]           
                                mod.ii <- equating.rasch( pars.data1.ii[ , c( 4 , 2) ] , pars.data1.ii[ , c(4, 3) ] ,
											alpha1=alpha1 , alpha2=alpha2)
                                l1 <- c(l1 , mod.ii$B.est$Stocking.Lord )
                                    }
                    res1[ nn , "linkerror"] <-  sqrt( ( N.units - 2 ) / ( N.units -1 ) * sum( ( l1 - res1[ nn , "shift" ]  )^2 ) )
                                   }
            # display progress
                if (display == TRUE){
                    cat( paste( nn , " " , sep="" ) ) ; 
					utils::flush.console()
                    if ( nn %% 10 == 0){ cat("\n") }
                                }
                    }
        cat("\n")
        linkerror <- sqrt( ( N.units - 1 ) / N.units * sum( ( res1[,2] - mod1$B.est$Stocking.Lord )^2 ) )
        se.sd <- sqrt( ( N.units - 1 ) / N.units * sum( ( res1[,3] - mod1$descriptives$SD )^2 ) )
        if (se.linkerror == TRUE){ 
                        se.linkerror <- sqrt( ( N.units - 1 ) / N.units * sum( ( res1[,4] - linkerror )^2 ) )
                                } else { se.linkerror <- NA }
       
        res <- data.frame(  "N.items" = N.items ,   "N.units" = N.units , 
                            "shift" = mod1$B.est$Stocking.Lord , "SD" = mod1$descriptives$SD ,
                                "linkerror.jackknife" = linkerror , "SE.SD.jackknife" = se.sd , 
                                "se.linkerror.jackknife" = se.linkerror )
        res <- list( "pars.data" = pars.data , "itemunits" = itemunits , "descriptives" = res )
        return( res )
        }
#.................................................................................................

