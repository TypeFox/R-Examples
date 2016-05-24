 

#--------------------------------------------------------------------------
# Reliability from a multidimensional nonlinear SEM for dichotomous data
greenyang.reliability <- function( object.tetra , nfactors){ 
    cat("Reliability Estimation Based on a Nonlinear SEM\n\n")
    cat("Green & Yang (2009, Psychometrika). Reliability of summed item scores\n") 
    cat("  using structural equation modeling: An alternative to coefficient alpha\n\n")
    mod.omega1 <- psych::omega( m = object.tetra$rho , nfactors = 1) 
    # reliability for one factor
    rel1 <- reliability.nonlinearSEM( facloadings = matrix( mod.omega1$schmid$sl[,1] , ncol=1 ) ,
                                thresh = object.tetra$tau )$omega.rel
    # reliability for f factors
    mod.omega <- psych::omega( m = object.tetra$rho , nfactors = nfactors) 
    rel1h <- reliability.nonlinearSEM( facloadings = matrix( mod.omega$schmid$sl[,1] , ncol=1 ) ,
                                thresh = object.tetra$tau )
    relf <- reliability.nonlinearSEM( facloadings = mod.omega$schmid$orthog ,
                                thresh = object.tetra$tau )
	#'''''''''
	# calculate Omega Hierarchical Asymptotic
        pthresh <- relf$pthresh
		I <- length(pthresh)
        # create matrix of multiplied facloadings (expected correlation)
        rho.exp <- matrix( 0 , I , I )
        # colnames(rho.exp) <- rownames(rho.exp) <- rownames(facloadings)
        # reliability matrix
        rel.matrix3 <- rel.matrix2 <- rel.matrix <- rho.exp   
        for (ii1 in 1:I){
            for (ii2 in 1:ii1){ 
                #        ii1 <- 2
                #        ii2 <- 4
                rho.exp[ii1,ii2] <- rel1h$rho.exp[ii1,ii2]
                rho.exp[ii2,ii1] <- rho.exp[ii1,ii2] 
                r1 <- rho.exp[ii1,ii2]
                rel.matrix[ii1,ii2] <- mvtnorm::pmvnorm( c(-Inf,-Inf) , pthresh[c(ii1,ii2)] ,
                                         corr = matrix( c( 1 , r1 , r1 , 1) ,2 ,2 ) ) - stats::pnorm( pthresh[ii1] ) * pnorm( pthresh[ii2] )
                rel.matrix[ii2,ii1] <- rel.matrix[ii1,ii2]
				# multidimensional analysis
                rho.exp[ii1,ii2] <- relf$rho.exp[ii1,ii2]
                rho.exp[ii2,ii1] <- rho.exp[ii1,ii2] 
                r1 <- rho.exp[ii1,ii2]				
                rel.matrix2[ii1,ii2] <- mvtnorm::pmvnorm( c(-Inf,-Inf) , pthresh[c(ii1,ii2)] ,
                                         corr = matrix( c( 1 , r1 , r1 , 1) ,2 ,2 ) ) - 
										 stats::pnorm( pthresh[ii1] ) * stats::pnorm( pthresh[ii2] )
                rel.matrix2[ii2,ii1] <- rel.matrix2[ii1,ii2]
				rel.matrix3[ii1,ii2] <- rel.matrix2[ii1,ii2]
				rel.matrix3[ii2,ii1] <- rel.matrix3[ii1,ii2]				
				if (ii1 == ii2 ){ 
					r1 <- 1				
					rel.matrix3[ii1,ii2] <- mvtnorm::pmvnorm( c(-Inf,-Inf) , pthresh[c(ii1,ii2)] ,
											 corr = matrix( c( 1 , r1 , r1 , 1) ,2 ,2 ) ) - 
											     stats::pnorm( pthresh[ii1] ) * stats::pnorm( pthresh[ii2] )
					rel.matrix3[ii2,ii1] <- rel.matrix3[ii1,ii2]
										}
                                }
                        }   		
        # calculation of reliability
        omega.relha <- sum( rel.matrix ) / sum( rel.matrix2 )	
       rel1h <- sum( rel.matrix ) / sum( rel.matrix3 )	
	#'''''''''
	# eigenvalue decomposition	
	eigenval.rho <- base::svd( object.tetra$rho )$d
	#'''''''''								
#	rel1h <- rel1h$omega.rel
	relf <- relf$omega.rel
    dfr <- data.frame( "coefficient" = c( "omega_1" , "omega_h" , "omega_t",
			"omega_ha","ECV" , 
				"ExplVar" , "EigenvalRatio") , 
                    "dimensions" = c(1,nfactors , nfactors,nfactors,
						nfactors , NA , NA) ,
                    "estimate" = c( rel1 , rel1h , relf , omega.relha ,
					   mod.omega$omega.lim , 
					   round( 100*eigenval.rho[1]/sum(eigenval.rho),3) ,
					   round( eigenval.rho[1]/eigenval.rho[2],3)
								)
						)
 	dfr <- dfr[ c(1,3,2,4,5,6,7) , ]
    rownames(dfr)[1] <- c("Omega Total (1D)")
    rownames(dfr)[2] <- paste("Omega Total (",nfactors,"D)",sep="")
    rownames(dfr)[3] <- paste("Omega Hierarchical (",nfactors,"D)",sep="")
    rownames(dfr)[4] <- paste("Omega Hierarchical Asymptotic (",nfactors,"D)",sep="")
    rownames(dfr)[5] <- paste("Explained Common Variance (",nfactors,"D)",sep="")
    rownames(dfr)[6] <- "Explained Variance (First Eigenvalue)"	
    rownames(dfr)[7] <- "Eigenvalue Ratio (1st to 2nd Eigenvalue)"	
	
    dfr1 <- dfr
    dfr1[,"estimate"] <- round( dfr1[,"estimate"] , 3)
	cat( paste( rep( "-" , 70) , collapse="") )
	cat("\n\n")
    print(dfr1)
	cat("\n\n") 
	cat( paste( rep( "-" , 70) , collapse="") )
	cat("\n\nOutput from Hierarchical Factor Analysis (psych package)\n\n" )
	cat( paste( nfactors , "-dimensional model\n\n",sep="") )
    print(mod.omega)
	cat("\n\n") 
	cat( paste( rep( "." , 45) , collapse="") )
	cat( paste( "\n" , 1 , "-dimensional model\n\n",sep="") )
	mod.omega1 <- psych::omega( m = object.tetra$rho , nfactors = 1) 
	print(mod.omega1)
    invisible(dfr1)
    }
#---------------------------------------------------------------------------

