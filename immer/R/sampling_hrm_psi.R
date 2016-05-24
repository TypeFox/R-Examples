
###########################################################
# sampling psi parameters
sampling_hrm_psi <- function( dat , dat_ind , maxK , R , rater , pid , phi , psi ,
		prior , MHprop , I , xi , useRcpp , est_settings  ){
		eps <- 1E-20
		est.psi <- est_settings$est.psi
		MHprop$refresh_count$psi <- MHprop$refresh_count$psi + 1			
		minpsi <- .001
		N <- nrow(xi)
		if ( est.psi == "r"){ 
				ratio1 <- log( rep(1,R) )
							}			
		if ( est.psi == "e"){ 
				ratio1 <- log( 1 )
							}		
		if ( est.psi != "n" ){
			for (ii in 1:I){
			# ii <- 1
			if ( est.psi=="a"){
				psi_new <- psi_old <- psi
				psi_new[ii,] <- rtrnorm( N=R , mean=psi[ii,] , sd=MHprop$SD$psi[ii,] , 
									lower = rep(0,R) , upper=rep(Inf,R) )
							}			
			  if ( ( est.psi=="r" ) & ( ii == 1 ) ){
				psi_new <- psi_old <- psi
				psi_new[1:I,] <- matrix( rtrnorm( N=R , mean=psi[1,] , sd=MHprop$SD$psi[1,] , 
									lower = rep(0,R) , upper=rep(Inf,R) ) , nrow=I , ncol=R , byrow=TRUE )
								}	

			  if ( ( est.psi=="e" ) & ( ii == 1 ) ){
				psi_new <- psi_old <- psi
				psi_new[1:I,1:R] <- rtrnorm( N=1 , mean=psi[1,1] , sd=MHprop$SD$psi[1,1] , 
									lower = rep(0,1) , upper=rep(Inf,1) )
								}								
								
			  if ( ( est.psi=="i" )  ){
				psi_new <- psi_old <- psi
				psi_new[ii,] <-  rep( rtrnorm( N=1, mean=psi[ii,1] , sd=MHprop$SD$psi[ii,1] , 
									lower = 0 , upper= Inf ) , R )
								}	
								
				if ( ( est.psi=="a" ) | ( ii == 1 ) ){
					p_new <- stats::dnorm( psi_new[ii,] , mean = prior$psi$M[ii,] , sd = prior$psi$SD[ii,] ) 
					p_old <- stats::dnorm( psi_old[ii,] , mean = prior$psi$M[ii,] , sd = prior$psi$SD[ii,] ) 
									} else {
					p_new <- p_old <- rep(1,R)
									}
				
				
				ll_new <- probs_hrm( x= dat[,ii] , xi=xi[ pid , ii ] , phi = phi[ ii , rater ] , 
							   psi = psi_new[ii,rater ] , K=maxK[ii] , x_ind = dat_ind[,ii] , useRcpp)
				ll_new <- rowsum( log( ll_new + eps ) , rater )[,1]
				ll_old <- probs_hrm( x= dat[,ii] , xi=xi[ pid , ii ] , phi = phi[ ii , rater ] , 
							   psi = psi_old[ii,rater ] , K=maxK[ii] , x_ind = dat_ind[,ii] , useRcpp)
				ll_old <- rowsum( log( ll_old + eps ) , rater )[,1]
				ratio <- p_new * exp( ll_new - ll_old ) / ( p_old  )				
	
				if ( est.psi %in% c("r") ){ 
					ratio1 <- ratio1 + log(ratio) 
							}				

				if ( est.psi %in% c("e") ){ 
					ratio1 <- ratio1 + sum( log(ratio) )
							}				
							
							
				if ( est.psi == "a"){
					for (rr in 1:R){
					  if( is.na(ratio[rr] ) ){ ratio[rr] <- 0 }
						if ( ratio[rr] > stats::runif(1) ){
									MHprop$accept$psi[ii,rr] <- MHprop$accept$psi[ii,rr] + 1 
									psi[ii,rr] <- psi_new[ii,rr]
													 }
									}  # end rr
							}
				if ( est.psi == "i"){	
				  ratio <- exp(sum( log(ratio) ))
					  if( is.na(ratio ) ){ ratio <- 0 }
						if ( ratio > stats::runif(1) ){
									MHprop$accept$psi[ii,1:R] <- MHprop$accept$psi[ii,1:R] + 1 
									psi[ii,1:R] <- psi_new[ii,1:R]
													 }
							}							

								
						}  # end ii
				#------------		
				if ( est.psi == "r"){
					for (rr in 1:R){
					  ratio[rr] <- exp(ratio1[rr])
					  if( is.na(ratio[rr] ) ){ ratio[rr] <- 0 }
						if ( ratio[rr] > stats::runif(1) ){
									MHprop$accept$psi[ii,rr] <- MHprop$accept$psi[ii,rr] + 1 
									psi[ii,rr] <- psi_new[ii,rr]
													 }
									}  # end rr
							} 

				if ( est.psi == "e"){			
					ratio <- exp(ratio1)
					  if( is.na(ratio ) ){ ratio <- 0 }
						if ( ratio > stats::runif(1) ){
									MHprop$accept$psi[1:I,1:R] <- MHprop$accept$psi[1:I,1:R] + 1 
									psi[1:I,1:R] <- psi_new[1:I,1:R]
													 }
							}
							
				} # end est.psi != "n"

		res <- list( psi = psi , MHprop = MHprop )
        return(res)				
					}