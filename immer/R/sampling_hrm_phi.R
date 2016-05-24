
###########################################################
# sampling phi parameters
sampling_hrm_phi <- function( dat , dat_ind , maxK , R , rater , pid , phi , psi ,
		prior , MHprop , I , xi , useRcpp, est_settings  ){
		eps <- 1E-20
		MHprop$refresh_count$phi <- MHprop$refresh_count$phi + 1			
		est.phi <- est_settings$est.phi
		if ( est.phi == "r"){ 
				ratio1 <- log( rep(1,R) )
							}		
		if ( est.phi != "n" ){
			for (ii in 1:I){
			# ii <- 1				
			  if ( est.phi=="a"){
				phi_new <- phi_old <- phi
				phi_new[ii,] <- stats::rnorm( R , mean= phi[ii,] , sd = MHprop$SD$phi[ii,] )
								}
			  if ( ( est.phi=="r" ) & ( ii == 1 ) ){
				phi_new <- phi_old <- phi
				phi_new[1:I,] <- matrix( stats::rnorm( R , mean= phi[1,] , sd = MHprop$SD$phi[1,] ) ,
										nrow=I , ncol=R , byrow=TRUE )
								}																								
				if ( ( est.phi=="a" ) | ( ii == 1 ) ){
					p_new <- stats::dnorm( phi_new[ii,] , mean = prior$phi$M[ii,] , sd = prior$phi$SD[ii,] ) 
					p_old <- stats::dnorm( phi_old[ii,] , mean = prior$phi$M[ii,] , sd = prior$phi$SD[ii,] ) 
									} else {
					p_new <- p_old <- rep(1,R)
									}
				ll_new <- probs_hrm( x= dat[,ii] , xi=xi[ pid , ii ] , phi = phi_new[ ii , rater ] , 
							   psi = psi[ii,rater ] , K=maxK[ii] , x_ind = dat_ind[,ii] , useRcpp )
				ll_new <- rowsum( log( ll_new + eps ) , rater )[,1]
				ll_old <- probs_hrm( x= dat[,ii] , xi=xi[ pid , ii ] , phi = phi_old[ ii , rater ] , 
							   psi = psi[ii,rater ] , K=maxK[ii] , x_ind = dat_ind[,ii] , useRcpp  )
				ll_old <- rowsum( log( ll_old + eps ) , rater )[,1]
				ratio <- p_new * exp( ll_new - ll_old ) / ( p_old  )
				if ( est.phi=="r"){
					ratio1 <- log( ratio ) + ratio1
									}
			
				# estimate phi[i,r]
				if ( est.phi == "a" ){	
					for (rr in 1:R){
					  if( is.na(ratio[rr] ) ){ ratio[rr] <- 0 }
						if ( ratio[rr] > stats::runif(1) ){
									MHprop$accept$phi[ii,rr] <- MHprop$accept$phi[ii,rr] + 1 
									phi[ii,rr] <- phi_new[ii,rr]
													 }
									}  # end rr
								}  # end est.phi == "a"										
						}  # end ii
				####****
				# if est.phi == "r"
				if ( est.phi == "r" ){	
					for (rr in 1:R){
					  ratio[rr] <- exp(ratio1[rr]) 
					  if( is.na(ratio[rr] ) ){ ratio[rr] <- 0 }
						if ( ratio[rr] > stats::runif(1) ){
									MHprop$accept$phi[1:I,rr] <- MHprop$accept$phi[1:I,rr] + 1 
									phi[1:I,rr] <- phi_new[1:I,rr]
													 }
									}  # end rr
								}  # end est.phi == "r"														
						
					}
					
		res <- list( phi = phi , MHprop = MHprop )
        return(res)				
					}