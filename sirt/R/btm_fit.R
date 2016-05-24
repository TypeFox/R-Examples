

############################################
# item outfit and infit statistic
btm_fit <- function( probs , dat0 , ind1 , ind2 , TP){			
			# first individual
			X_exp1 <- probs[,1]*1 + 1/2*probs[,3]
			X_var1 <- 1*probs[,1] + 1/4*probs[,3]
			X_var1 <- X_var1 - X_exp1^2
			Z_1 <- ( dat0[,3] - X_exp1 ) / sqrt( X_var1 )
	
			# second individual
			X_exp2 <- probs[,2]*1 + 1/2*probs[,3]
			X_var2 <- 1*probs[,2] + 1/4*probs[,3]
			X_var2 <- X_var2 - X_exp2^2
			Z_2 <- ( 1 - dat0[,3] - X_exp2 ) / sqrt( X_var2 )
			
			# compute outfit statistic
			out1 <- rowsum( Z_1^2 , dat0[,1] )
			N1 <- rowsum( 1+0*Z_1 , dat0[,1] )
			out2 <- rowsum( Z_1^2 , dat0[,2] )
			N2 <- rowsum( 1+0*Z_1 , dat0[,2] )
			wvar1 <- rowsum( X_var1 , dat0[,1] )
			wvar2 <- rowsum( X_var2 , dat0[,2] )
			win1 <- rowsum( X_var1*Z_1^2 , dat0[,1] )
			win2 <- rowsum( X_var2*Z_1^2 , dat0[,2] )
			
			out <- btm_combine_tables( out1 , out2 , ind1 , ind2 , TP )					
			N <- btm_combine_tables( N1 , N2 , ind1 , ind2 , TP )					
			wvar <- btm_combine_tables( wvar1 , wvar2 , ind1 , ind2 , TP )			
			win <- btm_combine_tables( win1 , win2 , ind1 , ind2 , TP )
			
			res0 <- list( outfit = out / N , 
						  infit = win / wvar )
			return(res0)
				}
#########################################################################				
btm_combine_tables <- function( win1 , win2 , ind1 , ind2 , TP ){
			win <- rep( 0 , TP )
			win[ind1] <- win1[  , 1] 
			win[ind2] <- win[ind2] + win2[  , 1] 
			return(win)
					}	
#####################################################################					