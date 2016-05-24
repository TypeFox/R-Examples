
####################################################################
immer_identifiers_relabel <- function( dat , pid , rater ){

			dat0 <- data.frame( pid , rater , dat )
			
			
			pid_unique <- unique( paste(pid) )
			rater_unique <-  unique( paste( rater ) )			
			dat0$pid <- match( pid , pid_unique )
			dat0$rater <- match( rater , rater_unique )
			
			# item parameters
			item0 <- data.frame( "item" = colnames(dat) , 
						"N_Rat" = colSums( 1 - is.na( dat ) ) ,
						"M"= colMeans( dat , na.rm=TRUE )
								)
			# item and rater combinations
			rater_pars0 <- NULL
			I <- ncol(dat)
			R <- length(rater_unique)
			for (ii in 1:I){
				dfr <- data.frame( "item" = colnames(dat)[ii] , 
				    "rater" = rater_unique , 
					"rid" = 1:R , 
					"N_Rat" = stats::aggregate( 1 - is.na(dat[ , ii  ]) , list(rater ) , sum )[,2] ,
					"M" = stats::aggregate( dat[,ii] , list(rater ) , mean , na.rm=TRUE )[,2] 
									)
				rater_pars0 <- rbind( rater_pars0 , dfr )
					}
			
			res <- list( pid = dat0$pid , rater = dat0$rater ,
							dat = dat0[ , - c(1,2) , drop=FALSE ] ,
							pid_unique = pid_unique , 
							rater_unique = rater_unique , item0 = item0 ,
							rater_pars0 = rater_pars0						
									)
			return(res)		

				}
####################################################################				