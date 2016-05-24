
###########################################
# agreement statistics for two raters
immer_agree2 <- function( y , w = rep(1,nrow(y) ) ,
		symmetrize = FALSE , tol= c(0,1) ){		
	CALL <- match.call()
	res <- immer_unique_patterns( dat = y , w=w )
	y <- res$y
	w <- res$w
	#***
	# symmetrize frequency table
	if ( symmetrize ){
		y <- rbind( y , y[,c(2,1) ] )
		w <- c( w/2 , w/2 )
		res <- immer_unique_patterns( dat = y , w=w )
		y <- res$y
		w <- res$w    
			}

	w0 <- w / sum(w)

	#*** create frequency table
	categs <- unique( c( y[,1] , y[,2] ) )
	CC <- length(categs)

	agree_table <- matrix( 0 , nrow=CC , ncol=CC)
	for (ii in 1:CC){
	for (jj in 1:CC){
		agree_table[ii,jj] <- sum( ( y[,1] == categs[ii] ) * 
					( y[,2] == categs[jj] ) * w0 )
					}
			}
	rownames(agree_table) <- paste0( colnames(y)[1] , "_Cat" , categs)
	colnames(agree_table) <- paste0( colnames(y)[2] , "_Cat" , categs)

	#****
	# compute absolute agreement (with tolerance specified in vector)
	TT <- length(tol)
	agree_raw <- rep(NA,TT)
	names(agree_raw) <- paste0("tol" , tol)
	for (tt in 1:TT){
		# tt <- 1
		agree_raw[tt] <- sum( w0[ abs( y[,1] - y[,2] ) <= tol[tt] ]  )
					}

	#*** marginal probabilities
	marg <- matrix( 0 , nrow=3 , ncol=CC)
	rownames(marg) <- c( colnames(y) , "Mean" )
	colnames(marg) <- paste0( "Cat" , categs )
	marg[1,] <- rowSums( agree_table )
	marg[2,] <- colSums( agree_table )
	marg[3,] <- ( marg[1,] + marg[2,] ) / 2


	#*** overall agreement
	Pa <- sum( diag( agree_table ) )

	#**** chance agreement
	Pe <- c()
	# Scott's Pi
	Pe["pi"] <- sum( marg[3,]^2 )
	# Cohen's kappa
	Pe["kappa"] <- sum( marg[1,] * marg[2,] )
	# AC1 Gwet
	Pe["AC1"] <- sum( marg[3,] * ( 1 - marg[3,] ) ) / ( CC - 1 )

	#*** algorithm Aicken's alpha
	PAk <- marg[1,]
	PBk <- marg[2,]
	res0 <- agree_aicken( PAk=PAk , PBk=PBk , Pa=Pa )

	Pe["Aicken"] <- res0$Pe
	agree_stats <- ( Pa - Pe ) / ( 1 - Pe )
	agree_stats["Aicken"] <- res0$alpha
	PAH <- res0$PAH
	PBH <- res0$PBH
	# hard to classify probs
	PH <- matrix( 0 , nrow=2 , ncol=CC )
	colnames(PH) <- colnames(marg)
	rownames(PH) <- rownames(marg)[1:2]
	PH[1,] <- PAH
	PH[2,] <- PBH

	#*****
	agree_stats <- ( Pa - Pe ) / ( 1 - Pe )
		
	#-----
	# output
	res <- list( "agree_raw" = agree_raw ,
			"agree_stats" = agree_stats ,
			"agree_table" = agree_table ,
			"marg" = marg , 
			"Pe" = Pe , "Pa" = Pa , 
			"alpha" = res0$alpha	, "PH" = PH,
			"nobs" = sum(w) , "ncat" = CC , "tol" = tol , 
			"y" = y , "w" = w , 
			"CALL" = CALL
			)
	class(res) <- "immer_agree2"
	return(res)
		}
#############################################################