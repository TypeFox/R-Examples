
############################################################
BIFIE_by_helper_pureR <- function( 
		group_values , userfct , datalistM ,
		N , vars_index , wgt_ , wgtrep , Nimp , RR , fayfac ,
		group_index , userparnames
		){
	
	G <- length(group_values)
		
	h1 <- do.call( userfct , list( "X" = datalistM[1:N, vars_index] , "w" = wgt_ ) )
	NP <- length(h1)
	
	parsM <- matrix( NA , nrow=NP*G , ncol=Nimp )
	parsrepM <- matrix( NA , nrow=NP*G , ncol= Nimp*RR) 
	sumwgtM <- matrix( NA , nrow=G , ncol = Nimp )
	ncasesM <- matrix( NA , nrow=G , ncol = Nimp )	
	
	
	cat("|")
s1 <- Sys.time()	
	
	for (ii in 1:Nimp){
		# ii <- 1  # imputed dataset
		
		cat("-"); utils::flush.console();
		dat.ii <- datalistM[ 1:N + (ii-1)*N , ]
		
		for (gg in 1:G){	
			# gg <- 1	
			ind.gg <- which( dat.ii[ , group_index ] == group_values[gg] )
			ind.gg <- stats::na.omit(ind.gg)	
			dat1 <- dat.ii[ ind.gg , vars_index ]
			w1 <- wgt_[ ind.gg ]
			sumwgtM[gg,ii] <- sum(w1)
			ncasesM[gg,ii] <- length(w1)
			wgtrep1 <- wgtrep[ ind.gg , ]	
			h1 <- do.call( userfct , list( "X" = dat1 , "w" = w1 ) )
			parsM[ 1:NP + (gg-1)*NP , ii ] <- h1
			h1 <- sapply( 1:RR , FUN = function(rr){
					do.call( userfct , list( "X" = dat1 , "w" = wgtrep1[ , rr] ) )
									} )
			parsrepM[ 1:NP + (gg-1)*NP , 1:RR + (ii-1)*RR ] <- h1
					}
			}
		cat("|\n"); utils::flush.console()	
		
	# statistical inference	
    res0 <- .Call( "bifie_comp_vcov_within" , parsM , parsrepM , fayfac , 
				RR , Nimp , package="BIFIEsurvey" )
	u_diag <- res0$u_diag
	eps <- 1E-15
	qhat <- parsM
	u_diag <- array( u_diag , dim= c(NP*G , Nimp) )
	qbar <- rowMeans(qhat)
	var_w <- rowMeans(u_diag)
	var_b <- rowSums( ( parsM - qbar )^2 / ( Nimp - 1 + eps ) )
	df <- rubin_calc_df2( B= var_b , W= var_w  , Nimp , digits=2)
	var_t <- ( 1 + 1 / Nimp) * var_b + var_w
	fmi <-  ( 1+1/Nimp) * var_b / ( var_t + eps )           
	parsL <- list( pars = qbar , pars_se = sqrt( var_t ) ,
			pars_varWithin = var_w , pars_varBetween = var_b ,
			pars_fmi = fmi , df = df)
					
	# arrange output		
	res <- list( parsM = parsM , parsrepM = parsrepM , userfct = userfct ,
					ncasesM = ncasesM , sumwgtM = sumwgtM , N = N , NP = NP ,
					WW = RR	, parsL = parsL )

	return(res)
		}
############################################################