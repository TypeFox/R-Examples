
#####################################################
# Hierarchical rater model Patz et al. (2002)
immer_HRM <- function( dat , pid , rater , 
			      iter ,  burnin , N.save = 3000 , prior=NULL , 
				  est.a = FALSE , est.sigma = TRUE , 
				  est.mu = FALSE , est.phi = "a" , est.psi = "a" ,  
				  MHprop=NULL  ,
				  theta_like = seq( -10 , 10 , len=30) ){
			
		time <- list( "start" = Sys.time() )	
		useRcpp <- TRUE	
		CALL <- match.call()
			
		#************************************
		# rater and person identifiers
		res0 <- immer_identifiers_relabel( dat , pid , rater )		
		rater <- res0$rater
		pid <- res0$pid		
		dat <- res0$dat		
		pid_unique <- res0$pid_unique
		rater_unique <- res0$rater_unique
		item0 <- res0$item0	
		rater_pars0 <- res0$rater_pars0
			
		#************************************
		# inits						
		# inits theta parameters
		theta <- inits_theta_1dim( dat , pid , eps=.05 )							
		N <- length(theta)
		
		# parameters settings
		 # Options are 'n' (no estimation), 'e' (set all parameters equal to each other), 
		 # 'i' (item wise estmation), 'r' (rater wise estimation) and 
		 # 'a' (all parameters are estimated independently from each other). 		
		est_settings <- list( est.a = est.a , est.sigma = est.sigma , est.mu = est.mu ,
				est.phi = est.phi , est.psi = est.psi )

		 
		# inits item parameters
		res0 <- inits_itempars( dat , prior)					
		b <- res0$b
		a <- res0$a
		K <- res0$K
		I <- res0$I
		maxK <- res0$maxK

		
		# inits rater parameters HRM
		res0 <- inits_raterpars_hrm( rater , I , est_settings )
		phi <- res0$phi
		psi <- res0$psi
		R <- res0$R
		dat <- as.matrix(dat)
		dat_ind <- 1 - is.na(dat)
		xi_ind <- 1 * ( rowsum( 1 - is.na(dat) , pid ) > 0 )
		ND <- nrow(dat)
		
		#*********************************************
		# prior parameters			
	    prior <- prior_hrm( prior , b , a , phi , est_settings )	
		sigma <- sqrt( prior$sigma2$sig02 )
		mu <- prior$mu$M
		
		#**********************************************
		# objects for saving traces	


		BB1 <- iter - burnin
		save_list <- seq( burnin+1 , iter )
		if ( N.save > BB1 ){ BB <- BB1 }
		if ( N.save < BB1 ){
			h1 <- ceiling( BB1 / N.save )			
			BB <- BB1 / h1
			save_list <- seq( burnin+1 , iter , h1)
							}
		
		psiM <- phiM <- array( NA , dim=c(I,R,BB) )	
		bM <- array(NA , dim=c(I,K,BB) )
		aM <- matrix( NA , nrow= I , ncol=BB )
		sigmaM <- rep( NA , BB )
		muM <- rep( NA , BB )
		devM <- rep(NA, BB)
		person <- data.frame( "pid" = pid_unique , "NSamp"=0 , "EAP" = 0 , 
						"SD.EAP" = 0 )
		
		#**********************************************
		# Metropolis-Hastings tuning		
		MHprop <- MHprop_hrm( MHprop , b , a , phi , theta , iter , burnin )
		
		
		#********************************************
		#  ITERATIONS
		it <- 0
		bb <- 1
		
		eps <- 1E-20
		eps11 <- 1E-7
		dat <- as.matrix(dat)
		dat_ind <- as.matrix( dat_ind )
		maxcat <- max( maxK )
		b <- as.matrix(b)
		phi <- as.matrix(phi)
		psi <- as.matrix(psi)
		xi_ind <- as.matrix( xi_ind )


		
		for ( it in seq(1,iter ) ){
#  a0 <- Sys.time()

			#**** sample xsi
			xi <- sampling_hrm_xi( dat , theta , b , a , phi , psi , K  , pid , rater ,  ND ,
						 dat_ind , N , I , maxK , useRcpp , xi_ind )
					 
						 					 
# cat("samp xsi") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1
						 
			#**** sample b
			res0 <- sampling_hrm_b( xi , xi_ind  , b , a , maxK ,  prior , MHprop , I  , theta , useRcpp )
			b <- res0$b
			MHprop <- res0$MHprop

# cat("samp b") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1			
			
			#**** sample a
			res0 <- sampling_hrm_a( xi , xi_ind  , b , a , maxK ,  prior , MHprop , I , theta )
			a <- res0$a
			MHprop <- res0$MHprop
			
# cat("samp a") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1			
			
			#**** sample phi
			res0 <- sampling_hrm_phi( dat , dat_ind , maxK , R , rater , pid , phi , psi ,
					      prior , MHprop , I , xi , useRcpp , est_settings )
			phi <- res0$phi
			MHprop <- res0$MHprop

# cat("samp phi") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1
# cat(".............." , it ,  ".....\n")			
			#**** sample psi
			res0 <- sampling_hrm_psi( dat , dat_ind , maxK , R , rater , pid , phi , psi ,
					prior , MHprop , I , xi , useRcpp , est_settings )
			psi <- res0$psi
			MHprop <- res0$MHprop

# cat("samp psi") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1
			
			#**** sample theta
			mu_theta <- rep( mu ,N)
			SD_theta <- rep( sigma ,N)	
			res0 <- sampling_hrm_theta_1dim( theta , N , I , maxK , a , b , xi , xi_ind , 
							dat , dat_ind , pid , MHprop , mu_theta , SD_theta , useRcpp , eps=1E-20 )
			theta <- res0$theta
			MHprop <- res0$MHprop
# cat("samp theta") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1

			
			# sampling of mu
			mu <- sampling_hrm_mu_1dim( theta , prior , N )
			# sampling of sigma
			sigma <-  sqrt( immer_mcmc_draw_variance( 1 , w0 =prior$sigma2$w0 , 
							sig02=prior$sigma2$sig02 , n=N , sig2=var(theta) ) )

			
			#-------------- OUTPUT
			if ( it %in% save_list ){
				bM[ ,, bb ] <- b
				aM[,bb] <- a
				phiM[,,bb] <- phi
				psiM[,,bb] <- psi
				muM[bb] <- mu
				sigmaM[bb] <- sigma						
				person$NSamp <- person$NSamp + 1
				person$EAP <- person$EAP + theta
				person$SD.EAP <- person$SD.EAP + theta^2
				ND <- nrow(dat)
				#***
				# calculate deviance
			if (FALSE){				
				dev <- 1			
				theta0 <- theta[ pid ]				
				for (ii in 1:I){
					# ii <- 1
					K_ii <- maxK[ii]					
					pr_tot <- 0
					for (hh in seq(0,K_ii) ){
						# hh <- 1     # category hh for item ii and theta tt
						x0 <- rep(hh,ND)		
						pr1_hh <- probs_gpcm( x=x0 , theta= theta0 , b = as.numeric(b[ii,]) , a = a[ii] , 
									K = K_ii , useRcpp=FALSE)		
						pr2_x <- probs_hrm( x = dat[,ii] , xi = x0 , 
									phi= phi[ ii,rater ]  , psi = psi[ ii , rater ] ,
										K = K_ii , useRcpp=FALSE )
						pr_tot <- pr_tot + pr1_hh * pr2_x
								}
					dev <- ifelse( dat_ind[,ii] == 1 , dev*pr_tot , dev )
#					ll <- log( ll + 1E-200 )
#					dev <- dev + ll		
								}
					dev <- rowsum( log( dev + 1E-200 ) , pid )
					dev <- sum( dev[,1] )
# 					dev <- sum(dev)
					devM[bb] <- dev
					}
				bb <- bb + 1	
							}
				
			
			#------- update MH parameters
			if ( sum( MHprop$ITER_refreshing %in% it ) > 0  ){			
					MHprop <- MHprop_refresh( MHprop )
													}			
			
			if ( it %% 20 == 0 ){
				cat("Iteration " , it , "\n")
				utils::flush.console()
							}
			
			}
			
			
		#**********************************************
		# arrange output		

		# item parameters
		item <- data.frame( item0	,				    
						 "a" = rowMeans( aM ) 
							)
		for (kk in 1:K){
			item[ , paste0("b" , kk) ] <- rowMeans( bM[,kk,,drop=FALSE] )
					}
					
		# rater parameters
		rater_pars <- NULL		
		for (ii in 1:I){
#			ii <- 1
			dfr <- data.frame( 
						"phi" = apply( phiM[ii,,] , 1 , mean , na.rm=TRUE ) ,
						"psi" = apply( psiM[ii,,] , 1 , mean  , na.rm=TRUE ) 						
						)
			rater_pars <- rbind( rater_pars , dfr )			
					}			
		rater_pars <- data.frame( rater_pars0 , rater_pars )
		
		# person parameters
		person$EAP <- person$EAP / person$NSamp
		person$SD.EAP <- sqrt( ( person$SD.EAP - person$NSamp * person$EAP^2  ) / person$NSamp  )
		
		# EAP reliability
		EAP.rel <- var( person$EAP ) / ( var( person$EAP ) + mean( person$SD.EAP^2 ) )
		
		# information criteria
		ic <- list( N=N , I=I , R=R , ND = nrow(dat) , maxK=maxK , K = K )		
		time$end <- Sys.time()		
		
		# list with all traces
		traces <- list( b = bM , a = aM , phi=phiM , psi=psiM , mu=muM , 
							sigma=sigmaM , deviance = devM )
		# attr( traces , "NSamples" ) <- iter - burnin
		attr( traces , "NSamples" ) <- BB
		burnin -> attr(traces , "burnin")
	    iter -> attr( traces , "iter" )
		
		# collect traces and produce MCMC summary
		res11 <- immer_collect_traces( traces , est_settings  )
		summary.mcmcobj <- sirt::mcmc.list.descriptives( res11$mcmcobj )
		
		#******
		# extract estimated parameters
		est_pars <- list()
		est_pars$a <- item$a
		est_pars$b <- item[ , paste0("b" , seq(1,K) ) ]
		items <- paste0( item$item )
		phi <- matrix( NA , nrow=I , ncol=R)
		rownames(phi) <- items
		psi <- phi
		item_index <- match( paste(rater_pars$item) , items )
		phi[ cbind( item_index , rater_pars$rid ) ] <- rater_pars$phi
		psi[ cbind( item_index , rater_pars$rid ) ] <- rater_pars$psi
		est_pars$phi <- phi
		est_pars$psi <- psi
		est_pars$mu <- mean(muM)
		est_pars$sigma <- mean(sigmaM)		
				
		#****
		# compute likelihood and posterior
		res_ll <- loglik_HRM( dat , dat_ind , est_pars , theta_like ,
							rater , pid , maxK  )						

		#*****
		# add information criteria
		ic$dev <- -2*res_ll$ll
		ic <- immer_ic_hrm( ic , summary.mcmcobj )
			
		#*****
		# output list
		res <- list( person = person , item = item , rater_pars = rater_pars , 
						est_pars = est_pars , 
						sigma = est_pars$sigma , mu = est_pars$mu , 
						mcmcobj = res11$mcmcobj , summary.mcmcobj=summary.mcmcobj , 
						dat = dat , pid = pid , rater = rater , 
						EAP.rel = EAP.rel , ic = ic ,  
						f.yi.qk = res_ll$f.yi.qk , f.qk.yi = res_ll$f.qk.yi ,
						theta_like = theta_like , pi.k = res_ll$pi.k , 
						like = res_ll$ll , 						
						traces = traces ,
   					    MHprop = MHprop , prior = prior , est_settings = est_settings , 
						N.save = BB , 
						iter = iter , burnin=burnin , time=time , CALL = CALL)
		res$description <- "Function 'immer_HRM' | Hierarchical Rater Model (Patz et al., 2002)" 					
		class(res) <- "immer_HRM"
		return(res)				
			}
#############################################################################			