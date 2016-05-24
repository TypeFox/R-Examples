
#############################################################
# Multiple Choice DINA Model
# mcdina model (de la Torre, 2009)
mcdina <- function( dat , q.matrix , group =NULL , 
         itempars = "gr" , weights=NULL , skillclasses = NULL , 
		 zeroprob.skillclasses = NULL ,  reduced.skillspace=TRUE , 
		 conv.crit = 0.0001, dev.crit = .1 , maxit = 1000 , progress=TRUE ){
	# prepare data
	s1 <- Sys.time()
	cl <- match.call()
	dat <- as.matrix(dat)
	# zero/one entries. q.matrix from ordinary DINA model!!	
	res0 <- .mcdina.prepare.qmatrix( dat , q.matrix )

	dat <- res0$dat
	q.matrix0 <- q.matrix <- res0$q.matrix

	# handle polytomous attributes
	res1 <- mcdina.modify.qmatrix( q.matrix , skillclasses)
	q.matrix <- res1$q.matrix
	q.matrix0 <- res1$q.matrix0
	skillclasses <- res1$skillclasses
	skillclasses0 <- res1$skillclasses0
	maxmaxattr <- res1$maxmaxattr
	
    dat0 <- dat
	dat.resp <- 1* ( 1 - is.na(dat) )
	dat[ dat.resp == 0 ] <- 1
	dat_ <- dat - 1 	
	eps <- 10^(-10)
	I <- ncol(dat)	# number of items
	CC <- max( q.matrix[,2] )	# maximal number of categories
	K <- ncol(q.matrix)-2		# number of skills
	if (K<=3 ){ reduced.skillspace <- FALSE }
	# group identifier
	if ( is.null(group) ){ group <- rep(1,nrow(dat))}
	group0 <- group
	group0_unique <- sort( unique( group ) )
	group <- match( group , group0_unique )
	group <- group - 1 
	G <- length( unique(group) )
	# weights
	if ( is.null(weights) ){ weights <- rep(1,nrow(dat))}
	
	# define skill classes
	if ( is.null(skillclasses) ){
		skillclasses <- as.matrix( expand.grid( 
				as.data.frame( rbind( rep(0,K) , rep(1,K)  ) ) ) )
								}
	classes <- .matrixstring( matr=skillclasses , string="P" )
	rownames(skillclasses) <- classes
	TP <- nrow(skillclasses)
	# define specification of estimation of item parameters
	if ( mean( itempars == "gr" ) == 1 ){ 
			itempars <- rep( "gr" , I ) 
			}
	if ( ( mean( itempars == "gr" ) < 1 ) & ( length(itempars) != I ) ){ 
			itempars <- rep( "jo" , I ) 
			}
		
	# prepare latent responses
	res <- .mcdina.prep.test.latent.response( q.matrix , K , TP , skillclasses , classes )

	lc <- res$lc
	lr <- res$lr		
	itemstat <- res$itemstat
	itemstat$G <- G
	itemstat$partype <- itempars
	itemstat$N.pars <- itemstat$N.lr * (itemstat$N.cat - 1 )
	itemstat$N.pars <- ifelse( itemstat$partype == "gr" , 
			itemstat$N.pars*itemstat$G , itemstat$N.pars )
	# list of lr
	lc_list <- lr_list <- list(1:I)
	for (ii in 1:I){
		lr_list[[ii]] <- lr[ lr$item == ii , ]
		lc_list[[ii]] <- lc[ lc$item == ii , ]
					}
	# delta parameter inits
	res <- .mcdina.init.delta( lc , lr )	
	delta_ideal <- res$delta_ideal
	delta0 <- res$delta
	
	# delta parameters
	delta <- array( 0 ,  dim=c(I,CC,CC,G) )
	for (gg in 1:G){
		delta[,,,gg] <- delta0
				}
	# init probabilities
	probs <- array( 0 , dim=c(I , CC , TP , G ) )				
				
	# init latent class distribution
	pi.k <- rep( 1 / TP , TP )
	pi.k <- matrix( pi.k , nrow=TP , ncol=G )

	# counts latent responses
	lr_counts <- array(0 , dim=c(I,CC,G) )

	#*****************************
	# define reduced skillspace
	Z <- Z.skillspace <- NULL
	if ( reduced.skillspace ){
		A <- skillclasses
		attr.patt <- skillclasses
		maxAttr <- 1
		# combinations
		kombis <- utils::combn( K , 2 )	
		KK <- ncol(kombis)
		B <- NULL
		for (kk in 1:KK){
			B <- cbind( B , attr.patt[ , kombis[1,kk] ] * attr.patt[ , kombis[2,kk] ] )
					}
		 Z <- cbind( 1 , A , B )
		 ncolZ <- ncol(Z)
     	v1 <- c("Int" ,  paste("A",1:K , sep="") ) 		 
		v1 <- c(v1,apply( kombis , 2 , FUN = function(ll){ 
			paste( paste( "A" , ll , sep="") , collapse="_" ) } ))
		colnames(Z) <- v1	
	
		m1 <- which( maxAttr > 1 )
		if ( max(maxAttr) > 1 ){
		   Z1 <- Z[ , m1 , drop=FALSE ]^2
		   colnames(Z1) <- paste0( colnames(q.matrix)[m1] , "*2")
		   Z <- cbind( Z , Z1 )
						}
		if ( ! is.null(Z.skillspace) ){ 
				Z <- Z.skillspace
						}
		# check for equal columns
		Z <- Z[ , ! duplicated( t(Z) ) ]
		ncolZ <- ncol(Z)		 
			}	

	
	iter <- dev <- 0	
	max.par.change <- 1000
	devchange <- 100
	# display for progress
	disp <- "...........................................................\n"	

	
	#****************************************
	#************ begin algorithm ***********
    while ( ( iter < maxit ) & 
				( ( max.par.change > conv.crit ) | ( devchange > dev.crit  ) )
					){

# z0 <- Sys.time()
					
		#--- (0) collect old parameters
		dev0 <- dev
		delta0 <- delta

		#--- (1) calculate probabilities	
		for (gg in 1:G){ # gg <- 1
			for (ii in 1:I){    # ii <- 1
				lr.ii <- lr_list[[ii]]
				probs[ii,,,gg] <- delta[ ii ,  , lr.ii$lr_index , gg]
						}
					}

# cat("calc probs ") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1	
					
		#--- (2) calculate likelihood
		# probs_pcm_groups_Cpp <- 
		# function( dat_ , dat_resp_,  group_ , probs_,  CC_ ,  TP_ ){

		probs_ <- as.matrix( array( probs , dim=c(I , CC*TP*G) ) )			
		f.yi.qk <- probs_pcm_groups_Cpp( dat_=dat_ , dat_resp_=dat.resp ,  group_ = group  ,
					 probs_ = probs_ ,  CC_ = CC,  TP_ =TP )$fyiqk

# cat("calc like ") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1	
					 
		#--- (3) calculate posterior and expected counts
		# calccounts_pcm_groups_Cpp <- 
		# function( dat_,  dat_resp_,  group_, fyiqk_,  pik_,  CC_,  weights_ )

		res1 <- calccounts_pcm_groups_Cpp( dat_ = dat_ ,  dat_resp_ = dat.resp ,  group_=group , 
					fyiqk_ = f.yi.qk ,  pik_ = pi.k ,  CC_ =CC ,  weights_ =weights )
		n.ik <- array( res1$nik , dim = c( I , CC , TP , G ) )
		count_pik <- res1$count_pik	
		for (gg in 1:G){		
			pi.k[,gg] <- count_pik[,gg] / sum( count_pik[,gg] )
				}
		# set some probabilities of skill classes to zero
		if ( ! is.null(zeroprob.skillclasses ) ){
			pi.k[ zeroprob.skillclasses , ] <- 0
					}				
		LL <- res1$LL
		dev <- -2*LL
		f.qk.yi <- res1$fqkyi
# cat("calc post ") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1	
		
		#--- (4) log-linear smoothing of skill class distribution
		
		if (reduced.skillspace){
			for (gg in 1:G){
				ntheta <- pi.k[,gg]
				res <- gdina.reduced.skillspace( ntheta , Z , 
						  reduced.skillspace.method= 2 )		
				pi.k[,gg] <- res$pred.ntheta
						}
					}
# cat("calc smoothing distribution") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1						
					
		#--- (5) calculate updated item parameters
		res1 <- mcdina.est.item( n.ik , lr_list , lc_list , delta , I , G , eps ,
				itemstat , itempars , lr_counts)
		delta <- res1$delta
		lr_counts <- res1$lr_counts

# cat("calc item parameters") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1						
		
		#--- (11) convergence
		max.par.change <- max( abs( delta - delta0 ) )
		devchange <- abs( dev- dev0)
		iter <- iter + 1	
				
		#--- (99) display progress
		if (progress) {
			cat(disp)	
			cat("Iteration" , iter , "   " , paste( Sys.time() ) , "\n" )	   	
			cat( "Deviance = "  , round( dev , 5 ) )
			g11 <-  - ( dev - dev0 )
				if (iter >1){ 
					cat(" | Deviance change = " , round( -(dev-dev0) , 7) ) 				
				if (g11 < 0 ){ cat( "\n**** Deviances decreases! Check for nonconvergence.   ****\n") 
						}
					}
			cat("\n" )				
			cat("Maximum parameter change:" , round( max.par.change, 6), "\n")
			utils::flush.console() 			
				}
 	
		}
			
	#*************** end algorithm ***********	
	#*****************************************
	
	# include information criteria
	ic <- mcdina.calc.ic( dev, weights , itemstat , pi.k , G , I ,
				zeroprob.skillclasses , reduced.skillspace , Z )
		
	# include standard error
	se.delta <- mcdina.calc.se.delta( delta , n.ik , probs , lr_list , lc_list , 
				itemstat , I , G , itempars , lr_counts , CC )

	# labeling
	rownames(pi.k) <- classes	
	colnames(pi.k) <- paste0("Group." , group0_unique )
	
	# rename skill classes in case of polytomous attributes
	if (maxmaxattr > 1 ){
		skillclasses <- skillclasses0
		lc$Q <- .matrixstring(q.matrix0[,-c(1:2) ] , "Q" )
		q.matrix <- q.matrix0 
				}

	# item overview		
	item <- mcdina.collect.itempars( I , lc , itempars , itemstat , dat ,
		G , CC , delta , se.delta , group0_unique  )

	# skill pattern
	skill.patt <- mcdina.skill.patt( q.matrix , skillclasses , G , pi.k , group0_unique )

	# person classification	
	mle.class <- skillclasses[ max.col( f.yi.qk ) , ]	
	map.class <- skillclasses[ max.col( f.qk.yi ) , ]
	N11 <- nrow(mle.class)
	K11 <-  ncol(mle.class)
	K12 <- nrow(skillclasses)
	
	eap.class <- matrix( 0 , nrow= N11 , ncol= K11 )
	colnames(eap.class) <- colnames(mle.class)
	for (kk in 1:K11){
		# kk <- 4
		sckk <- matrix( skillclasses[,kk] , nrow=N11 , ncol=K12 , byrow=TRUE )
		eap.class[,kk] <- rowSums( sckk * f.qk.yi )
						}
	# output	
	res <- list( "item" = item , "posterior"=f.qk.yi , "like"=f.yi.qk , "ic"=ic , 
				"q.matrix" = q.matrix , "pik"=probs , 
				"delta"=delta , "se.delta"=se.delta ,  "itemstat" = itemstat , 
				"n.ik"=n.ik , "deviance" = dev, 
				"attribute.patt" = pi.k , "attribute.patt.splitted" = skillclasses , 
				"skill.patt" = skill.patt , 
				"MLE.class" = mle.class , "MAP.class" = map.class , "EAP.class" = eap.class , 
				"dat" = dat0 , "skillclasses"= skillclasses , "group"=group0 ,
				"lc" = lc , "lr" =lr , "iter"=iter , "itempars" = itempars , 
				"weights"=weights , "I" = nrow(dat) , "G"= G , "CC" = CC ,
				"loglike" = - dev / 2 , "AIC" = ic$AIC , "BIC" = ic$BIC , 
				"Npars" = ic$np )
	res$converged <- iter < maxit			
				
	res$control$weights <- weights
	res$control$group <- group				
				
	s2 <- Sys.time()		
	
	res$time <- list( "s1"=s1,"s2"=s2 , "timediff"=s2-s1)				
        cat("----------------------------------- \n")
        cat("Start:" , paste( s1) , "\n")
        cat("End:" , paste(s2) , "\n")
        cat("Difference:" , print(s2 -s1), "\n")
        cat("----------------------------------- \n")
	class(res) <- "mcdina"	
	res$call <- cl
	return(res)
		}
###########################################################