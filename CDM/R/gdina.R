

################################################################################
# GDINA Model
################################################################################

gdina <-
function( data, q.matrix, skillclasses=NULL , conv.crit = 0.0001, 
					dev.crit = .1 , maxit = 1000,
					linkfct = "identity" , Mj = NULL , 
					group = NULL , 
					method = "WLS" , 
					delta.init = NULL , 
					delta.fixed = NULL ,
					delta.designmatrix = NULL , 
					delta.basispar.lower = NULL , 
					delta.basispar.upper = NULL , 					
					delta.basispar.init = NULL , 
					zeroprob.skillclasses = NULL , 
					attr.prob.init = NULL , 
					reduced.skillspace=TRUE , 
					reduced.skillspace.method=2 , 
					HOGDINA = -1 , 
					Z.skillspace = NULL , 
                    weights = rep( 1, nrow( data ) ),  rule = "GDINA", 
                    progress = TRUE , 
					progress.item = FALSE , 
					increment.factor = 1.01 ,
					fac.oldxsi = 0 ,
					avoid.zeroprobs = FALSE , 
					seed = 0 , 		
					save.devmin=TRUE , calc.se = TRUE ,
					...
						){
                    
# data: a required matrix of binary response data, whereas the items are in the columns 
#       and the response pattern in the rows. NA values are allowed.
#
# q.matrix: a required binary matrix describing which attributes are required, coded by 1,
#       and which attributes are not required, coded by 0, to master the items, whereas the
#       attributes are in the columns and the items in the rows.
#
# method: WLS (using W matrix) or ULS (without a W matrix) estimation
#
# conv.crit: termination criterion of the iterations defined as the maximum change in parameter
#       estimates. Iteration ends if maximal parameter change is below this value.
#
# maxit: maximal number of iterations.
#
# zeroprob.skillclasses:  an optional vector of integers which indicates which skill classes should have
#							zero probability
#
# weights: an optional vector of weights for the response pattern. Noninteger weights allow for different
#       sampling schemes.
#
# weight.matrix: use weighting matrix in least squares estimation

# rule: an optional character string or vector of character strings specifying the model rule that is used. 
#       The character strings must be of "DINA" or "DINO". If a vector of character strings is specified, 
#       implying an itemwise condensation rule, the vector must be of length ncol(data). The default is the 
#       condensation rule "DINA" for all items.
#		See help: DINA, DINO, ACDM (=GDINA1), GDINA1, GDINA2
#		The saturated specification GDINA is the default.
#
# progress: an optional logical indicating whether the function should print the progress of iteration.

#	skillclasses <- NULL


max.increment <- .5
# increment.factor <- 1.05

if (progress){
    cat("---------------------------------------------------------------------------------\n")
		d1 <- packageDescription("CDM")
		cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n" )		
		}
    time1 <- list( "s1" = Sys.time() )
	cl <- match.call()
	
########################################################
# add item and attribute labels	if necessary
########################################################

	if ( is.null( colnames( data ) ) ){
			colnames(data) <- paste( "Item" , seq(1,ncol(data)) , sep="")
						}
	if ( is.null( colnames( q.matrix ) ) ){
			colnames(q.matrix) <- paste( "Attr" , seq(1,ncol(q.matrix)) , sep="")
						}

################################################################################
# check consistency of input (data, q.matrix, ...)                             #
################################################################################

#    clean <- check.input(data, q.matrix, conv.crit, maxit, constraint.guess,
#        constraint.slip, guess.init, slip.init, weights, rule, progress)   

#    if (is.character(clean)) return(clean)
		dat.items <- data
	
		# check of admissible rules
		admiss.rules <- c("GDINA" , "ACDM" , "DINA" , "DINO" ,
							"GDINA1" , "GDINA2" , "RRUM" )
		i1 <- which( ! ( rule %in% admiss.rules ) )
		if ( length(i1) > 0 ){
			cat("The following rules are not implemented in gdina: ")
			cat( paste( unique( rule[i1] ) , collapse= " " ) , "\n" )
			stop("Change your argument 'rule'")
				}
	# estimation of a reduced RUM model
    rrum.params <- rrum.model <- FALSE	
	if ( any( rule  == "RRUM" ) ){
		rule <- "ACDM" 
		linkfct <- "log"
		rrum.model <- TRUE
					}
	
	
#    dat.items <- clean$data; q.matrix <- clean$q.matrix; conv.crit <- clean$conv.crit;
#    maxit <- clean$maxit; constraint.guess <- clean$constraint.guess; 
#    constraint.slip <- clean$constraint.slip; guess.init <- clean$guess.init;
#    slip.init <- clean$slip.init; weights <- clean$weights; rule <- clean$rule;
#    progress <- clean$progress    

################################################################################
# model specification: DINA, DINO or itemwise specification of DINA or DINO    #
################################################################################


#****
# include specifications here
#****
	r1 <- "GDINA Model"


################################################################################
# multiple group estimation
################################################################################

	
	G <- 1
	if ( is.factor( group ) ){ group <- paste( group ) }
	if ( ! is.null( group) ){
	    group0 <- group
		groups <- sort( unique( group) )
		G <- length(groups)	
		group2 <- match( group , groups )
		groupre <- FALSE
		if ( any( group != group2 ) ){
				group <- group2
				groupre <- TRUE
								}
							}	
    group.stat <- NULL							
    if (G>1){							
		# group statistics
		a1 <- stats::aggregate( 1+0*group , list(group) , sum )
	#    a2 <- aggregate( group0 , list(group) , mean )
		a2 <- rep("",G)
		for (gg in 1:G){
			a2[gg] <- group0[ which( group == gg )[1]  ]
						}
		group.stat <- cbind( a2 , a1 )
		colnames(group.stat) <- c(  "group.orig" , "group" , "N"  )	
					}
							
				

							
###############################################################
# HOGDINA model
	if (HOGDINA >= 0){						
		reduced.skillspace <- FALSE
		theta.k <- seq( -6,6 , len=21 )
		wgt.theta <- stats::dnorm( theta.k )
		w1 <- wgt.theta / sum( wgt.theta )
		wgt.theta <- matrix( w1 , nrow=length(w1) , ncol=G)
				}
								
							
################################################################################
# display on R console                                                         #
################################################################################

    disp <- r1      
	if (progress){
		cat(disp,"\n")
		cat( " Link function:" , linkfct , "\n")
		if (G>1){ 
			cat(" Multiple group estimation with",G,"groups\n")
			if (groupre){ cat( "  Renumbered group identifier from 1 to",G,"\n") }
				}
				}
		s1 <- Sys.time()
	if (progress){		
		cat( "**", paste(s1), "\n"   )
		cat("---------------------------------------------------------------------------------\n")
		utils::flush.console()
			}

			
################################################################################
# definition of model parameters                                               # 
################################################################################

    I <- nrow(dat.items)   # number of persons
    J <- ncol(dat.items)   # number of items
    K <- ncol(q.matrix)       # number of attributes
#    L <- 2^K               # number of latent class pattern of attributes
    dat.items <- as.matrix( dat.items)
    q.matrix <- as.matrix( q.matrix)
            
	if ( length(rule) == 1){ 
			rule <- rep( rule , J )
						}


	if (HOGDINA >= 0){						
		b.attr <- a.attr <- matrix( 0 , nrow=K , ncol=G )
				}						
						
 a0 <- Sys.time()
						
################################################################################
# Initialization and missing data handling                                     #
################################################################################

    # initialize guessing and slipping parameters
    # without constraints, the default is set equal to .2 for all items
#    guess <- guess.init ; slip <- slip.init
    
    # missing data is coded by 9
    resp <- 1 - is.na(dat.items)
    dat.items[ resp == 0 ] <- 9
    
    # standardize weights such that the sum of defined weights is equal to the number of rows in the data frame
    weights <- nrow(dat.items)*weights / sum(weights )
	

#vv <- "init" ; a1 <- Sys.time() ; cat( vv , a1-a0 , "\n") ; a0 <- Sys.time()
	
################################################################################
# calculate item response patterns                                             #
################################################################################

    # string with item response patterns
#    item.patt.subj <- sapply( 1:I, FUN = function(ii){ paste( dat.items[ ii, ], collapse="" )  } )
	 item.patt.subj <- dat.items[,1]
	 for (jj in 2:J){
		item.patt.subj <- paste( item.patt.subj , dat.items[,jj] , sep="")
					}	
    # calculate frequency of each item response pattern
    item.patt <- table( item.patt.subj )
	
    # sort item response pattern according to their absolute frequencies
    six <- sort( item.patt, index.return=FALSE, decreasing=TRUE)
    # define data frame 'item.patt' with item response pattern and its frequency (weight)
    item.patt <- cbind( "pattern" = rownames(six), "freq" = as.numeric(as.vector( six ) ) )
	
    # calculate weighted frequency for each item response pattern
	if (G== 1){ 
 
		h1 <- rowsum( weights , item.patt.subj )	
		item.patt[,2] <- h1[ match( item.patt[,1] , rownames(h1) ) , 1]
							
		item.patt.freq <- as.numeric(item.patt[,2])
				}			
				
	if (G > 1){
		item.patt.freq <- matrix( 0 , nrow(item.patt) , G )
		for (gg in 1:G){ 
#			item.patt[,2] <- sapply( 1:( nrow(item.patt) ), FUN = function(kk){
#								sum( weights * ( item.patt[ kk, 1] == item.patt.subj  ) * (group == gg ) )
# 								} )  
		h1 <- rowsum( weights * (group == gg ), item.patt.subj )	
		item.patt[,2] <- h1[ match( item.patt[,1] , rownames(h1) ) , 1]
							
			item.patt.freq[,gg] <- as.numeric(item.patt[,2])
						}
			 }


# stop()			 
# vv <- "item response patterns" ; a1 <- Sys.time() ; cat( vv , a1-a0 , "\n") ; a0 <- Sys.time()

		 
################################################################################ 
# generate all attribute patterns                                              #
################################################################################


		# extract unique Q-matrix entries
		K <- ncol(q.matrix)
		q.entries <- as.list( 1:K )

		maxAttr <- rep(1,K)
		for (kk in 1:K){ 
#			q.entries[[kk]] <- sort(unique( q.matrix[,kk] ))
			q.entries[[kk]] <- sort(unique( c(0,q.matrix[,kk] )))
			maxAttr[kk] <- length( q.entries[[kk]] ) - 1
					}
					
		attr.patt <- as.matrix( expand.grid( q.entries ) )
		if ( ! is.null(skillclasses) ){   attr.patt <- skillclasses  }
		colnames(attr.patt) <- colnames(q.matrix)
		
	
		L <- nrow(attr.patt)

    # combine all attributes in an attribute pattern as a string
    attr.patt.c <- apply( attr.patt, 1, FUN = function(ll){ paste(ll,collapse="" ) } )

	# create designmatrix for reduced skill space
    if ( K < 4 | ( ! is.null( zeroprob.skillclasses ) ) | G > 1 ){ 
			reduced.skillspace <- FALSE 
				}
	if ( ! is.null(Z.skillspace) ){ 
			reduced.skillspace <- TRUE 
			Z.skillspace <- as.matrix(Z.skillspace)
			}
	Z <- NULL ; covbeta <- NULL ; beta <- NULL
	ncolZ <- nrow(attr.patt)-1


	
	if ( reduced.skillspace ){
		A <- attr.patt
		# combinations
		kombis <- combn( K , 2 )	
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
		 ncolZ <- ncol(Z)						
			}
			
# vv <- "attribute patterns" ; a1 <- Sys.time() ; cat( vv , a1-a0 , "\n") ; a0 <- Sys.time()			



################################################################################
# uniform prior distribution of all latent class patterns                      #
################################################################################

	if (G== 1){ 
		attr.prob <- rep( 1/L, L )
		if ( seed > 0 ){
		   set.seed(seed)
		   attr.prob <- stats::runif( L , 0 , 10/L )
		   attr.prob <- attr.prob / sum( attr.prob )
						}		   
			} else {
		attr.prob <- matrix( 1/L , L , G )
		if ( seed > 0 ){
		   set.seed(seed)
		   for (gg in 1:G){
			   a1 <- runif( L , 0 , 10/L )
			   attr.prob[,gg] <- a1 / sum( a1 )
							}
						}			
					}
			
					
	if ( ! is.null(attr.prob.init) ){
          attr.prob <- attr.prob.init
		  if (G==1){ 
				attr.prob <- as.vector( attr.prob) 
						}
					}

					
################################################################################
# create design matrices 
################################################################################	

	Mj.userdefined <- TRUE
	if ( is.null(Mj) ){ 
			Mj.userdefined <- FALSE
			Mj <- as.list( rep("noe" , J ) ) 
					}
	# Notation for Mj and Aj follows De La Torre (2011)
	Aj <- NULL
	Nattr.items <- rowSums(q.matrix >= 1)
	# list of necessary attributes per item
	necc.attr <- as.list( rep(NA,J) )
	# list of rows in attr.patt which correspond to attribute classes 
		# for one item
	attr.items <- NULL
	# list of indices of attribute patterns which should be
	#    aggregated for each item
	aggr.attr.patt <- NULL
	for (jj in 1:J){ 	# loop over items jj
		# jj <- 10
		nj1 <- necc.attr[[jj]] <- which( q.matrix[jj,] > 0 )
		if ( length(nj1)==0 ){ 
				stop( paste("Q matrix row " , jj , " has only zero entries\n" , sep="") ) 
						}	
		Aj1 <- Aj[[jj]] <- .create.Aj( Nattr.items[jj] )
		if ( ! Mj.userdefined ){ 
				Mj[[jj]] <- .create.Mj( Aj[[jj]] , rule = rule[jj] )	
						}
 		l1 <- as.list( 1 )
		l2 <- rep(0,L)	
		for (zz in seq(1,nrow(Aj1)) ){  # begin row zz
 			Aj1zz <- outer( rep(1,nrow(attr.patt)) , Aj1[zz,] )
			apzz <- attr.patt[ , nj1 ]
			apzz <- 1 * ( apzz >= q.matrix[ rep(jj,L) ,nj1] )
			l1[[zz]] <- which( rowMeans( apzz == Aj1zz  ) == 1)
			l2[ l1[[zz]] ] <- zz
						  }   # end row zz
		attr.items[[jj]] <- l1
		aggr.attr.patt[[jj]] <- l2		
						}	# end item jj

    #******						
	# indices for Mj
	Mj.index <- matrix( 0 , J , 6 )
	for (jj in 1:J){
			Mj.index[jj,1] <- ncol( Mj[[jj]][[1]] )	
			Mj.index[jj,4] <- nrow( Aj[[jj]])
					}
	Mj.index[,3] <- cumsum( Mj.index[,1] )
	Mj.index[,2] <- c(1,Mj.index[-J,3] + 1 )
	Mj.index[,6] <- cumsum( Mj.index[,4] )	
	Mj.index[,5] <- c(1,Mj.index[-J,6] + 1 )	
	# compute designmatrix of aggregation of pattern
	aggr.patt.designmatrix <- matrix( 0 , L , max(Mj.index[,6]) )
	for (jj in 1:J){
#		jj <- 2
		Mj.index.jj <- Mj.index[jj,]
		for (vv in seq(1,Mj.index.jj[4]) ){
			aggr.patt.designmatrix[ , Mj.index.jj[5] - 1 + vv ] <- 1  * ( aggr.attr.patt[[jj]] == vv )
								}				
				}

				
###############################################################################
# initial item parameters
###############################################################################

	delta <- delta.init
	if ( is.null( delta.init ) ){
		
		#****
		# identity link
		if (linkfct == "identity" ){ 	
			for ( jj in 1:J){
				N1jj <- ncol(Mj[[jj]][[1]])
				l1 <- rep(0,N1jj)
				if ( seed == 0 ){
					dd1 <- .2 ; dd2 <- .6 
						} else {
					dd1 <- stats::runif( 1 , 0 , .4 )
					dd2 <- stats::runif( 1 , 0 , 1 - dd1 - .1 )				
						}
				l1[1] <- dd1
				l1[2:N1jj] <- rep( dd2 / (N1jj - 1) , N1jj - 1 )
				delta[[jj]] <- l1
							}
						}
		#*****
		# logit link
		if (linkfct == "logit" ){ 	
			for ( jj in 1:J){
				N1jj <- ncol(Mj[[jj]][[1]])
				l1 <- rep(0,N1jj)
				if ( seed == 0 ){
					dd1 <- -1 ; dd2 <- 1
						} else {			
					dd1 <- stats::runif( 1 , -2 , 0 )					
					dd2 <- stats::runif( 1 , 0 , 2 )					
						}
				l1[1] <- dd1
				l1[N1jj] <- dd2
				delta[[jj]] <- l1
				
							}
						}
		#*****
		# log link
		if (linkfct == "log" ){ 	
			for ( jj in 1:J){
				N1jj <- ncol(Mj[[jj]][[1]])
				l1 <- rep(0,N1jj)
				if ( seed == 0 ){
					dd1 <- -1.5 ; dd2 <- .75
						} else {
					dd1 <- stats::runif( 1 , -3 , -1 )					
					dd2 <- stats::runif( 1 , .25 , 1 )
						}
				l1[1] <- dd1
				l1[N1jj] <- dd2
				delta[[jj]] <- l1
							}
						}	
				}
	###########################
	# import inits delta basis parameter
	if ( ! is.null( delta.basispar.init ) ){
		u.delta <- delta.designmatrix %*% delta.basispar.init
		for (jj in 1:J){
			delta[[jj]] <- u.delta[ seq( Mj.index[jj,2] , Mj.index[jj,3] ) , 1]
						}
					}		
    #----------------------------------------
	# compute inverse matrices for least squares estimation
	invM.list <- list( 1:J )
	for (jj in 1:J){
		Mjjj <- Mj[[jj]][[1]]
#		invM.list[[jj]] <- solve( t(Mjjj) %*% Mjjj	)
		invM.list[[jj]] <- solve( crossprod(Mjjj))
				}
	
		
	
				

# print("a700")

 #vv <- "Mj / Aj" ; a1 <- Sys.time() ; cat( vv , a1-a0 , "\n") ; a0 <- Sys.time()			
		if ( fac.oldxsi>= 1){ fac.oldxsi <- 0 }
djj_old <- as.list( 1:J )
				
################################################################################
# some prelimaries for EM algorithm                                            #  
################################################################################

    # split item response pattern in a data frame with items as columns
    spl <- sapply( as.vector(item.patt[,1]), FUN = function(ii){ strsplit( ii, split = NULL) } ) 
    item.patt.split <- matrix( rep( 0, length(spl) * J ), ncol=J )
    for (ll in 1:length(spl) ){
        item.patt.split[ ll, ] <- as.numeric( spl[[ll]] )
        }

    # response pattern matrix: each observed entry corresponds to a 1, each unobserved entry to a 0
    resp.patt <- 1* ( item.patt.split != 9 )

    # number of item response patterns
    IP <- nrow(item.patt.split)           
   
    iter <- 1 # Iteration number
    likediff <- 1 # Difference in likelihood estimates
    loglike <- 0 # init for log-Likelihood
    
    # init value for maximum parameter change in likelihood maximization
    max.par.change <- 1000
    devchange <- 1000
	
	# response patterns
	cmresp <- colMeans( resp.patt )
    some.missings <- if( mean(cmresp) < 1){ TRUE } else { FALSE }
	
    # calculations for expected counts
	# response indicator list
    resp.ind.list <- list( 1:J )
	for (i in 1:J){ resp.ind.list[[i]] <- which( resp.patt[,i] == 1)  }

# print("a800")
	
	# this matrix is needed for computing R.lj
	if (G==1 ){
		ipr <- item.patt.split * item.patt.freq*resp.patt
				}
       
    disp <- "...........................................................\n"		

# print("B000")


# vv <- "item patt" ; a1 <- Sys.time() ; cat( vv , a1-a0 , "\n") ; a0 <- Sys.time()			
    
	#********************************
	# extract parameters with minimal deviances

    dev.min <- 10^99
	R.lj.gg <- I.lj.gg <- NULL
	suffstat_probs <- as.list(1:J)
	
	
	devchange <- 0	
	#*********************************

################################################################################
# BEGIN OF THE ITERATION LOOP                                                  #
################################################################################
    
    while ( ( iter <= maxit ) & 
				( ( max.par.change > conv.crit ) | ( devchange > dev.crit  ) )
					){
# a0 <- Sys.time()
################################################################################
# STEP I:                                                                      #
# calculate P(X_i | alpha_l):                                                  # 
# probability of each item response pattern given an attribute pattern         #
################################################################################

	if ( progress ){
		  cat(disp)	
		  cat("Iteration" , iter , "   " , paste( Sys.time() ) , "\n" )	   
				}

		pj1 <- matrix( 0 , nrow = J , ncol = L )
		# calculate P(X_j | alpha_l )
		for (jj in 1:J){
#			jj <- 3
			ajj <- ( aggr.attr.patt[[jj]] )
			mjjj <- Mj[[jj]][[1]][ ajj , ]
			djj <- matrix( delta[[jj]] , L , length(delta[[jj]]) , byrow=TRUE )
			pj1[jj,] <- rowSums( mjjj * djj )
			if (linkfct == "logit"){
				pj1[jj,] <- stats::plogis( pj1[jj,] )
									}
			if (linkfct == "log"){
				pj1[jj,] <- exp( pj1[jj,] )
#				if (max(pj1[jj,]>1){
#					pj1[jj,] <- pj1[jj,] / max( pj1[jj,] )
#								}
									}									
									}
# Revalpr("round(pj1,4)")
									
									
    # restrict probabilities in calculations									
	eps <- 10^(-10)
	pj1[ pj1 < 0 ] <- eps
	pj1[ pj1 > 1] <- 1 - eps	
	
# cat( "\n Step 1a (calculate P(X_j|alpha_l)\n" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
									
	pjM <- array( NA , dim=c(J,2,L) )
	pjM[,1,] <- 1 - pj1
	pjM[,2,] <- pj1
	h1 <- matrix( 1 , nrow=IP , ncol=L )
	
    p.xi.aj <- calc_posterior.v2(rprobs= pjM , gwt=h1 , resp=item.patt.split , 
								 nitems= J , 
                                 resp.ind.list=resp.ind.list , normalization=FALSE , 
                                 thetasamp.density= NULL , snodes=0 )[["hwt"]]	
#    p.xi.aj <- res.hwt[["hwt"]]  			

	if ( ! is.null(zeroprob.skillclasses) ){
		p.xi.aj[ , zeroprob.skillclasses ] <- 0
								}	
	
# cat( "\n Step 1 (calc individual likelihood)\n" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1	
	
################################################################################
# STEP II:                                                                     #
# calculate P(  \alpha_l | X_i ):                                              #
# posterior probability of each attribute pattern given the item response pattern
################################################################################

                                           
    # posterior probabilities  P( \alpha_l | X_i ) 
	if (G== 1){ 
		p.aj.xi <- outer( rep(1,IP), attr.prob ) * p.xi.aj 
			 } else {
			 p.aj.xi <- array( 0 , c( IP , L , G ) )
		for (gg in 1:G){
			p.aj.xi[,,gg] <- outer( rep(1,IP), as.vector(attr.prob[,gg]) ) * p.xi.aj
						}
				}
			 
	if (G == 1){ 
		if ( ! is.null( zeroprob.skillclasses ) ){
			p.aj.xi[ , zeroprob.skillclasses ] <- 0 
						}
		p.aj.xi <- p.aj.xi / rowSums( p.aj.xi )
		# calculate marginal probability P(\alpha_l) for attribute alpha_l
		if (! reduced.skillspace ){
			attr.prob <- colSums( p.aj.xi * item.patt.freq / I )
							}
				}
    if ( G > 1 ){ 					
			if ( ! is.null( zeroprob.skillclasses ) ){
			 for (gg in 1:G){ 
					p.aj.xi[ , zeroprob.skillclasses , gg ] <- 0 
							}
						}	
		for( gg in 1:G){
			p.aj.xi[,,gg] <- p.aj.xi[,,gg] / rowSums( p.aj.xi[,,gg] )
			Igg <- sum( item.patt.freq[,gg] )
			attr.prob[,gg] <- colSums( p.aj.xi[,,gg] * item.patt.freq[,gg] / Igg )
						}
					}

 # cat( "\n Step 2 (calc P(alpha|xi) \n" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1					


#######################################################################
# STEP II0: higher order GDINA model
#######################################################################

if (HOGDINA >= 0){
    for (gg in 1:G){ # gg <- 1
		if (G==1){ ap.gg <- attr.prob } else {
			ap.gg <- attr.prob[,gg] }
		res <- .attr.rpf( attr.patt , attr.prob=ap.gg , theta.k , wgt.theta[,gg] , HOGDINA )
		if (G==1){ attr.prob <- res$attr.prob } else {
			attr.prob[,gg] <- res$attr.prob }
		a.attr[,gg] <- res$a.attr
		b.attr[,gg] <- res$b.attr 
					}
				}

				
#######################################################################
# STEP IIa: reduction of skill space					
#######################################################################

# This currently only works in case of a single group
	if (reduced.skillspace){
		ntheta <- colSums( outer( item.patt.freq , rep( 1 , L) )*p.aj.xi )
		res <- gdina.reduced.skillspace( ntheta , Z , 
			      reduced.skillspace.method= reduced.skillspace.method )		
		beta <- res$beta
		attr.prob <- res$attr.prob
		pred.ntheta <- res$pred.ntheta
			}

#  cat( "\n Step 2a (reduced skillspace) \n" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1				
 
################################################################################
# STEP III:                                                                    #
# calculate I_{lj} and R_{lj}                                                  #
# for a derivation see De La Torre (2008, Journal of Educational and           #
# Behavioral Statistics)                                                       #
# I_{lj} ... expected frequency of persons in attribute class l for item j     #
#               (in case of no missing data I_{lj} = I_l for all items j       #
# R_{lj} ... expected frequency of persons in attribute class l for item j     #
#               which correctly solve item j                                   #
################################################################################

               
	if (G == 1){ 	
		R.lj <- I.lj <- matrix( 0 , nrow=J , ncol=L )
		if ( some.missings ){
			I.lj <- crossprod( item.patt.freq*resp.patt	, p.aj.xi )
					} else {
			I.lj <- matrix( t( item.patt.freq ) %*% p.aj.xi , nrow=J , 
							ncol=L , byrow=TRUE )
					}
												
		R.lj <- crossprod(ipr ,  p.aj.xi )
		colnames(I.lj) <- colnames(R.lj) <- attr.patt.c
		rownames(I.lj) <- rownames(R.lj) <- colnames(data)

				}	# end one group

				
	if (G > 1){ 					 
		R.lj.gg <- I.lj.gg <- array( 0 , c( J, L , G ) )
		for (gg in 1:G){ 	

		if ( some.missings ){
				I.lj.gg[,,gg] <- crossprod( item.patt.freq[,gg]*resp.patt , p.aj.xi[,,gg] )
					} else {
			I.lj.gg[,,gg] <- crossprod( item.patt.freq[,gg]*resp.patt , p.aj.xi[,,gg] )
					}

		    R.lj.gg[,,gg] <- crossprod( item.patt.split  * item.patt.freq[,gg] * resp.patt , p.aj.xi[,,gg] )
			colnames(I.lj.gg) <- colnames(R.lj.gg) <- attr.patt.c
			rownames(I.lj.gg) <- rownames(R.lj.gg) <- colnames(data)
						}
		# calculate I.lj and R.lj
		I.lj <- I.lj.gg[,,1]
		R.lj <- R.lj.gg[,,1]
		for (gg in 2:G){ 
			I.lj <- I.lj + I.lj.gg[,,gg]
			R.lj <- R.lj + R.lj.gg[,,gg]
						}
				}				

# cat( "\n Step 3 (calculate expected counts) \n" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1					
	
################################################################################
# STEP IV:                                                                     #
# M Step																	   # 
# GDINA Model																   #
################################################################################

	# calculation of expected counts
	R.ljM <- R.lj %*% aggr.patt.designmatrix
	I.ljM <- I.lj %*% aggr.patt.designmatrix


	eps <- 10^(-10)
	eps2 <- 10^(-10)

   max.increment <- max.increment / increment.factor
	
	delta.new <- NULL
	for (jj in 1:J){ 	# begin item
	#	jj <- 2 
		Ajjj <- Aj[[jj]]
		Mjjj <- Mj[[jj]][[1]]
		Rlj.ast <- R.ljM[ jj, Mj.index[jj,5]:Mj.index[jj,6] ]
		Ilj.ast <- I.ljM[ jj, Mj.index[jj,5]:Mj.index[jj,6] ]
		pjjj <- Rlj.ast / ( Ilj.ast + eps2 )
		suffstat_probs[[jj]] <- pjjj
		
		
		if (linkfct == "logit" ){ 
				pjjj[ pjjj > 1-eps ] <- 1 - eps
				pjjj[ pjjj < eps ] <- eps
				pjjj <- stats::qlogis( pjjj ) 
				# maxval <- 5 ;  pjjj <- squeeze.cdm( pjjj , c(-maxval , maxval ) )
								}
		#*****
	if (linkfct == "log" ){ 
				pjjj[ pjjj < eps ] <- eps
				pjjj <- log( pjjj ) 				
								}

								
		Wj <- diag( Ilj.ast )

		
		if ( avoid.zeroprobs ){
			 ind <- which( Ilj.ast  < 10^(-10)  )
			 if ( length(ind) > 0 ){
				 Wj <- diag( Ilj.ast[-ind] )
				 Mjjj <- Mjjj[ - ind , ]
				 pjjj <- pjjj[ - ind  ]
						}
					}
		
	    if ( ( rule[jj] == "GDINA" )| ( method == "ULS" ) ){ 
				invM <- invM.list[[jj]] 
#				delta.jj <- invM %*% t(Mjjj) %*% pjjj				
				delta.jj <- invM %*% crossprod(Mjjj ,pjjj)
							} else { 
#				invM <- solve( t(Mjjj) %*% Wj %*% Mjjj + diag( rep( eps2 , ncol(Mjjj) )) )
				invM <- solve( crossprod(Mjjj , Wj ) %*% Mjjj + diag( rep( eps2 , ncol(Mjjj) )) )				
#				delta.jj <- invM %*% t(Mjjj) %*% Wj %*% pjjj
				delta.jj <- tcrossprod( invM , Mjjj ) %*% Wj %*% pjjj
								}
		djj <- delta.jj[,1]
		djj.change <- djj - delta[[jj]]
		if (linkfct == "identity" & (iter > 3) ){ 
#		if ( (iter > 3) ){ 
			step.change <- .20
# 			djj.change <- ifelse( abs(djj.change) > step.change ,
#									step.change*sign(djj.change) , djj.change )
 			djj.change <- ifelse( abs(djj.change) > step.change ,
									djj.change / 2 , djj.change )
							}

									
			djj <- delta[[jj]] + djj.change
		if ( linkfct == "identity"){
				if ( sum(djj) > 1 ){ 	djj <- djj / ( sum( djj ) )  }											
									}

		#######################################################################
		iter_min <- 10
#		if (linkfct == "log" & iter > iter_min ){ 
#			if ( rule[jj] == "ACDM" ){
#				if ( sum( djj ) > 0 ){
#					djj <- djj - sum(djj )
#									}
#								}
#								}				
								
								
		djj <- ifelse ( is.na(djj) , delta[[jj]] , djj )
		
#		if ( fac.oldxsi > 0 & (iter > 1 ) ){ 
#				djj <- ( 1 - fac.oldxsi ) * djj + fac.oldxsi * djj_old[[jj]]
#				djj_old[[jj]] <- djj				
#								}
		
		# control
        djj.change <- djj - delta[[jj]]		
		while( max(abs(djj.change)) > max.increment ){
#					djj.change <- djj.change / 2 
				djj.change <- ifelse( abs(djj.change) > max.increment , djj.change / 2 , djj.change )
						}
		djj <- delta[[jj]] + djj.change						

		if ( rrum.model & (iter > 10) ){

	#---
	#  RRUM parametrization
	#  log( P(X=1) ) = b0 + b1*alpha1 + b2 * alpha2 
	#  RRUM:
	#  P(X=1) = pi * r1^( 1- alpha1) * r2^(1-alpha2)
	#  => log( P(X=1) ) = log[ pi * r1 * r2 * r1^(-alpha1) * r2^(-alpha2) ]
	#                   = log( pi ) + log(r1) + log(r2) + -log(r1)*alpha1 + -log(r2) * alpha2
	#  => b1 = -log(r1) and r1 = exp( -b1 )
	#  => log(pi) = b0 + b1 + b2 and pi = exp( b0 + b1 + b2 )		
			
			d1 <- djj
# d01 <- d1	
#			d1 <- ifelse( d1 < 0 , .1 , d1 )						
	        sum_d1 <- sum(d1)
			if ( sum_d1 > 0 ){
                d1 <- d1 - sum(d1)
							}

			d1_samp <- stats::runif( length(d1) , 0 , .01 )
			d1[-1] <- ifelse( d1[-1] < 0 , d1_samp[-1] , d1[-1] )						
							
	        sum_d1 <- sum(d1)
			if ( sum_d1 > 0 ){
                d1 <- d1 - sum(d1)
							}														
			djj <- d1										
#			d1 <- djj[-1] 
#			d1 <- ifelse( d1 < 0 , 0.01 , d1 )
#			djj[-1] <- d1			
#			if ( djj[1] > 0 ){
#				djj[1] <- 0
#								}						


						}


		delta.new[[jj]] <- djj
		if ( (fac.oldxsi > 0 ) & (iter>3)){
		    fac.oldxsi1 <- fac.oldxsi * ( devchange >= 0 )
			delta.new[[jj]] <- fac.oldxsi1*delta[[jj]] + ( 1 - fac.oldxsi1 ) * delta.new[[jj]]
						}

		# fix delta parameter here!!
		if ( ! is.null( delta.fixed ) ){
			delta.fixed.jj <- delta.fixed[[jj]]
			if ( ! is.na( delta.fixed.jj)[1] ){
					delta.new[[jj]] <- delta.fixed.jj
									}
							}

					}		# end item
	#.............................................................					

	
	
	##########################################################################
	# estimation with a design matrix for delta parameters
	##########################################################################
	if ( ! is.null( delta.designmatrix ) ){ 
		u.delta.new <- unlist( delta.new )
		# calculate basis parameter of delta
		delta.basispar <- solve( t( delta.designmatrix) %*% delta.designmatrix ) %*% 
								t(delta.designmatrix) %*% u.delta.new
		if ( ! is.null( delta.basispar.lower )){						
			delta.basispar <- ifelse( delta.basispar < delta.basispar.lower , 
										delta.basispar.lower , delta.basispar )
											}
		if ( ! is.null( delta.basispar.upper )){						
			delta.basispar <- ifelse( delta.basispar > delta.basispar.upper , 
										delta.basispar.upper , delta.basispar )
									}
		delta.new1 <- ( delta.designmatrix %*% delta.basispar )[,1]
		for (jj in 1:J){
			delta.new[[jj]] <- delta.new1[ seq( Mj.index[jj,2] , Mj.index[jj,3] ) ]
						}
									}

									
# cat( "\n Step 4 (m step item parameters) \n" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1
# stop("here")									
	#################################################
	
    # calculate the updated likelihood    
	p.xi.aj[ p.xi.aj > 1 ] <- 1-10^(-30)
	p.xi.aj[ p.xi.aj < 0 ] <- 10^(-30)	
	if (G==1){ 	
		l1 <- rowSums( p.xi.aj * outer( rep(1,IP), attr.prob )  ) + 10^(-30) 
		l1[ l1 < 0 ] <- 10^(-30)
			}
	if (G>1){
		l1 <- matrix( 0 , IP , G )
		for (gg in 1:G){ 
			l1[,gg] <- rowSums( p.xi.aj * outer( rep(1,IP), attr.prob[,gg] )  ) + 10^(-30) 
			l1[ l1[,gg] < 0 ,gg] <- 10^(-30)
					}
				}
	like.new <- sum( log( l1 ) * item.patt.freq ) 
    likediff <- abs( loglike - like.new )
	loglikeold <- loglike
    loglike <- like.new
    # maximum parameter change
    max.par.change <- max( abs( unlist( delta.new ) - unlist( delta ) ) )
	if ( linkfct %in% c("logit","log") ){
				max.par.change <- max( abs( plogis(unlist( delta.new )) -
										plogis( unlist( delta ) )) )
								}
	
    # define estimates which are updated in this iteration
    delta <- delta.new
	
    if (progress) {  
	if (progress.item){ 
			g1 <- unlist( lapply( delta , FUN = function(ll){ paste( round(ll,4) , collapse= " " ) } ))
			g1 <- matrix( paste( colnames(data) , g1 ) , ncol=1)
			print(g1)
			}
		cat( "Deviance = "  , round( -2*like.new , 5 ) )
        devchange <- g11 <- 2*(like.new-loglikeold)		
			if (iter >1){ cat(" | Deviance change = " , round( 2*(like.new-loglikeold), 7) ) }
			cat("\n" )
			if (g11 < 0 ){ cat( "**** Deviances decreases! Check for nonconvergence.   ****\n") }
		cat("Maximum parameter change:" , round( max.par.change, 6), "\n") 			
			}
    
    flush.console() # Output is flushing on the console
    iter <- iter + 1 # new iteration number                                    
	devchange <- abs( 2*(like.new-loglikeold) )

# cat( "\n Step 5 likelihood \n" ) ; a1 <- Sys.time() ; print(a1-a0) ; a0 <- a1									
# stop()    
	
	#******************************
	# update parameters at minimal deviance
	dev <- -2*like.new
	
	if (save.devmin){
	
	if ( dev < dev.min ){
		iter.min <- iter-1	
		delta.min <- delta
		dev.min <- dev
		p.aj.xi.min <- p.aj.xi
		p.xi.aj.min <- p.xi.aj
		R.lj.min <- R.lj
		I.lj.min <- I.lj		
		attr.prob.min <- attr.prob
		loglike.min <- loglike
				}		
				}
	#********************************
		
	}

	
	
################################################################################
# END OF THE ITERATION LOOP                                                    #
################################################################################

	#***************************************
	# use parameters with minimal deviance
	iterused <- iter - 1
	if (save.devmin){
		iter.min -> iter	
		delta.min -> delta
		dev.min -> dev
		p.aj.xi.min -> p.aj.xi
		p.xi.aj.min -> p.xi.aj
		R.lj.min -> R.lj
		I.lj.min -> I.lj		
		attr.prob.min -> attr.prob
		loglike.min -> loglike		
		}
	#****************************************


    # calculate posterior probability for each attribute pattern
	if (G==1){
	
		# set likelihood for skill classes with zero probability to zero
	 if ( ! is.null(zeroprob.skillclasses) ){
		p.xi.aj[ , zeroprob.skillclasses ] <- 0
								}
#		pattern <- cbind( 
		pattern <- data.frame( 
						freq = round(as.numeric(item.patt[,-1]),3),
						mle.est = attr.patt.c[ max.col( p.xi.aj ) ], 
						mle.post = rowMaxs( p.xi.aj ) / rowSums( p.xi.aj ), 
						map.est = attr.patt.c[ max.col( p.aj.xi ) ], 
						map.post = rowMaxs( p.aj.xi )
						)
				}

				
	if (G>1){
		ind1 <- match( item.patt.subj , item.patt[,1] )
		l1 <- attr.patt.c[ max.col( p.xi.aj ) ]
		pattern <- data.frame( "mle.est" = l1[ind1] ) 
		l1 <- rowMaxs( p.xi.aj ) / rowSums( p.xi.aj )
		pattern$mle.post <- l1[ind1]
		pattern$map.est <- NA
		pattern$map.post <- NA
		for (gg in 1:G){
			# gg <- 1	
			ind.gg <- which( group2 == gg )
			ind2.gg <- match( item.patt.subj[ind.gg] , item.patt[  , 1] )
			l1 <- attr.patt.c[ max.col( p.aj.xi[,,gg] ) ]
			pattern$map.est[ind.gg] <- l1[ind2.gg]
			l1 <-  rowMaxs( p.aj.xi[,,gg] )
			pattern$map.post[ind.gg] <- l1[ind2.gg]		
					}
			}



    # calculate posterior probabilities for all skills separately
	if (G==1){
		attr.postprob <- p.aj.xi %*% attr.patt
		colnames( attr.postprob ) <- paste("post.attr",1:K, sep="")
		pattern <- cbind( pattern,  attr.postprob )
			}
	
	#####################################################
	# itemwise standard error calculation
	
	varmat.delta <- varmat.palj <-  NULL
	se.delta <- NULL	

	
	delta.summary <- NULL

	if (G == 1){ 	
		PAJXI <-  p.aj.xi
				}
	if (G>1){			
		a1 <- outer( rep(1,nrow(attr.prob) ) , colSums( item.patt.freq ) ) / sum( item.patt.freq)
		attr.prob.tot <- rowSums( attr.prob * a1 )
		PAJXI <- outer( rep(1,IP), attr.prob.tot ) * p.xi.aj
		#	p.aj.xi <- p.aj.xi / rowSums( p.aj.xi )
		PAJXI <- PAJXI / rowSums(PAJXI)
			}

	# matrix form of item.patt.freq
	if (G==1){ item.patt.freq <- matrix( item.patt.freq , ncol=1 ) }
	freq.pattern <- rowSums( item.patt.freq )
	
    # if ( calc.se ){
	
	for (jj in 1:J){	
		se.jj <- NA
		if ( calc.se ){
#	 cat("........",jj,".,,,\n")
			#	jj <- 1		# Item jj
				Ajjj <- Aj[[jj]]
				Mjjj <- Mj[[jj]][[1]]
				Rlj.ast <- stats::aggregate( R.lj[jj,] , list( aggr.attr.patt[[jj]]) , sum )
				Ilj.ast <- stats::aggregate( I.lj[jj,] , list( aggr.attr.patt[[jj]]) , sum )
				pjjj <- Rlj.ast[,2] / Ilj.ast[,2]
				Mjj2 <- Mj[[jj]][[2]]
				apjj <- aggr.attr.patt[[jj]] 			
				Mjjj <- Mjjj[ sort(unique(apjj)) , ]

				# M1 <- max( apjj )
				M1 <- length( unique(apjj) )
				p.ajast.xi <- matrix( 0 , nrow=IP , ncol = M1 )
				for (kk in 1:M1){
					pg1 <-  PAJXI[ , apjj == kk  ]					
					if ( is.vector(pg1)){ 
								p.ajast.xi[,kk] <- pg1 
									} else {
								p.ajast.xi[,kk] <- rowSums( pg1 ) 
										}
								}	
				Rlj.ast <- stats::aggregate( R.lj[jj,] , list( aggr.attr.patt[[jj]]) , sum )
				Ilj.ast <- stats::aggregate( I.lj[jj,] , list( aggr.attr.patt[[jj]]) , sum )
				pjjj <- Rlj.ast[,2] / Ilj.ast[,2]		
				pjjjM <- outer( rep(1,IP) , pjjj ) + 10^(-20)		
				nM <- ncol(pjjjM) 
				x1 <- outer( item.patt.split[,jj] , rep(1,nM) )
				r1 <- outer( resp.patt[,jj] * item.patt.freq , rep(1,ncol(pjjjM) ) )
				# Formula (17) for calculating the standard error	
				mat.jj <- p.ajast.xi * ( x1 - pjjjM) / ( pjjjM * ( 1 - pjjjM ) )	
				infomat.jj <- matrix( 0 , nM , nM )
				for (kk1 in 1:nM){
					for (kk2 in kk1:nM){ 
						# kk1 <- 1
						# kk2 <- 1
#							infomat.jj[kk2,kk1] <- infomat.jj[kk1,kk2] <-  
#											sum( mat.jj[,kk1] * mat.jj[,kk2]  )
						#@@ARb (2012-07-20) correction
						# frequency weights must be taken into account
						hh1 <- sum( mat.jj[,kk1] * mat.jj[,kk2] * freq.pattern * 
											resp.patt[,jj] * item.patt.split[,jj] )
						infomat.jj[kk2,kk1] <- infomat.jj[kk1,kk2] <-  hh1
										}
									}
				a1 <- NULL
		if ( avoid.zeroprobs ){
			 ind <- which( is.na(diag(infomat.jj) ))
			 if ( length(ind) > 0 ){
				infomat.jj <- infomat.jj[-ind, -ind]	
						}			 
			 
					}				
		
#				try( a1 <- solve( infomat.jj ) )
				a1 <- try( solve( infomat.jj + diag( eps2 , ncol(infomat.jj) ) ) )
				if ( is(a1 , "try-error") ){ 
						cat( "Item" , colnames(data)[jj] , "Singular item parameter covariance matrix\n")
						a1 <- NA*infomat.jj 
							}
				varmat.palj[[jj]] <- Ijj <- a1
				Wj <- diag( Ilj.ast[,2] )	
			
		if ( avoid.zeroprobs ){
			 ind <- which( Ilj.ast[,2]  < 10^(-10)  )
			 if ( length(ind) > 0 ){
				 Wj <- diag( Ilj.ast[-ind,2] )
				 Mjjj <- Mjjj[ - ind , ]
				 pjjj <- pjjj[ - ind  ]
						}
					}

				if ( ( method == "ULS" ) ){ 			
				    x1 <- t(Mjjj) %*% Mjjj	
				    diag(x1) <- diag(x1) + 10^(-8)						
					Wjjj <- solve( x1 ) %*% t(Mjjj)
					
									} else {
					x1 <- t(Mjjj) %*% Wj %*% Mjjj									
					diag(x1) <- diag(x1) + 10^(-8)					
					Wjjj <- solve( x1 ) %*% t(Mjjj) %*% Wj
											}
				if ( linkfct == "logit" ){
					pjjj.link <- 1 / ( ( pjjj * ( 1 - pjjj ) ) + eps2 )
					pjjj.link <- diag( pjjj.link )
				    Wjjj <- Wjjj %*% pjjj.link
						}
				if ( linkfct == "log" ){
					pjjj.link <- 1 /  ( pjjj  + eps2 )
					pjjj.link <- diag( pjjj.link )
					Wjjj <- Wjjj %*% pjjj.link
						}
				varmat.delta[[jj]] <- Wjjj %*% Ijj %*% t(Wjjj)		
				se.jj <- sqrt( diag(varmat.delta[[jj]] )  ) 
								}
						
				delta.summary.jj <-
					data.frame( "link" = linkfct , "item" = colnames(data)[jj] , 
								"itemno" = jj , 
								"type" = Mj[[jj]][2] , 
								"rule" = rule[jj] , 
								"est" = delta[[jj]] , 
								"se" = se.jj
								)
								
		# fix delta parameter here!!
		if ( ! is.null( delta.fixed ) ){
			delta.fixed.jj <- delta.fixed[[jj]]
			if ( ! is.na( delta.fixed.jj)[1] ){
					delta.summary.jj$se <- 0
									}
							}								
								
				colnames(delta.summary.jj)[4] <- "partype"					
				delta.summary <- rbind( delta.summary , delta.summary.jj )
				
				
				}
				

	
	delta.summary$partype.attr <- paste(delta.summary$partype)
	if (calc.se){		
	for (jj in 1:J){
		ind.jj <- which( delta.summary$itemno == jj )
		qjj <- which( q.matrix[ jj , ]	> 0 )
		pgjj <- pajj <- paste(delta.summary$partype.attr[ind.jj])
		cjj <- paste(colnames(q.matrix)[qjj])
		NN <- length(pajj)
		pajj <- gsub( "|" , "-" , pajj )
		pajj <- gsub( "=" , "-" , pajj )
		for (nn in 1:NN){
			st1 <- as.numeric(unlist( strsplit( paste(pajj[nn]) , "-" ) ))
			st1 <- st1[ ! is.na( st1 ) ]
			st1 <- st1[ st1 > 0 ]
			pgjj[nn] <- paste( cjj[ st1 ] , collapse="-" )
						}
		delta.summary$partype.attr[ind.jj] <- pgjj
					}
			}
							
	
	# compute RRUM parametrization if model is specified
	if (rrum.model){
		rrum.params <- .rrum.param( delta.summary , q.matrix )
				}
				
    # attribute pattern
	if (G==1){ 
		attr.prob <- matrix( attr.prob, ncol=1)
		colnames( attr.prob ) <- "class.prob"		
			}
	if (G>1){
		colnames( attr.prob ) <- paste( "class.prob.group" , 1:G , sep="")
				}
		rownames( attr.prob ) <- attr.patt.c

	
	mA <- max( maxAttr)
	
	if (G==1){   
	    sp <- NULL 
		# pattern for separate skills
#		for (kk in 0:( mA -1) ){
		for (kk in 0:mA ){
	#		kk <- 0
			skill.patt <- matrix(apply( matrix( rep(  attr.prob, K ), ncol=K) * 
					(attr.patt==kk), 2, sum ),ncol=1)
	#		rownames(skill.patt) <- paste("Skill_", colnames(q.matrix),sep="")
			rownames(skill.patt) <- colnames(q.matrix)
			colnames(skill.patt) <- paste0("skill.prob" ,kk )
			sp <- cbind( sp , skill.patt )
						}
            skill.patt <- sp
			for (kk in 1:K){ 
				ind.kk <- setdiff( 1:mA , 1 + q.entries[[kk]] )
				if ( length(ind.kk) > 0 ){ 	skill.patt[ kk ,ind.kk ] <- NA 	}	
							}
				}
	if (G>1){   
	sp <- NULL
	for (kk in 0:( mA ) ){	
	# kk <- 0
		skill.patt <- matrix( 0 , K , G )
		for (gg in 1:G){
		skill.patt[,gg] <- matrix(apply( matrix( rep(  attr.prob[,gg], K ), ncol=K) * 
									( attr.patt == kk ), 2, sum ),ncol=1)
						}
		#		rownames(skill.patt) <- paste("Skill_", colnames(q.matrix),sep="")
		rownames(skill.patt) <- colnames(q.matrix)
		colnames(skill.patt) <- paste0( "skill.prob" , kk , ".group"  , 1:G )				
		sp <- cbind( sp , skill.patt )
						}
        skill.patt <- sp

			for (kk in 1:K){ 
			    v1 <- rep(1:mA,each=G)
				ind.kk <- setdiff( v1 , rep(1 + q.entries[[kk]],each=G) )
				ind.kk <- which( v1 %in% ind.kk )
				if ( length(ind.kk) > 0 ){ 	skill.patt[ kk ,ind.kk ] <- NA 	}	
							}
		
				}
				
				
		#####################################################################		
		# calculation of the AIC und BIC        
		bb <- 0
	#    if ( is.null( constraint.guess ) == F ){  bb <- bb + nrow(constraint.guess) }
	#    if ( is.null( constraint.slip ) == F ){  bb <- bb + nrow(constraint.slip) }
	
		Nipar <- length( unlist( delta) )
		if ( ! is.null( delta.designmatrix ) ){ 
			Nipar <- ncol(delta.designmatrix ) 
			}
		if ( ! is.null( delta.fixed) ){
				Nipar <- Nipar - sum(1 - is.na( unlist( delta.fixed )) )										
							}
		
		
		Nskillpar <- G*ncolZ - length( zeroprob.skillclasses )	
		
		if (HOGDINA==1){ Nskillpar <- 2*K*G }
		if (HOGDINA==0){ Nskillpar <- K*G }
		Npars <- Nipar  - bb + Nskillpar
		II <- sum( item.patt.freq )
		aic <- -2*loglike + 2 * Npars  
		bic <- -2*loglike + Npars*log(II)
		#         ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
		caic <- -2*loglike + ( log(II) + 1 ) * Npars

	if (G==1){	
		rownames( p.aj.xi ) <- rownames( pattern ) # output rownames posterior probabilities    
		pattern <- data.frame(pattern) # convert pattern to numeric format
		for (vv in seq(1,ncol(pattern))[ -c(2,4) ] ){
						pattern[,vv ] <- as.numeric( paste( pattern[,vv] ) ) }
		
		# subject pattern
		item.patt.subj <- data.frame( "case" = 1:(nrow(data) ), 
									   "pattern" = item.patt.subj, 
                                       "pattern.index" = match( item.patt.subj, item.patt[,1] )
												)
											
		# attribute pattern (expected frequencies)
		attr.prob0 <- attr.prob
		attr.prob <- data.frame( attr.prob )
		attr.prob$class.expfreq <-  attr.prob[,1] * nrow(data) 
		
		#*****
		# modify output (ARb 2012-06-05)
		pattern <- pattern[ item.patt.subj$pattern.index , ]	
		pattern[,1] <- paste( item.patt.subj$pattern )
		colnames(pattern)[1] <- "pattern"
		p.aj.xi <- p.aj.xi[ item.patt.subj$pattern.index , ]
		rownames(p.aj.xi) <- pattern$pattern
		p.xi.aj <- p.xi.aj[ item.patt.subj$pattern.index , ]
		rownames(p.xi.aj) <- pattern$pattern
	
		#*****				
				}
		if (G==1){
			posterior <- p.aj.xi
					}
		if (G>1){
			ind <- match( item.patt.subj , item.patt[,1] )			
			# colnames(pattern)[1] <- "pattern"		
			p.xi.aj <- p.xi.aj[ ind , ]
			rownames(p.xi.aj) <- pattern$pattern
			p.aj.xi <- p.aj.xi[ ind , , ]
			rownames(p.aj.xi) <- pattern$pattern		
		
			ND <- dim(p.aj.xi)
			posterior <- matrix( 0 , nrow=ND[1] , ncol=ND[2] )
		    for (gg in 1:G){
				ind.gg <- which( group == gg )
				posterior[ ind.gg , ] <- p.aj.xi[ ind.gg , , gg ]
							}
			attr.prob0 <- attr.prob							
				
				}


				
	############################################
	# item fit [ items , theta , categories ] 
	# # n.ik [ 1:TP , 1:I , 1:(K+1) , 1:G ]
	n.ik <- array( 0 , dim=c(L , J , 2 , 1 ) )
	n.ik[  , , 2 , 1 ] <- t(R.lj)
	n.ik[  , , 1 , 1 ] <- t(I.lj-R.lj)
	pi.k <- array( 0 , dim=c(L,1) )
	
	if (G>1){	# item fit only in multiple group case
		g1 <- colSums( item.patt.freq )
		g1 <- g1 / sum(g1) 
		for (gg in 1:G){		
			pi.k[,1] <- pi.k[,1] + attr.prob[,gg] * g1[gg]
				}
			} 
		
	if (G==1){	# item fit only in one group case
		pi.k[,1] <- attr.prob$class.prob
			} 
		probs <- aperm( pjM , c(3,1,2) )
		itemfit.rmsea <- itemfit.rmsea( n.ik , pi.k , probs )$rmsea	
		names(itemfit.rmsea) <- colnames(data)
	#********************************************************	
	# calculate model implied probabilities	

#	if ( calc.se ){ 
		probitem <- gdina.probitem( Mj, Aj , delta , rule , linkfct , delta.summary )
#					} else {
#		probitem <- NULL
#					}
	
	
	# labels likelihood
	colnames(p.xi.aj) <- paste(rownames(attr.prob))
	

	#************** OUTPUT **********************************
	if (progress){
		cat("---------------------------------------------------------------------------------\n")
			}
	iter <- iterused
    res <- list( coef = delta.summary , "item" = delta.summary , 
				delta = delta , se.delta = se.delta , 
			    "probitem"=probitem , 
				"itemfit.rmsea" = itemfit.rmsea , 
				"mean.rmsea" = mean(itemfit.rmsea) ,	
				loglike = loglike, deviance = -2*loglike , G = G , N = colSums( as.matrix(item.patt.freq) ) , 
				AIC = aic, BIC = bic, CAIC = caic , Npars  = Npars , 
				Nipar=Nipar  , Nskillpar = Nskillpar ,
				Nskillclasses = L , 	
				varmat.delta = varmat.delta ,  varmat.palj = varmat.palj ,
                 posterior = posterior , "like" = p.xi.aj, "data" = data, "q.matrix" = q.matrix,
                 pattern = pattern , attribute.patt = attr.prob, skill.patt = skill.patt,
                 "subj.pattern" = item.patt.subj, "attribute.patt.splitted" = attr.patt, 
				 "pjk" = pjM , 
				 Mj = Mj , Aj = Aj , "rule"=rule , "linkfct"=linkfct ,
				 delta.designmatrix = delta.designmatrix , 
				 "reduced.skillspace" = reduced.skillspace , 
				 "Z.skillspace" = if(reduced.skillspace){ Z } else { NULL } , 
#				 "delta.index" = if( reduced.skillspace){ sum( abs( ntheta-pred.ntheta) )/ (2*sum(ntheta)) } ,
				 beta = beta , covbeta = covbeta , 
				 "display" = disp,
                 "item.patt.split" = item.patt.split, "item.patt.freq" = item.patt.freq,
                 "model.type" = r1 , "iter" = iter , 
				 "iterused" = iterused ,
				 "rrum.model" = rrum.model ,
				 "rrum.params"= rrum.params ,
				 "group.stat" = group.stat , 
				 "NAttr" = maxAttr , 
#				 "q.matrix" = q.matrix ,
				 "HOGDINA" = HOGDINA ,
				 "seed"= seed  ,
				 iter = iter , 
				 "converged" = iter < maxit 
				 )
		 
	if (HOGDINA>=0) { 
	    colnames(a.attr) <- paste0( "a.Gr" , 1:G )
		colnames(b.attr) <- paste0( "b.Gr" , 1:G )
		rownames(b.attr) <- rownames(a.attr) <- colnames(q.matrix)
		res$a.attr <- a.attr 
		res$b.attr <- b.attr
		res$attr.rf <- cbind( b.attr , a.attr )
				}						
	# computation time
    time1$s2 <- Sys.time()
	res$time <- time1
	# res$time$timediff <- print(res$time$s2 - res$time$s1)	
	res$time$timediff <- res$time$s2 - res$time$s1	
	if ( progress ){
		print(res$time$s2 - res$time$s1)	
						}
	
	# control parameter
	control <- list( skillclasses=skillclasses , q.matrix=q.matrix, conv.crit = conv.crit , 
					dev.crit = dev.crit , maxit = maxit ,
					linkfct = linkfct , Mj = Mj , Aj = Aj , 
					group = group , 
					method = method , 
					delta.designmatrix = delta.designmatrix , 
					delta.basispar.lower = delta.basispar.lower , 
					delta.basispar.upper = delta.basispar.upper , 					
					delta.basispar.init = delta.basispar.init , 
					zeroprob.skillclasses = zeroprob.skillclasses , 
					reduced.skillspace=reduced.skillspace , 
					HOGDINA = HOGDINA , 
					Z.skillspace = Z.skillspace , 
                    weights = weights ,  rule = rule ,
					I.lj=I.lj , R.lj=R.lj , I.lj.gg = I.lj.gg , 
					R.lj.gg = R.lj.gg , aggr.patt.designmatrix=aggr.patt.designmatrix ,
					Mj.index=Mj.index , method=method ,
					aggr.attr.patt=aggr.attr.patt,IP=IP,
					p.aj.xi=p.aj.xi,item.patt.split=item.patt.split,
					resp.patt=resp.patt,freq.pattern=freq.pattern ,
					item.patt.freq=item.patt.freq,invM.list=invM.list ,
					suffstat_probs = suffstat_probs , 					
					increment.factor = increment.factor ,
					fac.oldxsi = fac.oldxsi ,
					avoid.zeroprobs = avoid.zeroprobs ,
					attr.prob = attr.prob0	,
					delta.fixed = delta.fixed
						) 	
	res$control <- control	
	
	
    # create parameter table						
	res$partable <- gdina.partable(res)
	
	# polychoric correlations
	res$polychor <- CDM.calc.polychor( res )	
	res$call <- cl
    class(res) <- "gdina"
    return(res)
}
##################################################################



