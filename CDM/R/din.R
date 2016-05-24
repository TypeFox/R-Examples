
#************************************************************************
# Function rowProds2 now included in din
# more efficient than rowProds
rowProds2 <- function(matr){
	y <- matr[,1]
	for (ii in 2:dim(matr)[2]){
		y <- y * matr[,ii] }
	return(y)
		}
#...................................................................


################################################################################
# Main function for parameter estimation in cognitive diagnosis models         #
################################################################################

din <-
function( data, q.matrix, skillclasses = NULL , conv.crit = 0.001, dev.crit = 10^(-5) , maxit = 500,
                    constraint.guess = NULL, constraint.slip = NULL,
                    guess.init = rep(.2, ncol(data) ), slip.init = guess.init,
                    guess.equal = FALSE , slip.equal = FALSE , 
					zeroprob.skillclasses = NULL , 
                    weights = rep( 1, nrow( data ) ),  rule = "DINA", 
					wgt.overrelax = 0 , 
					wgtest.overrelax = FALSE , 
					param.history = FALSE , 
					seed = 0 , 
                    progress = TRUE , guess.min=0 , slip.min=0, guess.max=1 , slip.max=1){
                    
# data: a required matrix of binary response data, whereas the items are in the columns 
#       and the response pattern in the rows. NA values are allowed.
#
# q.matrix: a required binary matrix describing which attributes are required, coded by 1,
#       and which attributes are not required, coded by 0, to master the items, whereas the
#       attributes are in the columns and the items in the rows.
#
# conv.crit: termination criterion of the iterations defined as the maximum change in parameter
#       estimates. Iteration ends if maximal parameter change is below this value.
#
# maxit: maximal number of iterations.
#
# constraint.guess: an optional matrix of fixed guessing parameters. The first column of this 
#       matrix indicates the items whose guessing parameters are fixed and the second column 
#       the values the guessing parameters are fixed to.
#
# constraint.slip: an optional matrix of fixed slipping parameters. The first column of this matrix 
#       indicates the items whose guessing parameters are fixed and the second column the values 
#       the guessing parameters are fixed to.
#
# guess.init: an optional initial vector of guessing parameters. Guessing parameters are bounded between
#       0 and 1.
#
# slip.init: an optional initial vector of guessing parameters. Slipping parameters are bounded between
#       0 and 1.
#
# guess.equal: an optional logical indicating if all guessing parameters are equal to each other
#
# slip.equal: an optional logical indicating if all slipping parameters are equal to each other
#
# zeroprob.skillclasses:  an optional vector of integers which indicates which skill classes should have
#							zero probability
#
# weights: an optional vector of weights for the response pattern. Noninteger weights allow for different
#       sampling schemes.
#
# wgt.overrelax ... convergence weight (overrelaxed EM) => w = 0 standard EM

# rule: an optional character string or vector of character strings specifying the model rule that is used. 
#       The character strings must be of "DINA" or "DINO". If a vector of character strings is specified, 
#       implying an itemwise condensation rule, the vector must be of length ncol(data). The default is the 
#       condensation rule "DINA" for all items.
#
# progress: an optional logical indicating whether the function should print the progress of iteration.

	if ( progress){
		cat("---------------------------------------------------------------------------------\n")
				}
	cl <- match.call()
	
################################################################################
# check consistency of input (data, q.matrix, ...)                             #
################################################################################

	I <- ncol(data) 
	if ( is.null( colnames(data) ) ){
		colnames(data) <- paste0("I" , 1:I )
					}

	# different inits if seed larger than zero
    if ( seed > 0 ){
		set.seed(seed) 
		slip.init <- stats::runif( I , 0 , .4 )
		guess.init <- stats::runif( I , 0 , .4 )
					}
					
    clean <- check.input(data, q.matrix, conv.crit, maxit, constraint.guess,
        constraint.slip, guess.init, slip.init, weights, rule, progress)   

		
    if (is.character(clean)) return(clean)
    
    dat.items <- clean$data; q.matrix <- clean$q.matrix; conv.crit <- clean$conv.crit;
    maxit <- clean$maxit; constraint.guess <- clean$constraint.guess; 
    constraint.slip <- clean$constraint.slip; guess.init <- clean$guess.init;
    slip.init <- clean$slip.init; weights <- clean$weights; rule <- clean$rule;
    progress <- clean$progress    

	
################################################################################
# model specification: DINA, DINO or itemwise specification of DINA or DINO    #
################################################################################
	if ( length(rule) == 1){
		if ( rule == "DINA" ){ r1 <- "DINA MODEL"  }
		if ( rule == "DINO" ){ r1 <- "DINO MODEL" }
			} else { r1 <- "Mixed DINA & DINO Model" }
				
################################################################################
# display on R console                                                         #
################################################################################

    disp <- r1      
    s1 <- Sys.time()
	if (progress){
	    cat(disp,"\n")
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
    L <- 2^K               # number of latent class pattern of attributes
    dat.items <- as.matrix( dat.items)
    q.matrix <- as.matrix( q.matrix)
            
################################################################################
# Initialization and missing data handling                                     #
################################################################################

    # initialize guessing and slipping parameters
    # without constraints, the default is set equal to .2 for all items
    guess <- guess.init ; slip <- slip.init
    
    # missing data is coded by 9
    resp <- 1 - is.na(dat.items)
    dat.items[ resp == 0 ] <- 9
    
    # standardize weights such that the sum of defined weights is equal to the number of rows in the data frame
    weights <- nrow(dat.items)*weights / sum(weights )

################################################################################
# calculate item response patterns                                             #
################################################################################

    # string with item response patterns
#    item.patt.subj <- sapply( 1:I, FUN = function(ii){ paste( dat.items[ ii, ], collapse="" )  } )	
	item.patt.subj <- dat.items[,1]
	for ( jj in 2:J){
		item.patt.subj <- paste( item.patt.subj  , dat.items[,jj] , sep="")
				}
	
	
    # calculate frequency of each item response pattern
	a2 <- rowsum( rep(1,I) , item.patt.subj) 
	item.patt <- a2[,1]	
	
    # sort item response pattern according to their absolute frequencies
#    six <- sort( item.patt, index.return=F, decreasing=T)
	six <- item.patt	
    # define data frame 'item.patt' with item response pattern and its frequency (weight)
#    item.patt <- cbind( "pattern" = rownames(six), "freq" = as.numeric(as.vector( six ) ) )
    item.patt <- cbind( "pattern" = names(six), "freq" = as.numeric(as.vector( six ) ) )


    # calculate weighted frequency for each item response pattern
#     item.patt[,2] <- sapply( 1:( nrow(item.patt) ), FUN = function(kk){
#                         sum( weights * ( item.patt[ kk, 1] == item.patt.subj  ) )
#                         } )   

	h1 <- rowsum( weights , item.patt.subj )	
	item.patt[,2] <- h1[ match( item.patt[,1] , rownames(h1) ) , 1]
    item.patt.freq <- as.numeric(item.patt[,2])
	
################################################################################ 
# generate all attribute patterns                                              #
################################################################################

    attr.patt <- matrix( rep( 0, K*L) , ncol=K)
    h1 <- 2
	if (K>=2 ){
		for(ll in 1:(K-1) ){
			lk <- combn( 1:K, ll ) 
#			lk
			for ( jj in 1:( ncol(lk) ) ){ 
				attr.patt[ h1, lk[,jj] ] <- 1
				h1 <- h1 + 1
				}
			}
		}
    attr.patt[ L, ] <- rep( 1, K )
	if ( ! is.null( skillclasses) ){ 
		attr.patt <- skillclasses
		L <- nrow(attr.patt)
			}		
		
    # combine all attributes in an attribute pattern as a string
    attr.patt.c <- apply( attr.patt, 1, FUN = function(ll){ paste(ll,collapse="" ) } )

################################################################################
# uniform prior distribution of all latent class patterns                      #
################################################################################

    attr.prob <- rep( 1/L, L )
	if ( seed > 0 ){
		attr.prob <- stats::runif( L , 0 , 10/ L )
		attr.prob <- attr.prob / sum( attr.prob )
				}
            
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

    # constraints for guessing and slipping parameters
    if ( is.null( constraint.slip ) == FALSE ){  slip[ constraint.slip[,1] ] <- constraint.slip[,2] }
    if ( is.null( constraint.guess ) == FALSE ){  guess[ constraint.guess[,1] ] <- constraint.guess[,2] }
            
    
    iter <- 1 # Iteration number
    likediff <- 1 # Difference in likelihood estimates
    loglike <- 0 # init for log-Likelihood
    
    # init value for maximum parameter change in likelihood maximization
    dev.change <- max.par.change <- 1000
    
    # calculate for each item how many attributes are necessary for solving the items
    # according to the specified DINA or DINO rule                                  
    comp <- ( rowSums(q.matrix)  )*( rule=="DINA")   + 1* ( rule == "DINO" )               
	# need number of components for I.lj ...
	compL <- outer(  comp, rep(1,L) )
	attrpatt.qmatr <- t(( attr.patt %*% t(q.matrix) ) )                      
	
	# compute latent response
    latresp <- apply( attr.patt, 1, FUN = function(attr.patt.ll){
		# form attribute pattern (Lx1) in matrix (J x 1 )
		attr.patt.ll <- outer( rep(1,J), attr.patt.ll )
		ind <- 1*(rowSums(q.matrix * attr.patt.ll)  >= comp )   	
		ind } )

	# response patterns
	cmresp <- colMeans( resp.patt )
    some.missings <- if( mean(cmresp) < 1){ TRUE } else { FALSE }
	
    # calculations for expected counts
	# this matrix is needed for computing R.lj
	ipr <- item.patt.split * item.patt.freq*resp.patt
#	ips.ipf.rp <- item.patt.freq * resp.patt*item.patt.split
#	ips.ipf.rp <- ips.ipf.rp[ , rep( 1:J , L ) ]	
#	ipf.rp <- item.patt.freq * resp.patt
#	ipf.rp <- ipf.rp[ , rep( 1:J , L ) ]	
	# response indicator list
    resp.ind.list <- list( 1:J )
	for (i in 1:J){ resp.ind.list[[i]] <- which( resp.patt[,i] == 1)  }
	# parameter history
	if ( param.history ){
		likelihood.history <- rep(NA, maxit )
		slip.history <- guess.history <- matrix( NA , nrow=maxit , ncol=ncol(data) )
				}
	
################################################################################
# BEGIN OF THE ITERATION LOOP                                                  #
################################################################################
    
    while ( iter <= maxit & 
			( ( max.par.change > conv.crit ) | ( dev.change > dev.crit ) )
								){

								
################################################################################
# STEP I:                                                                      #
# calculate P(X_i | alpha_l):                                                  # 
# probability of each item response pattern given an attribute pattern         #
################################################################################

 z0 <- Sys.time()                    
				
	slipM <- matrix( slip , nrow= nrow(latresp) , ncol=ncol(latresp))
	guessM <- matrix( guess , nrow= nrow(latresp) , ncol=ncol(latresp))
    pj <- (1 - slipM )*latresp + guessM * ( 1 - latresp )
	pjM <- array( NA , dim=c(J,2,L) )
	pjM[,1,] <- 1 - pj
	pjM[,2,] <- pj
	h1 <- matrix( 1 , nrow=IP , ncol=L )

# cat("probs") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1		
	
    res.hwt <- calc_posterior.v2(rprobs= pjM , gwt=h1 , resp=item.patt.split , 
								 nitems= J , 
                                 resp.ind.list=resp.ind.list , normalization=FALSE , 
                                 thetasamp.density= NULL , snodes=0 )	
    p.xi.aj <- res.hwt[["hwt"]]  	
	# set likelihood for skill classes with zero probability to zero
	if ( ! is.null(zeroprob.skillclasses) ){
		p.xi.aj[ , zeroprob.skillclasses ] <- 0
								}

# cat("calc.like") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1										
	   
################################################################################
# STEP II:                                                                     #
# calculate P(  \alpha_l | X_i ):                                              #
# posterior probability of each attribute pattern given the item response pattern
################################################################################
                                           
    # posterior probabilities  P( \alpha_l | X_i ) 
    p.aj.xi <- outer( rep(1,IP), attr.prob ) * p.xi.aj 
    
	attr.prob0 <- attr.prob
	
	if ( ! is.null( zeroprob.skillclasses ) ){
		p.aj.xi[ , zeroprob.skillclasses ] <- 0 }
	
    p.aj.xi <- p.aj.xi / rowSums( p.aj.xi )

    # calculate marginal probability P(\alpha_l) for attribute alpha_l
    attr.prob <- colSums( p.aj.xi * item.patt.freq / I )

# cat("calc.post") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1											
	
	#***
	# include structural constraints for high-dimensional skill vectors here!
	#***

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

      
#    I.lj <- sapply( 1:L, FUN = function(ll){ # skill pattern ll
    # I_{lj} = \sum_i [ N_i * P( \alpha_l | X_i ) ]
    # I_{lj} is the sum of the product of item reponse pattern and the conditional
    # probability of being in attribute class l given item response pattern
#    colSums( item.patt.freq * resp.patt * p.aj.xi[,ll] )    
#                            } )

				
		R.lj <- I.lj <- matrix( 0 , nrow=J , ncol=L )
			
		if ( some.missings ){
#			I.lj <- t( item.patt.freq*resp.patt	) %*% p.aj.xi
			I.lj <- crossprod( item.patt.freq*resp.patt	, p.aj.xi )
					} else {
			I.lj <- matrix( t( item.patt.freq ) %*% p.aj.xi , nrow=J , 
							ncol=L , byrow=TRUE )
					}

					
#		R.lj <- t(ipr) %*% p.aj.xi	
		# ipr <- item.patt.split * item.patt.freq*resp.patt	
		R.lj <- crossprod(ipr , p.aj.xi	)
		colnames(I.lj) <- colnames(R.lj) <- attr.patt.c
		rownames(I.lj) <- rownames(R.lj) <- colnames(data)
			

# cat("calc R and I ") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1											
			
################################################################################
# STEP IV:                                                                     #
# calculate R_j0, I_j0, I_j1, R_j1                                             #
# R_{j1} ... expected frequency of students which correctly solve item j and   #
#               possess all necessary attributes for this item                 #
# I_{j1} ... expected frequency of students which correctly solve the item     #
# I_{j0}, R_{j0} ... expected frequencies of students which incorrectly        #
#                      solve the item                                          #
################################################################################

    I.j0 <- rowSums( ( attrpatt.qmatr  < compL  ) * I.lj ) 
    I.j1 <- rowSums( ( attrpatt.qmatr  >= compL ) * I.lj )             
    R.j0 <- rowSums( ( attrpatt.qmatr  <  compL ) * R.lj )
    R.j1 <- rowSums( ( attrpatt.qmatr  >= compL )  * R.lj )

# cat("calc suff stats m step") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1									
	
################################################################################
# STEP V:                                                                      #
# M-Step: update slipping and guessing parameters.                             #
# The guessing parameter 'guess.new' can be claculated as R.j0 / I.j0          #
#   -> correct solution for item j if not all necessary attributes are possessed
################################################################################

    I.j0[ I.j0 == 0] <- 0.05
    I.j1[ I.j1 == 0] <- 0.05           
    guess.new <- R.j0 / I.j0
    slip.new <- ( I.j1 - R.j1 ) / I.j1
    # equal guessing and slipping parameter
    if (guess.equal ){ guess.new <- rep( sum( R.j0 ) / sum( I.j0) , J ) }
    if (slip.equal ){ slip.new <- rep( sum( I.j1 - R.j1 ) / sum( I.j1) , J ) }

	# overrelaxed convergence
	if ( wgt.overrelax > 0 ){
		if ( wgtest.overrelax & ( iter > 2) ){
			l1 <- c( guess.new , slip.new )
			l2 <- c( guess , slip )
			l3 <- c( guess.old , slip.old )
			lambda <- sum( abs( l1 - l2 ) ) / sum( abs( l3 - l2 ) ) 
			wgt.overrelax <- lambda / ( 2 - lambda )
			wgt.overrelax[ wgt.overrelax <= 0 ] <- 10^(-5)
			wgt.overrelax[ wgt.overrelax >= .98 ] <- .98
					}
		guess.new <- guess.new + wgt.overrelax * ( guess.new - guess )
		slip.new <- slip.new + wgt.overrelax * ( slip.new - slip )
		guess.new[ guess.new < 0 ] <- 0
		slip.new[ slip.new < 0 ] <- 0
		if ( wgtest.overrelax ){ 	
			slip.old <- slip
			guess.old <- guess
					}
			}
    # constrained slipping and guessing parameter
    if ( is.null( constraint.slip ) == F ){  slip.new[ constraint.slip[,1] ] <- constraint.slip[,2] }
    if ( is.null( constraint.guess ) == F ){  guess.new[ constraint.guess[,1] ] <- constraint.guess[,2] }

	
	#*** adjustment parameters
	if ( guess.min > 0 ){
		guess.new <- ifelse( guess.new < guess.min , guess.min , guess.new )
				}
	if ( guess.max < 1 ){
		guess.new <- ifelse( guess.new < guess.max , guess.max , guess.new )
				}
	if ( slip.min > 0 ){
		slip.new <- ifelse( slip.new < slip.min , slip.min , slip.new )
				}
	if ( slip.max < 1 ){
		slip.new <- ifelse( slip.new < slip.max , slip.max , slip.new )
				}				
				

    # calculate the updated likelihood                                               
    like.new <- sum( log( rowSums( p.xi.aj * outer( rep(1,IP), attr.prob )  ) + 10^(-300) ) * item.patt.freq )
    likediff <- abs( loglike - like.new )
    loglike <- like.new
    dev.change <- abs( likediff / loglike )
	
    # maximum parameter change
    max.par.change <- max( abs( guess.new - guess ), abs( slip.new - slip ) ,
				abs( attr.prob - attr.prob0 ) )

# cat("calc pars and like") ; z1 <- Sys.time(); print(z1-z0) ; z0 <- z1											
				
    # define estimates which are updated in this iteration
    guess <- guess.new
    slip <- slip.new
    if (progress) {  
		cat( "Iter. ",iter, " :", 
			substring( paste(Sys.time()), first=11 ), ", ", " loglike = ", 
			round( like.new, 7), 
			" / max. param. ch. : ", round( max.par.change, 6) , 
			" / relative deviance change : ", round( dev.change, 6)  			
			) 
			
			if( wgt.overrelax > 0){ 
				cat(" / overrelax. parameter =" , round( wgt.overrelax , 4 ))
					}
			cat("\n", sep="") 
				}
	if ( param.history ){ 				
			likelihood.history[iter] <- like.new
			slip.history[iter,] <- slip
			guess.history[iter,] <- guess		
						}
					
    flush.console() # Output is flushing on the console
    iter <- iter + 1 # new iteration number                                    
}

################################################################################
# END OF THE ITERATION LOOP                                                    #
################################################################################

#    if(any(guess > 1 - slip)) 
#        warning(paste("Negative item discrimination index in item",ifelse(length(which(guess > 1 - slip))>1,"s "," "),
#    paste(which(guess > 1 - slip), collapse=", "),". See Help-files to fix this problem.", sep=""))

	# set likelihood for skill classes with zero probability to zero
	if ( ! is.null(zeroprob.skillclasses) ){
		p.xi.aj[ , zeroprob.skillclasses ] <- 0
								}
    # calculate posterior probability for each attribute pattern
    pattern <- data.frame( freq = round(as.numeric(item.patt[,-1]),3),
                    mle.est = attr.patt.c[ max.col( p.xi.aj ) ], 
                    mle.post = rowMaxs( p.xi.aj ) / rowSums( p.xi.aj ), 
                    map.est = attr.patt.c[ max.col( p.aj.xi ) ], 
                    map.post = rowMaxs( p.aj.xi )
                    )

    # calculate posterior probabilities for all skills separately
    attr.postprob <- p.aj.xi %*% attr.patt
    colnames( attr.postprob ) <- paste("post.attr",1:K, sep="")
    pattern <- cbind( pattern,  attr.postprob )

################################################################################
# estimation of the standard errors for slipping and guessing parameters       #
# guessing parameters (DINA and DINO)                                          #
# NOTE: In this calculation of standard errors, sample weights were not        #
#   taken into account. Use for example the bootstrap to do some inference.    #
################################################################################

    se.guess <- sapply( 1:J, FUN = function(jj){
                indA0.jj <- rowSums( q.matrix[jj,] * attr.patt ) < comp[jj]
                l1 <- rowSums( p.aj.xi * outer( resp.patt[,jj], rep(1,L) ) *
                        ( outer( item.patt.split[,jj ], rep(1,L) ) - guess[ jj ] ) * outer( rep(1,IP), indA0.jj )
                                    )
                ( sum( item.patt.freq * l1^2  ) / guess[jj]^2 / (1-guess[jj])^2 )^(-.5)
                    } )
    # equal guessing parameter
    if ( guess.equal ){ se.guess <- rep( sqrt( 1 / sum( 1/ se.guess^2  ) ) , J ) } 
    
    
    
    # constrained guessing parameter
    if ( ! is.null( constraint.guess )  ){  se.guess[ constraint.guess[,1] ] <- NA }
    guess <- data.frame( est = guess, se = se.guess                 )      
    
    # slipping parameter (DINA and DINO)
    se.slip <- sapply( 1:J, FUN = function(jj){
                indA0.jj <- rowSums( q.matrix[jj,] * attr.patt ) >= comp[jj]
                l1 <- rowSums( p.aj.xi * outer( resp.patt[,jj], rep(1,L) ) *
                        ( outer( item.patt.split[,jj ], rep(1,L) ) - 1 + slip[ jj ] ) * outer( rep(1,IP), indA0.jj )
                                    )
                ( sum( item.patt.freq * l1^2  ) / slip[jj]^2 / (1-slip[jj])^2 )^(-.5)
                    } )

    # equal slipping parameter
    if ( slip.equal ){ se.slip <- rep( sqrt( 1 / sum( 1/ se.slip^2  ) ) , J ) }
    
    # constrained slipping parameter
    if ( is.null( constraint.slip ) == F ){  se.slip[ constraint.slip[,1] ] <- NA }        
    slip <- data.frame( est = slip, se = se.slip )
    
    # attribute pattern
    attr.prob <- matrix( attr.prob, ncol=1)
    colnames( attr.prob ) <- "class.prob"
    rownames( attr.prob ) <- attr.patt.c
    
    # pattern for seperate skills
    skill.patt <- matrix(apply( matrix( rep(  attr.prob, K ), ncol=K) * attr.patt, 2, sum ),ncol=1)
#    rownames(skill.patt) <- paste("Skill_", colnames(q.matrix),sep="")
    rownames(skill.patt) <- colnames(q.matrix)
    colnames(skill.patt) <- "skill.prob" 
    
    # calculation of the AIC und BIC        
    bb <- 0
    if ( is.null( constraint.guess ) == F ){  bb <- bb + nrow(constraint.guess) }
    if ( is.null( constraint.slip ) == F ){  bb <- bb + nrow(constraint.slip) }
	if( guess.equal){ bb <- bb + J - 1 }
	if( slip.equal){ bb <- bb + J - 1 }
	# collect number of parameters
	pars <- data.frame( "itempars" = 2*J - bb )
	# number of skill classes
	pars$skillpars <- L - 1 - length(  zeroprob.skillclasses )
	Np <- pars$itempars + pars$skillpars
    aic <- -2*loglike + 2*( Np )
    bic <- -2*loglike + ( Np )*log(I)

    rownames( p.aj.xi ) <- rownames( pattern ) # output rownames posterior probabilities    
    pattern <- data.frame(pattern) # convert pattern to numeric format
    for (vv in seq(1,ncol(pattern))[ -c(2,4) ] ){
                    pattern[,vv ] <- as.numeric( paste( pattern[,vv] ) ) }
    
    # subject pattern
	# changed item.patt.subj$pattern.index (ARb 2012-06-05)
    item.patt.subj <- data.frame( "case" = 1:(nrow(data) ), 
                                   "pattern" = item.patt.subj, 
                                    "pattern.index" = match( item.patt.subj, item.patt[,1] )
                                            )
    
    # attribute pattern (expected frequencies)
    attr.prob <- data.frame( attr.prob )
    attr.prob$class.expfreq <-  attr.prob[,1] * nrow(data) 

	#*****
	# modify output (ARb 2012-06-05)
	pattern <- pattern[ item.patt.subj$pattern.index , ]	
	pattern[,1] <- paste( item.patt.subj$pattern )
	colnames(pattern)[1] <- "pattern"
	p.aj.xi__ <- p.aj.xi
	p.aj.xi <- p.aj.xi[ item.patt.subj$pattern.index , ]
	rownames(p.aj.xi) <- pattern$pattern
	colnames(p.aj.xi) <- rownames(attr.prob)
	p.xi.aj <- p.xi.aj[ item.patt.subj$pattern.index , ]
	rownames(p.xi.aj) <- pattern$pattern
	colnames(p.xi.aj) <- colnames(p.aj.xi)

	#########################################
	# item fit [ items , theta , categories ] 
	# # n.ik [ 1:TP , 1:I , 1:(K+1) , 1:G ]
	n.ik <- array( 0 , dim=c(L , J , 2 , 1 ) )
	n.ik[  , , 2 , 1 ] <- t(R.lj)
	n.ik[  , , 1 , 1 ] <- t(I.lj-R.lj)
	pi.k <- array( 0 , dim=c(L,1) )
	pi.k[,1] <- attr.prob$class.prob
	probs <- aperm( pjM , c(3,1,2) )
	itemfit.rmsea <- itemfit.rmsea( n.ik , pi.k , probs )$rmsea
	
	#*****
	
    datfr <-  data.frame( round( cbind( guess , slip  ), 3 ) )
    colnames(datfr) <- c("guess", "se.guess", "slip", "se.slip" )
    rownames(datfr) <- colnames( dat.items )
    datfr <- data.frame( "type" = rule, datfr )
	datfr$rmsea <- itemfit.rmsea
	names(itemfit.rmsea) <- colnames(data)	
    s2 <- Sys.time()
	
	#******
	# parameter table for din object
	res.partable <- din.partable( guess , slip , attribute.patt = attr.prob , data=data ,
					rule= paste0(datfr$type) ,
					guess.equal , slip.equal , constraint.guess , constraint.slip ,
					zeroprob.skillclasses ,
					attribute.patt.splitted = attr.patt  )	
					
    partable <- res.partable$partable
	vcov.derived <- res.partable$vcov.derived	
	
    if (progress){ cat("---------------------------------------------------------------------------------\n") }
	# coefficients
	p1 <- partable[ partable$free , ]
	p1 <- p1[ ! duplicated(p1$parindex) , ]	
	p11 <- p1$value
	names(p11) <- p1$parnames

	
    res <- list( coef=p11 , 
				item = datfr, guess = guess, slip = slip  , 
				"IDI" = round(1 - guess[,1] - slip[,1] ,4)  ,
				"itemfit.rmsea" = itemfit.rmsea , 
				"mean.rmsea" = mean(itemfit.rmsea) , 
				loglike = loglike, AIC = aic, BIC = bic, 
				"Npars" = pars ,
                 posterior = p.aj.xi, "like" = p.xi.aj, 
				 "data" = data, "q.matrix" = q.matrix,
                 pattern = pattern , attribute.patt = attr.prob, skill.patt = skill.patt,
                 "subj.pattern" = item.patt.subj, "attribute.patt.splitted" = attr.patt, 
                 "display" = disp, "item.patt.split" = item.patt.split, 
                 "item.patt.freq" = item.patt.freq, "model.type" = r1 , 
				 "rule" = rule , "zeroprob.skillclasses" = zeroprob.skillclasses , 
				 "weights" = weights , "pjk" = pjM , "I" = I , 
				 "I.lj"=I.lj , "R.lj" = R.lj , "partable" = partable ,
				 "vcov.derived" = vcov.derived , 				 
				 "seed" = seed , 
				 "start.analysis" = s1 , "end.analysis" = s2 ,
				 "iter" = iter 	,
				 "converged" = iter < maxit 
					) 
	res$timediff <- s2 - s1
	if (progress){ print(s2-s1) }
	if (param.history){
		param.history <- list( "likelihood.history" = likelihood.history , 
				"slip.history" = slip.history , 
				"guess.history" = guess.history )
		res$param.history <- param.history
					}					
	# control parameters
	control <- list( q.matrix=q.matrix , skillclasses = skillclasses , conv.crit = conv.crit , 
					dev.crit = dev.crit , maxit = maxit ,
                    constraint.guess = constraint.guess , constraint.slip = constraint.slip ,
                    guess.init = guess.init , slip.init = slip.init ,
                    guess.equal = guess.equal , slip.equal = slip.equal , 
					zeroprob.skillclasses = zeroprob.skillclasses , 
                    weights = weights ,  rule = rule , 
					wgt.overrelax = wgt.overrelax , wgtest.overrelax = wgtest.overrelax ,
					latresp=latresp , resp.ind.list=resp.ind.list 
						)	
    res$control <- control		
    res$call <- cl	
    class(res) <- "din"					
    return(res)
}
