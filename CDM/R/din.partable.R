
##############################################
# parameter table for the din model
din.partable <- function( guess ,  slip , attribute.patt , data , rule ,
		guess.equal , slip.equal , constraint.guess , constraint.slip ,
		zeroprob.skillclasses , attribute.patt.splitted ){
	# parameters
	J <- nrow(guess)
	L <- nrow(attribute.patt)
	items <- colnames(data)
	K <- ncol(attribute.patt.splitted)
	
	#**************
	# create parameter table
	partable <- data.frame( "partype" = c( rep( c("guess","slip") , J) , 
						  rep("probs" , L)  , rep("margprobs" , K ) )
							) 
	partable$parindex <- c( 1:J , J + 1:J , 2*J + 1:(L-1) , 0 , rep( 0 , K ))
	partable$item <- c( rep(1:J , each=2 ) , rep(0,L+K) )
	partable$item.name <- c( rep(colnames(data),each=2) , rep("",L+K) )
	partable$skillclass <- c( rep(0,2*J) , 1:(L) , rep(0,K))
	partable$varyindex <- c( rep(1:J , each=2 ) , 1:(L-1) , 0 , rep(0,K))
	m1 <- paste0( rep( items , each=2 ) , "_" , c( "guess" , "slip" ) )
	partable$parnames <- c( m1 , paste0( "prob_class" , 1:L ) ,
				paste0( "prob_skill" , 1:K ))
	dfr <- cbind( guess$est , slip$est )
	dfr <- matrix( t(dfr) , nrow=1 , byrow=TRUE )
	
	# marginal skill probabilities
	margskills <- colSums( attribute.patt.splitted * attribute.patt[,1] )

	partable$value <- c( dfr[1,] , attribute.patt[,1] , margskills )
	partable$fixed <- FALSE
	partable$free <- partable$parindex > 0
	partable$rule <- c( rep(rule,2) , rep("",L+K) )
	partable$totindex <- 1:(nrow(partable))
	
	#**************************
	# include item parameter constraints
	if (guess.equal){
		p1 <- partable[ partable$partype == "guess" , ]
		partable[ p1$totindex , "parindex" ] <- p1$parindex[1]
		partable[ p1$totindex , "parnames" ] <- "all_guess"		
				}
	if (slip.equal){
		p1 <- partable[ partable$partype == "slip" , ]
		partable[ p1$totindex , "parindex" ] <- p1$parindex[1]	
		partable[ p1$totindex , "parnames" ] <- "all_slip"			
				}				
	if ( ! is.null(constraint.slip) ){
		p1 <- partable[ ( partable$partype == "slip" ) &
						( partable$item %in% constraint.slip[,1] ) , ]
		partable[ p1$totindex , "fixed" ] <- TRUE				
		partable[ p1$totindex , "free" ] <- FALSE
		partable[ p1$totindex , "parindex" ] <- 0
						}
	if ( ! is.null(constraint.guess) ){
		p1 <- partable[ ( partable$partype == "guess" ) &
						( partable$item %in% constraint.guess[,1] ) , ]
		partable[ p1$totindex , "fixed" ] <- TRUE				
		partable[ p1$totindex , "free" ] <- FALSE
		partable[ p1$totindex , "parindex" ] <- 0
						}	
	if ( ! is.null(zeroprob.skillclasses) ){
		p1 <- partable[ ( partable$partype == "probs" ) &
						( partable$skillclass %in% zeroprob.skillclasses ) , ]
		partable[ p1$totindex , "fixed" ] <- TRUE				
		partable[ p1$totindex , "free" ] <- FALSE
		partable[ p1$totindex , "parindex" ] <- 0
						}																				
	#*********************************
	# include parameter transformation matrix
	
	estpars <- unique( partable[ partable$parindex > 0 , "parnames" ] )
	allpars <- unique( partable$parnames )
	MP <- length(allpars)
	FP <- length(estpars)
	A <- matrix( 0 , nrow=MP , ncol=FP) 
	rownames(A) <- allpars
	colnames(A) <- estpars
	
	# free parameters
	a1 <- match( estpars , allpars )
	A[ cbind( a1 , 1:FP ) ] <- 1
	# probabilities of last class
	probs_names <- partable[ partable$partype == "probs" , "parnames" ]
	v1 <- probs_names[ length(probs_names ) ]
	v2 <- intersect( setdiff( probs_names , v1 ) , estpars )
	A[ v1 , v2 ] <- - 1	
	# marginal skill probabilities
    rownames(attribute.patt.splitted) <- probs_names
	colnames(attribute.patt.splitted) <- 
			partable[ partable$partype == "margprobs" , "parnames" ]
	attribute.patt.splitted <- ( attribute.patt.splitted - 1 ) 
	a1 <- t(attribute.patt.splitted)
	a1 <- a1[ , intersect( estpars , colnames(a1) ) ]
	A[ rownames(a1) , colnames(a1) ] <- a1
	#*********************************
	# introduce new parameter index
	partable$parindex <- match( partable$parindex , sort(unique(partable$parindex) ) ) - 1
	res <- list( "partable" = partable , "vcov.derived" = list("A"=A) )
	return(res)
		}
###########################################################		
	
	
	