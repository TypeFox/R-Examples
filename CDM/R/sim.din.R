################################################################################
# utility function for data simulation from a CDM model                        #
################################################################################

sim.din <- function( N=0, q.matrix, guess = rep(.2, nrow(q.matrix) ), slip = guess, 
	mean = rep(0, ncol(q.matrix) ), Sigma = diag( ncol(q.matrix) ), rule = "DINA",
	alpha = NULL){

# N: a required numeric value specifying the number of requested response pattern. 
#
# q.matrix: a required binary matrix describing which attributes are required, coded by 1,
#		and which attributes are not required, coded by 0, to master the items, whereas the
#		attributes are in the columns and the items in the rows.
#
# guess: an optional vector of guessing parameters. Default is 0.2 for each item.
#
# slip: an optional vector of guessing parameters. Default is 0.2 for each item.
#
# mean: a numeric vector of length nrow(q.matrix) indicating the latent trait of the respondents. 
#		Default is rep(0, length = nrow(q.matrix)). 
#
# Sigma: a matrix of dimension nrow(q.matrix) times nrow(q.matrix) specifying the covariance matrix. 
#		Default is diag( 1, nrow(q.matrix)).
#
# rule: an optional character string or vector of character strings specifying the model rule that is used. 
#		The character strings must be of "DINA" or "DINO". If a vector of character strings is specified, 
#		implying an itemwise condensation rule, the vector must be of length nrow(q.matrix). The default is the 
#		condensation rule "DINA" for all items.

	# simulated normal variates
	if (N>0){
		normsim <- mvtnorm::rmvnorm( N, mean, Sigma)
		# dichotomous variates
		dichsim <- 1 * ( normsim > 0 )
				}
	if ( ! is.null( alpha) ){ 
			dichsim <- alpha 
			N <- nrow(alpha)
			}

# print(dichsim)
			
    # number of possessed attributes, of those which are necessary for this item
    poss.attr <- dichsim %*% t( q.matrix )

    # calculate for each item how many attributes are necessary for solving the items
	# according to the specified DINA or DINO rule        
	ness.attr <- ( rowSums(q.matrix)  )*( rule=="DINA") + 1* ( rule == "DINO" )  
    # latent response
	eta.pp <- poss.attr >= outer( rep(1,N), ness.attr )
	# simulating responses according DINA rule
    R <- matrix( eta.pp * stats::rbinom( N*nrow(q.matrix), size= 1, prob = 1 - outer( rep(1,N), slip ) ) + 
                                ( 1 - eta.pp) * stats::rbinom( N*nrow(q.matrix), size= 1, prob = outer( rep(1,N), guess ) ), ncol = nrow(q.matrix) )
	colnames(R) <- paste( "I" , substring( 1000 + 1:( nrow(q.matrix) ) , 2) , sep="")								
   
	res <- list( "dat" = R , "alpha" = dichsim ) 
	return(res)
}
 