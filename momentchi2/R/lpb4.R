#need error message to make sure number of terms greater than 3?
# if () warning("number of terms less than three: stop")
#or stop("number of terms les than three - migth give  NaNs")

#call this 'lpbcoeff' - will work for any p, with default p=4
#maybe will force p=4
#only make lpb4 for first 4 (or 8?) cumulants.
#need a function to compute coeffs for lpb to get started (at least three terms)
#then can update after that - include in next version of code?


#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#LPB method below, SW, HBE and WF at the end 

#Steps for computing cdf using (Lindsay, 2000) method

#step 0.1: Obtain coefficients d_i for H = sum_i^n d_i w_i^2
#		 These determine the distribution

#step 0.2: Decide on p, the number of support points to use for method
#		   The more support points, the more accurate (and more computationally intensive)	

#step 1: Determine/compute the moments/cumulants m_1(H), ... m_2p(H)
#			First compute cumulants - sums of powers of coefficients
#			cf: (Wood, 1989)
#			Then use the cumulants to compute moments (using a recursive formula):
#			mu_n = kappa_n + \sum_{m=1}^{n-1} (n-1)choose(m-1) kappa_n mu_{n-m}


#step 2.1: generate the matrices delta_i(x)

#step 2.2: Find lambdatilde_1, the unique root of det (delta_1(x))
#			This does not require bisection - just a rational expression in
#			terms of the moments

#step 3:	Use bisection method (R uniroot) to find lambdatilde_2 
#			in [0, lambdatilde_1)		
#			Find lambdatilde_i+1 in [0, lambdatilde_i) for i = 2,3,...p
#
#			End result: we have lambdatilde_p

#step 4: should have this method from step 2 already, but compute
#			deltastar_i(lambdatilde_p) for i =1,2,...2p-1

#step 5.1: use the deltastar_i(lambdatilde_p) from step 4 to generate
#			matrix Stilde(lambdatilde_p, t)
#
#step 5.2:	Then need to diagonalise/use linear algerba trick in paper
#			to get polynomial coefficients (from det) in "coeff_vec"
#
#step 5.3	use Re(polyroot(coeff_vec)) to obtain roots of polynomial
#			denoted mu_vec = (mu_1, ..., mu_p)	

#step 6:	Generate Vandermonde matrix using mu_vec
#			and vector using deltastar_i's, to solve for
#			pi_vec = (pi_1, ..., pi_p)

#step 7: 	? compute the linear combination (using pi_vec)
#			of the i gamma cdfs using parameters lambdatilde_p and mu_i
#
#			This should be the final answer


#at the moment, only working when length(coeffvec)==2 exactly; less (1 term) doesn't work, and 3 or more terms 
#doesn't work either. 
#So, logically, the error must be related to coeffvec.
#However, the only method that calls coeffvec seems to be correct (computing the cumulants)
#and computing the moments is also correct.

# quantile_vec - the vector of quantile value used for computing the cdf,
# coeffvec - the coefficients of the weighted sum of chi-squared terms,
# p - the p "support points" used for computing the cdf; the larger the p, the more accurate (and computationally
#	  intensive) the algorithm is.
# note: in Lindsay et al. (2000), a value of p=4 is recommended. p=2 corresponds to SW.
# default p set to 4 


#' Lindsay-Pilla-Basak method
#'
#' Computes the cdf of a positively-weighted sum of chi-squared random variables with the Lindsay-Pilla-Basak (LPB4) method using four support points. Note that the coefficient vector must be of length at least four.
#' @inheritParams sw
#' @details Note that the coefficient vector must of length at least four. In some cases when the coefficient vector was of length two or three, the algorithm would be unable to find roots of a particular equation during an intermediate step, and so the algorithm would produce \code{NaN}s. If the coefficient vector is of length less than four, the Hall-Buckley-Eagleson method is used (and a warning is displayed).
#' @keywords distribution
#' @references
#' \itemize{
#'  \item B. G. Lindsay, R. S. Pilla, and P. Basak. Moment-based approximations of distributions using mixtures: Theory and applications. \emph{Annals of the Institute of Statistical Mathematics}, 52(2):215-230, 2000.
#' }
#' @export
#' @examples
#' #Examples taken from Table 18.6 in N. L. Johnson, S. Kotz, N. Balakrishnan. 
#' #Continuous Univariate Distributions, Volume 1, John Wiley & Sons, 1994.
#'
#' lpb4(c(1.5, 1.5, 0.5, 0.5), 10.203)            # should give value close to 0.95
#' lpb4(coeff=c(1.5, 1.5, 0.5, 0.5), x=10.203)    # specifying parameters
#' lpb4(c(1.5, 1.5, 0.5, 0.5), c(0.627, 10.203))  # x is a vector, output approx 0.05, 0.95
#' lpb4(c(0.5, 0.3, 0.2), 2.708)                  # length(coeff) < 4, warning, uses hbe()

lpb4 <- function(coeff, x){
    if ( (missing(x)) || (missing(coeff)) )
        stop("missing an argument - need to specify \"coeff\" and \"x\"")

    if (checkCoeffsArePositiveError(coeff))
        stop(getCoeffError(coeff))

    if (checkXvaluesArePositiveError(x))
        stop(getXvaluesError(x))

	#check if there is less than 4 elements - if so, stop
	if (length(coeff) < 4){
        #stop("warning: less than four coefficients - LPB4 method may return NaN, so stopping.")
        warning("Less than four coefficients - LPB4 method may return NaN: running hbe instead.")
        return( hbe(coeff, x) )
	}
    #end of error checking
	
	#----------------------------------------------------------------#
	#step 0: decide on parameters for distribution and support points p
	#specified to be 4 for this version of the function
	p <- 4
	
	#----------------------------------------------------------------#
	#step 1: Determine/compute the moments m_1(H), ... m_2p(H)
	
	#compute the first 2p moments for Q = sum coeff chi-squared	
	moment_vec <- get_weighted_sum_of_chi_squared_moments(coeff, p)
	
	#----------------------------------------------------------------#
	#Step 2.1: generate matrices delta_i(x)
	#functions created:
	#deltaNmat_applied
	#and
	#det_deltamat_n
	
	#Step 2.2: get lambdatilde_1 - this method is exact (no bisection), solves determinant equation
	lambdatilde_1 <- get_lambdatilde_1(moment_vec[1], moment_vec[2])
	
	#----------------------------------------------------------------#
	#Step 3:	Use bisection method (R uniroot) to find lambdatilde_2 
	#and others up to lambdatilde_p, for tol=bisect_tol
	#all we need is the final lambdatilde_p, not the intermediate values lambdatilde_2, lambdatilde_3, etc
	bisect_tol <- 1e-9
	lambdatilde_p <- get_lambdatilde_p(lambdatilde_1, p, moment_vec, bisect_tol)
	
	#----------------------------------------------------------------#
	#Step 4:
	#Calculate delta_star_lambda_p
	#can already do this using methods in Step 2.1 
	#----------------------------------------------------------------#
	
	#----------------------------------------------------------------#
	#Step 5:
	
	#Step 5.1: use the deltastar_i(lambdatilde_p) from Step 4 to generate
	#			M_p, which will be used to create matrix Stilde(lambdatilde_p, t)
	M_p <- deltaNmat_applied(lambdatilde_p, moment_vec, p)
	
	#Step 5.2:	Compute polynomial coefficients of the modified M_p matrix (as in paper).
	mu_poly_coeff_vec <- get_Stilde_polynomial_coefficients(M_p)
	
	#step 5.3	use Re(polyroot(coeff_vec)) to obtain roots of polynomial
	#			denoted mu_vec = (mu_1, ..., mu_p)	
	mu_roots <- Re(polyroot(mu_poly_coeff_vec))
	
	#----------------------------------------------------------------#
	
	#Step 6:	Generate Vandermonde matrix using mu_vec
	#			and vector using deltastar_i's, to solve for
	#			pi_vec = (pi_1, ..., pi_p)
	pi_vec <- generate_and_solve_VDM_system(M_p, mu_roots)
	
	#----------------------------------------------------------------#
	
	#Step 7: 	Compute the linear combination (using pi_vec)
	#			of the i gamma cdfs using parameters lambdatilde_p and mu_i 
	#			(but need to create scale/shape parameters carefully)
	#	
	#			This is the final answer
	
	mixed_p_val_vec <- get_mixed_p_val_vec(x, mu_roots, pi_vec, lambdatilde_p)
	
	#We have our final answer, and so return it
	return(mixed_p_val_vec)
	#cat("done\n")
}



#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#Methods for computing steps

#Methods for computing steps

#Methods for computing steps


#-----------------------------------------------------------------------------#
#Step 1:

#hides the computation of the cumulants, by just talking about moments
get_weighted_sum_of_chi_squared_moments <- function(coeffvec, p){
	cumulant_vec <- get_cumulant_vec_vectorised(coeffvec, p)
	moment_vec <- get_moments_from_cumulants(cumulant_vec)
	return (moment_vec)
}


#get the cumulants kappa_1, kappa_2, ..., kappa_2p
get_cumulant_vec_vectorised <- function(coeffvec, p){
	index <- c(1:(2*p))
	#cumulants <- 2^(v-1) * factorial(v-1)
	#now cumulants are multiplied by sum_of_powers of coeffvec
	#moment_vec <- cumulants * vapply(X=v, FUN=sum_of_powers, FUN.VALUE=rep(0,1), v=coeffvec)	
	cumulant_vec <- 2^(index-1) * factorial(index-1) * vapply(X=index, FUN=sum_of_powers, FUN.VALUE=rep(0,1), v=coeffvec)
	return(cumulant_vec)
}


#returns the sum of the elements raised to a power
sum_of_powers <- function(index, v){
	return(sum(v^index))
}


#get the moment vector from the cumulant vector 
#have removed one for loop (vectorised), but can't rmeove the other one
get_moments_from_cumulants <- function(cumulant_vec){
	#start off by assigning it to cumulant_vec, since moment[n] = cumulant[n] + {other stuff}
	moment_vec <- cumulant_vec
	#check if more than 1 moment required
	if (length(moment_vec)>1){
		#can't get rid of this for loop, since updates depend on previous moments
		for (n in 2:length(moment_vec)){
			#can vectorise this part, I think
			moment_vec[n] <- moment_vec[n] + update_moment_from_lower_moments_and_cumulants(n, moment_vec, cumulant_vec) 					
		}#end of for
	}#end of if
	return(moment_vec) 
}#end of get_moments


#returns the sum of the additional terms/lower products of moments and cumulants
#used in the computation of moments
update_moment_from_lower_moments_and_cumulants <- function(n, moment_vec, cumulant_vec){
	m <- c(1:(n-1))
	sum_of_additional_terms <- sum(choose(n-1, m-1) * cumulant_vec[m] * moment_vec[n-m])
	return(sum_of_additional_terms)
} 


#-----------------------------------------------------------------------------#
#Step 2.1: get lambdatilde_1

#no need to use bisection method - can get lambdatilde_1 directly 
get_lambdatilde_1 <- function(m1, m2){
	return ( m2/(m1^2) - 1 )
}


#-----------------------------------------------------------------------------#
#Step 2.2: generate delta_mat_N and det_delta_mat_N

#compute the delta_N matrix - vectorised using lapply and mapply
deltaNmat_applied <- function(x, m_vec, N){
	Nplus1 <- N+1
	#want moments 0, 1, ..., 2N
	m_vec <- c(1, m_vec[1:(2*N)])
	
	#these will be the coefficients for the x in (1+c_1*x)*(1+c_2*x)*...
	#want coefficients 0, 0, 1, 2, .., 2N-1 - so 2N+1 in total 
	coeff_vec <- c(0, 0:(2*N-1))*x + 1
	#not necessary to initialise, could use length(m_vec) below, but we do it anyway for readability
	prod_x_terms_vec <- rep(0, 2*N+1)
	#this computes the terms involving lambda in a vectorised way
	prod_x_terms_vec <- 1/vapply(c(1:length(prod_x_terms_vec)), FUN=get_partial_products, FUN.VALUE=c(0), vec=coeff_vec)
	
	#going to use mapply over matrix indices i, j
	i_vec <- c(1:Nplus1)
	j_vec <- c(1:Nplus1)
	#not necessary to initialise
	#delta_mat <- matrix(0, Nplus1, Nplus1)
	delta_mat <- mapply( get_index_element, i=i_vec, 
			MoreArgs=list(j=j_vec, vec1 = m_vec, vec2 = prod_x_terms_vec), SIMPLIFY="matrix")	
	return(delta_mat)
}


#get_partial_products gets prod[1:index] 
get_partial_products <- function(index, vec){
	return( prod(vec[1:index]) )
}


#this function in deltaNmat_applied computes the index from i and j, and then returns the appropriate product
#of vec1 and vec2 
#(in deltaNmat_applied, these vectors are the moment vector and the vector of the products of the (1+N*lambda)^(-1) terms)
get_index_element <- function(i, j, vec1, vec2){
	index <- i +j -1
	return(vec1[index] * vec2[index])
}


#Simply uses above matrix generation function
det_deltamat_n <- function(x, m_vec, N){
	return(  det(deltaNmat_applied(x, m_vec, N))  )
}


#-----------------------------------------------------------------------------#
#Step 3: get lambdatilde_p
#uses det_delta_mat_n and uniroot

#get lambdatilde_p by using bisection method repeatedly. 
#Need lambdatilde_1 to start
#Need to use R function uniroot
get_lambdatilde_p <- function(lambdatilde_1, p, moment_vec, bisect_tol){
	lambdatilde_vec <- rep(0, p)
	lambdatilde_vec[1] <- lambdatilde_1 
	bisect_tol <- 1e-9
	
	#check that p>1
	if(p > 1){
		for (i in 2:p){
			root <- uniroot(det_deltamat_n, c(0, lambdatilde_vec[i-1]), m_vec=moment_vec, N=i, tol=bisect_tol)
			lambdatilde_vec[i] <- root$root		
		}#end of for		
	}
	#now distinguish lambdatilde_p
	lambdatilde_p <- lambdatilde_vec[p]
	return(lambdatilde_p)
}


#-----------------------------------------------------------------------------#
#Step 5.2: Compute polynomial coefficients for mu polynomial

#We could use the linear algebra trick described in the Lindsay paper, but want to avoid
#dealing with small eigenvalues. Instead, we simply compute p+1 determinants.
#
#This method replaces last column with the base vectors (0, ..., 0 , 1, 0, ... 0) 
#to compute the coefficients, and so does not need to compute 
#any eigen decomposition, just (p+1) determinants
get_Stilde_polynomial_coefficients <- function(M_p){
	#number of rows, number of coefficients ( ==(p+1) )
	n <- dim(M_p)[1]
	index <- c(1:n)
	mu_poly_coeff_vec <- vapply(X=index, FUN=get_ith_coeff_of_Stilde_poly, FUN.VALUE=rep(0,1), mat=M_p)
	return(mu_poly_coeff_vec)
}


#generate a base vector of all zeros except for 1 in ith position
get_base_vector <- function(n, i){
	base_vec <- rep(0, n)
	base_vec[i] <- 1
	return(base_vec)
} 


#get the ith coefficient by computing determinant of appropriate matrix
get_ith_coeff_of_Stilde_poly <- function(i, mat){
	n <- dim(mat)[1]
	base_vec <- get_base_vector(n, i)
	mat[, n] <- base_vec
	return(det(mat))
}


#-----------------------------------------------------------------------------#
#Step 6:Generate van der monde (VDM) matrix and solve the system VDM * pi_vec = b

#generates the VDM matrix and solves the linear system. 
#uses R's built in solve function - there may be a better VDM routine (as cited in Lindsay)
generate_and_solve_VDM_system <- function(M_p, mu_roots){
	#easiest way to get rhs vector is to just take first column of M_p
	b_vec <- get_VDM_b_vec(M_p)
	#print(b_vec)
	#generate Van der Monde matrix; just powers of mu_roots
	VDM <- generate_van_der_monde(mu_roots)
	#print(VDM)
	
	#cat("pi_vec:\n")
	#use R's solve function to solve the linear system
	#there may be better routines for this, but such an implementation is deferred until later
	#NB: If p is too large (p>10), this can yield an error (claims the matrix is singluar).
	#A tailor-made VDM solver should fix this.
	pi_vec <- solve(VDM, b_vec)
	return(pi_vec)	
}


#simply takes the last column, and removes last element of last column
get_VDM_b_vec <- function(mat){
	b_vec <- mat[, 1]
	b_vec <- b_vec[-length(b_vec)]
	return(b_vec)
}


#generates the van der monde matrix from a vector
generate_van_der_monde <- function(vec){
	p <- length(vec)
	vdm <- matrix(0, p, p)
	for (i in 0:(p-1)){
		vdm[i+1, ] <- vec^i
	}
	return(vdm)
}


#-----------------------------------------------------------------------------#
#Step 7: Here we use mu_vec, pi_vec and lambdatilde_p to compute the composite pgamma values 
#		 and combine them into the ifnal pvalue

#get_mixed_p_val - weight sum of pgammas
#now compute for a vector of quantiles - assume the vector of quantiles is very long,
#while p < 10 (so vectorise over length of quantiles)
get_mixed_p_val_vec <- function(quantile_vec, mu_vec, pi_vec, lambdatilde_p){
	#First compute the composite pvalues
	p <- length(mu_vec)
	
	#For pgamma, we need to specify the shape and scale parameters
	#shape alpha = 1/lambda
	alpha <- 1/lambdatilde_p
	#NB: scale beta = mu/alpha, as per formulation in Lindsay paper
	beta_vec <- mu_vec/alpha
	
	#we could probablu vectorise this, but this is simpler
	#we use the pgamma to compute a vector of pvalues from the vector of quantiles, for a given distribution
	#we then scale this by the appropriate pi_vec value, and add this vector to a 0 vector, and repeat
	#finally, each component of the vector is a pi_vec-scaled sum of pvalues
	partial_pval_vec <- rep(0, length(quantile_vec))
	for (i in 1:p){		
		partial_pval_vec <- partial_pval_vec + pi_vec[i] * pgamma(quantile_vec, shape=alpha, scale = beta_vec[i])		
	}
	return(partial_pval_vec)
}

#computes pgamma of the appropriate gamma function
compute_composite_pgamma <- function(index, qval, shape_val, scale_vec){
	return(  pgamma(qval, shape=shape_val, scale = scale_vec[index])  )
}

