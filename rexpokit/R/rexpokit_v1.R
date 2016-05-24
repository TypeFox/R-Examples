#' @include fermat.R
#' @include rexpokit-package.R
 
#sourcedir = '/Dropbox/_njm/'
#source3 = '_genericR_v1.R'
#source(paste(sourcedir, source3, sep=""))
#roxygenize()

# Re-source this R code after editing, without reinstalling from scratch:
# sourcedir = "/Dropbox/_njm/__packages/rexpokit_setup/"
# source8 = 'rexpokit_v1.R'
# source(paste(sourcedir, source8, sep=""))



# Original source:
# 
#sourcedir = "/Dropbox/_njm/"
#source8 = '_matrix_utils_v1.R'
#source(paste(sourcedir, source8, sep=""))

# for e.g. calc_loglike
# sourcedir = '/Dropbox/_njm/'
# source3 = '_R_tree_functions_v1.R'
# source(paste(sourcedir, source3, sep=""))


#######################################################
# EXPOKIT-RELATED FUNCTIONS
#######################################################

#' EXPOKIT dgpadm matrix exponentiation on Q matrix
#'
#' This function exponentiates a matrix via the EXPOKIT padm function
#' (designed for small dense matrices) and wrapper function 
#' \code{wrapalldgpadm_} around dmexpv.\cr
#'
#' From EXPOKIT:\cr
#'
#' \code{*     Computes exp(t*H), the matrix exponential of a general matrix in }\cr
#' \code{*     full, using the irreducible rational Pade approximation to the   }\cr
#' \code{*     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ), }\cr
#' \code{*     combined with scaling-and-squaring.                              }\cr
#'
#' If \code{Qmat} is NULL (default), a default matrix is input.\cr
#'
#' @param Qmat an input Q transition matrix
#' @param t one or more time values to exponentiate by
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' @return \code{tmpoutmat} the output matrix. \code{wrapalldmexpv_} produces
#' additional output relating to accuracy of the output matrix etc.; these can be
#' obtained by a direct call of wrapalldmexpv_.
#' @seealso \code{\link{mat2coo}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples	# Example:
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' 
#' # Exponentiate each with EXPOKIT's dgpadm (good for small dense matrices)
#' for (t in tvals)
#' 	{
#' 	Pmat = expokit_dgpadm_Qmat(Qmat=Qmat, t=t, transpose_needed=TRUE)
#' 	cat("\n\nTime=", t, "\n", sep="")
#' 	print(Pmat)
#' 	}
#' 
expokit_dgpadm_Qmat <- function(Qmat=NULL, t=2.1, transpose_needed=TRUE)
	{
	defaults = '
	Qmat=NULL
	t = 2.1
	transpose_needed=TRUE
	'
	
	# Check if Qmat is blank
	if (is.null(Qmat))
		{
		# Default Qmat
		cat("\nWARNING: expokit_dgpadm_Qmat() was provided a Qmat with value NULL.  Example Qmat provided instead\n")
		Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
		}
	# Check if t is blank
	if (is.null(t))
		{
		# Default Qmat
		stop("\nSTOP ERROR: expokit_dgpadm_Qmat() was provided a t (time or times list) with value NULL.  \n")
		}
	

	
	# FOR DGPADM
	# ideg = 
	# "(input) the degree of the diagonal Pade to be used.
	# a value of 6 is generally satisfactory."
	ideg = as.integer(6)

	# Order (numrows/numcols) of the matrix
	# "(input) order of H."
	m = as.integer(nrow(Qmat))

	# output matrix
	res = double(length=m*m)

	# Prepare input matrix
	matvec = Qmat
	if (transpose_needed == TRUE)
		{
		tmatvec = t(matvec)
		H = as.numeric(tmatvec)
		} else {
		H = as.numeric(matvec)
		}
	
	# "H(ldh,m)  : (input) argument matrix."
	# (ldh = numrows and m is numcols, or something)
	ldh = m
	
	# No tolin PADM
	# tol or t?  should be t
	# tol = as.double(1)
	
	# lwsp = length of wsp, the workspace
	# "wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1."
	lwsp = as.integer(4*m*m+ideg+1)
	wsp = double(length=lwsp)
	
	# "ipiv(m)   : (workspace)"
	ipiv = integer(length=m)
	
	#  "iexph     : (output) number such that wsp(iexph) points to exp(tH)"
	#  "            i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)"
	iexph = as.integer(0)
	
	# "ns        : (output) number of scaling-squaring used."
	ns = as.integer(0)
	
	# "iflag     : (output) exit flag."
	# "	*                      0 - no problem"
	# "	*                     <0 - problem"
	iflag = as.integer(0)
	
	# Run the function:
	res <- .C("wrapdgpadm_", as.integer(ideg), as.integer(m), as.double(t), as.double(H), as.integer(ldh), as.double(wsp), as.integer(lwsp), as.integer(ipiv), as.integer(iexph), as.integer(ns), as.integer(iflag))
	
	output = res[[6]]
	output_Pmat_is = seq(res[[9]], res[[9]]+m*m-1, by=1)
	output_Pmat = output[output_Pmat_is]
	output_Pmat = matrix(output_Pmat, nrow=m, byrow=TRUE)
	#print(output_Pmat)
	
	return(output_Pmat)
	}







#' EXPOKIT dmexpv matrix exponentiation on Q matrix
#'
#' This function converts a matrix to COO format and exponentiates
#' it via the EXPOKIT dmexpv function (designed for sparse matrices)
#' and wrapper functions \code{wrapalldmexpv_} around dmexpv.
#'
#' From EXPOKIT:\cr
#' \code{*     The method used is based on Krylov subspace projection}\cr
#' \code{*     techniques and the matrix under consideration interacts only}\cr
#' \code{*     via the external routine 'matvec' performing the matrix-vector} \cr
#' \code{*     product (matrix-free method).}\cr
#' \code{*}\cr
#' \code{*     This is a customised version for Markov Chains. This means that a}\cr
#' \code{*     check is done within this code to ensure that the resulting vector} \cr
#' \code{*     w is a probability vector, i.e., w must have all its components }\cr
#' \code{*     in [0,1], with sum equal to 1. This check is done at some expense}\cr
#' \code{*     and the user may try DGEXPV which is cheaper since it ignores }\cr
#' \code{*     probability constraints.}\cr
#'
#' COO (coordinated list) format is a compressed format that is\cr
#' required for EXPOKIT's sparse-matrix functions (like dmexpv and\cr
#' unlike EXPOKIT's padm-related functions.\cr
#'
#' COO (coordinated list) format is described here:\cr
#'
#' \url{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}\cr
#'
#' If \code{Qmat} is NULL (default), a default matrix is input.\cr
#'
#' @param Qmat an input Q transition matrix
#' @param t one or more time values to exponentiate by
#' @param inputprobs_for_fast If NULL (default), the full probability matrix (Pmat) is returned. However, the full
#' speed of EXPOKIT on sparse matrices will be exploited if inputprobs_for_fast=c(starting probabilities). In this case
#' these starting probabilities are input to \code{myDMEXPV} directly, as \code{v}, and \code{w}, the output probabilities,
#' are returned.
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' @param transform_to_coo_TF Should the matrix be tranposed to COO?  COO format is required
#' for EXPOKIT's sparse-matrix functions (like dmexpv and unlike the padm-related 
#' functions. Default TRUE; if FALSE, user must put a COO-formated matrix in \code{Qmat}. Supplying the
#' coo matrix is probably faster for repeated calculations on large matrices.
#' @param coo_n If a COO matrix is input, \code{coo_n} specified the order (# rows, equals # columns) of the matrix.
#' @param anorm \code{dmexpv} requires an initial guess at the norm of the matrix. Using the
#' @param check_for_0_rows If TRUE or a numeric value, the input Qmat is checked for all-zero rows, since these will crash the FORTRAN wrapalldmexpv function. A small nonzero value set to check_for_0_rows or the default (0.0000000000001) is input to  off-diagonal cells in the row (and the diagonal value is normalized), which should fix the problem.
#' R function \code{\link{norm}} might get slow with large matrices. If so, the user
#' can input a guess manually (\code{Lagrange} seems to just use 1 or 0, if I
#' recall correctly).
#' @return \code{tmpoutmat} the output matrix. \code{wrapalldmexpv_} produces
#' additional output relating to accuracy of the output matrix etc.; these can be
#' by a direct call of dmexpv.
#' @seealso \code{\link{mat2coo}}
#' @seealso \code{\link{expokit_dmexpv_wrapper}}
#' @seealso \code{\link{expokit_wrapalldmexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples	# Example:
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' 
#' # Exponentiate each with EXPOKIT's dmexpv (should be fast for large sparse matrices)
#' for (t in tvals)
#' 	{
#' 	Pmat = expokit_dmexpv_Qmat(Qmat=Qmat, t=t, transpose_needed=TRUE)
#' 	cat("\n\nTime=", t, "\n", sep="")
#' 	print(Pmat)
#' 	}
#' 
expokit_dmexpv_Qmat <- function(Qmat=NULL, t=2.1, inputprobs_for_fast=NULL, transpose_needed=TRUE, transform_to_coo_TF=TRUE, coo_n=NULL, anorm=NULL, check_for_0_rows=TRUE)
	{
	defaults = '
	Qmat=NULL
	t = 2.1
	inputprobs_for_fast=NULL
	transpose_needed=TRUE
	transform_to_coo_TF=TRUE
	coo_n=NULL
	anorm=NULL
	check_for_0_rows=FALSE
	check_for_0_rows=1e-15
	'
	
	matvec = Qmat
	
	# Check if Qmat is blank
	if (is.null(matvec))
		{
		# Default Qmat
		cat("\nWARNING: expokit_dmexpv_Qmat() was provided a Qmat with value NULL.  Example Qmat provided instead\n")
		matvec = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
		}
	# Check if t is blank
	if (is.null(t))
		{
		# Default Qmat
		stop("\nSTOP ERROR: expokit_dmexpv_Qmat() was provided a t (time or times list) with value NULL.  \n")
		}


	# Zero rows will crash the FORTRAN wrapalldmexpv function, and
	# therefore crash R.  This is annoying.
	if (is.null(inputprobs_for_fast))
		{ # Return the full Pmat (slow)
		#######################################################
		# Return the Pmat
		#######################################################
		
		# Zero rows will crash the FORTRAN wrapalldmexpv function, and
		# therefore crash R.  This is annoying.
		if (check_for_0_rows != FALSE)
			{
			# If not false, check_for_0_rows is either TRUE or numeric

			# Get T/F for rows with all zeros
			rows_w_all_zeros_TF = findrows_w_all_zeros(matvec)
			
			# If all FALSE, do nothing
			if (all(rows_w_all_zeros_TF == FALSE))
				{
				# Do nothing
				pass = 1
				} else {
				# indices of TRUE
				rows_allzero_indices = seq(1, length(rows_w_all_zeros_TF), 1)[rows_w_all_zeros_TF]

				# Here, you have to input a small value for each zero
				if (is.numeric(check_for_0_rows))
					{
					check_for_0_rows = check_for_0_rows
					} else {
					# 1e-15 appears to be the precision limit of the FORTRAN code
					check_for_0_rows = 1e-15
					}
				# Input the given value into all zeros
				newrowvals = rep(check_for_0_rows, ncol(matvec))
				matvec[rows_allzero_indices, ] = newrowvals
				diagonal_val = -1 * sum(newrowvals[-1])
				matvec[rows_allzero_indices, rows_allzero_indices] = diagonal_val
				
				cat("\nWARNING: ", sum(rows_w_all_zeros_TF), " rows of the Q matrix Qmat had all zeros. This will crash .Call('wrapalldmexpv_', ...)\nand therefore expokit_wrapalldmexpv_tvals() run with the inputprobs_for_fast=NULL option (producing a full Pmat),\nand therefore R.  Replacement value for 0:  check_for_0_rows=", check_for_0_rows, ".\n", sep="")
				}
			}
		}
		
	
	# Count the number of NON-zeros (nz)
	# and input the matrix size
	if (transform_to_coo_TF == TRUE)
		{
		# COO format
		# http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29

		# number of non-zeros
		nz  = sum(matvec != 0)
		
		# These vectors have that length
		ia  = integer(length=nz)
		ja  = integer(length=nz)
		a   = double(length=nz)	
		n=nrow(matvec)
		
		} else {
		n = coo_n
		# (And make a regular matrix from COO)

		# number of non-zeros
		# Assumes that coo-formatted matrix columns are
		# ia, ja, a
		nz  = sum(matvec[,"a"] != 0)
		}

	# ideg = degree of polynomial, 6 is usually considered sufficient
	ideg = as.integer(6)
	#n=nrow(Qmat)	# don't do here, possibly coo-formatted
	m=n-1
	# t=as.numeric(2.1)
	
	# v should have as many elements as n; first element = 1 (?)
	if (is.null(inputprobs_for_fast))
		{
		# Input space-fillers, these get written over by wrapall
		v=double(n)
		v[1] = 1
		# Input the input probabilities, these get used directly by myDMEXPV/myDGEXPV
		} else {
		v = double(n)
		v = inputprobs_for_fast
		}
	
	# w is the same length
	w = double(length=n)
	tol=as.numeric(0.01)
	
	# lwsp = length of wsp
	# wsp = workspace to hold various variables, cells of the matrix, etc.
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 4*(m+2)^2+ideg+1)
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 5*(m+2)^2+ideg+1)
	lwsp = as.integer(n*(m+2)+5*(m+2)^2+ideg+1)
	
	#lwsp = 100
	wsp = double(length=lwsp)
	
	# length of iwsp
	liwsp = max(m+2, 7)
	iwsp = integer(length=liwsp)
	
	#matvec = matrix(data=Q, nrow=n, byrow=TRUE)
	matvec = matvec
	
	# Don't transform if already coo
	if ((transform_to_coo_TF == TRUE) && (transpose_needed == TRUE))
		{
		tmatvec = t(matvec)
		} else {
		tmatvec = matvec
		}
	#rowSums(tmatvec)
	#colSums(tmatvec)
	
	# This might (?) get slow with large matrices -- doesn't seem to
	if ((exists("anorm") == FALSE) || is.null(anorm))
		{
		# Use the 1-norm or one-norm
		if (transform_to_coo_TF==FALSE && transpose_needed==FALSE)
			{
			tmpQmat1 = coo2mat(matvec, n=coo_n)
			tmpQmat2 = t(tmpQmat1)
			anorm = as.numeric(norm(tmpQmat2, type="O"))
			} else {
			anorm = as.numeric(norm(matvec, type="O"))
			}
		}
	
	# The itrace flag, if set to 1, results in dmexpv printing some details of
	# the function's run to screen.
	itrace = 0
	iflag = 0	
	
	# Make the input COO matrix
	# COO = coordinate list format, useful for sparse matrices with lots of zeros:
	# http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29
	# ia = rownum in the matrix
	# ja = colnum in the matrix
	# a  = value of that cell
	
	if (transform_to_coo_TF == TRUE)
		{		
		tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmatvec)
		} else {
		tmpmat_in_REXPOKIT_coo_fmt = matvec
		}
	# Either way, store the rows/columns in the input variables for FORTRAN
	ia = tmpmat_in_REXPOKIT_coo_fmt[,"ia"]
	ja = tmpmat_in_REXPOKIT_coo_fmt[,"ja"]
	a = tmpmat_in_REXPOKIT_coo_fmt[,"a"]

	# Run the wrapper function	
	if (is.null(inputprobs_for_fast))
		{
		######################################
		# Return the full Pmat (slow)
		######################################
		
		# Create the space for res (the returned Pmat)
		res = double(length=n*n)
		
		res2 <- .C("wrapalldmexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res))
	
		# wrapalldmexpv_ returns all kinds of stuff, list item 18 is the P matrix
		# However, this may be an inefficient use of the dmexpv sparse matrix capabilities (Hansen)
		# Try mydmexpv_ to just get the ancestral probabilities (w, res2[[5]])
		output_Pmat = matrix(res2[[18]], nrow=n, byrow=TRUE)
		
		return(output_Pmat)
		} else {
		######################################
		# Instead of returning the full Pmat (slow), just return the output probabilities (fast)
		######################################
		
		# Be sure to input the input probabilities
		v = inputprobs_for_fast
		
		res2 <- .C("mydmexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz))
		
		# w, list item #5, contains the output probabilities
		w_output_probs = res2[[5]]
		
		return(w_output_probs)
		}
	}
























#' EXPOKIT dmexpv wrapper function
#'
#' This function wraps the .C call to EXPOKIT for the dmexpv function.  Only the output probability
#' matrix is returned.
#'
#' @param n number of rows in Q matrix
#' @param m n-1
#' @param t the value to exponentiate the rate matrix by (often e.g. a time value)
#' @param v variable to store some results in; should have n elements (and perhaps start with 1)
#' @param w same length as v
#' @param tol tolerance for approximations; usually set to 0.01
#' @param anorm the norm of the Q matrix
#' @param lwsp length of workspace (wsp); for dmexpv, lwsp=n*(m+2)+5*(m+2)^2+ideg+1
#' @param wsp workspace to store some results in; should be a double with lwsp elements
#' @param liwsp length of integer workspace; for dmexpv, liwsp=m+2
#' @param iwsp integer workspace
#' @param itrace option, set to 0
#' @param iflag option, set to 0
#' @param ia i indices of Qmat nonzero values
#' @param ja j indices of Qmat nonzero values
#' @param a nonzero values of Qmat (ia, ja, a are columns of a COO-formatted Q matrix)
#' @param nz number of non-zeros in Qmat
#' @param res space for output probability matrix (n x n)
#'
#' EXPOKIT needs the input matrix to be transposed compared to normal.
#' COO format is required for EXPOKIT.
#' @return \code{tmpoutmat} the output matrix for the (first) input t-value
#' @seealso \code{\link{expokit_dmexpv_Qmat}}
#' @seealso \code{\link{expokit_wrapalldmexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' 
#' # DMEXPV and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
#' # DMEXPV is specifically designed for Markov chains and so may be slower, but more accurate.
#' 
#' # DMEXPV, single t-value
#' expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
#' expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=2)
#' 
#' # This function runs a for-loop itself (sadly, we could not get mapply() to work
#' # on a function that calls dmexpv/dgexpv), returning a list of probability matrices.
#' 
#' # DMEXPV functions
#' list_of_P_matrices_dmexpv = expokit_wrapalldmexpv_tvals(Qmat=Qmat, 
#' tvals=tvals, transpose_needed=TRUE)
#' list_of_P_matrices_dmexpv
#' 
expokit_dmexpv_wrapper <- function(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)
	{
	res2 = NULL
	
	res2 <- .C("wrapalldmexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res))
	
	output_Pmat = matrix(res2[[18]], nrow=n, byrow=TRUE)
	
	return(output_Pmat)
	}




#' EXPOKIT dmexpv wrapper function, return just output probs
#'
#' This function wraps the .C call to EXPOKIT for the dmexpv function.  Only the output probabilities
#' not the Pmat probability matrix, are returned.
#'
#' @param n number of rows in Q matrix
#' @param m n-1
#' @param t the value to exponentiate the rate matrix by (often e.g. a time value)
#' @param v variable to store some results in; should have n elements (and perhaps start with 1)
#' @param w same length as v
#' @param tol tolerance for approximations; usually set to 0.01
#' @param anorm the norm of the Q matrix
#' @param lwsp length of workspace (wsp); for dmexpv, lwsp=n*(m+2)+5*(m+2)^2+ideg+1
#' @param wsp workspace to store some results in; should be a double with lwsp elements
#' @param liwsp length of integer workspace; for dmexpv, liwsp=m+2
#' @param iwsp integer workspace
#' @param itrace option, set to 0
#' @param iflag option, set to 0
#' @param ia i indices of Qmat nonzero values
#' @param ja j indices of Qmat nonzero values
#' @param a nonzero values of Qmat (ia, ja, a are columns of a COO-formatted Q matrix)
#' @param nz number of non-zeros in Qmat
#'
#' EXPOKIT needs the input matrix to be transposed compared to normal.
#' COO format is required for EXPOKIT.
#' @return \code{w_output_probs} the output probabilities (= \code{myDMEXPV} variable \code{w}, or the fifth output
#' in the output from .Call("mydmexpv_", ...), given the (first) input t-value.
#' @seealso \code{\link{expokit_dmexpv_Qmat}}
#' @seealso \code{\link{expokit_wrapalldmexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' 
#' # DMEXPV and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
#' # DMEXPV is specifically designed for Markov chains and so may be slower, but more accurate.
#' 
#' # DMEXPV, single t-value
#' expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
#' expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=2)
#' 
#' # This function runs a for-loop itself (sadly, we could not get mapply() to work
#' # on a function that calls dmexpv/dgexpv), returning a list of probability matrices.
#' 
#' # DMEXPV functions
#' list_of_P_matrices_dmexpv = expokit_wrapalldmexpv_tvals(Qmat=Qmat, 
#' tvals=tvals, transpose_needed=TRUE)
#' list_of_P_matrices_dmexpv
#' 
expokit_mydmexpv_wrapper <- function(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz)
	{
	res2 = NULL
	
	# This must be mydmexpv_, not myDMEXPV_ !!!!
	
	res2 <- .C("mydmexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz))
	
	w_output_probs = matrix(res2[[5]], ncol=n, byrow=TRUE)
	
	return(w_output_probs)
	}
	


#' EXPOKIT dgexpv wrapper function, return just output probs
#'
#' This function wraps the .C call to EXPOKIT for the dgexpv function.  Only the output probabilities
#' not the Pmat probability matrix, are returned.
#'
#' @param n number of rows in Q matrix
#' @param m n-1
#' @param t the value to exponentiate the rate matrix by (often e.g. a time value)
#' @param v variable to store some results in; should have n elements (and perhaps start with 1)
#' @param w same length as v
#' @param tol tolerance for approximations; usually set to 0.01
#' @param anorm the norm of the Q matrix
#' @param lwsp length of workspace (wsp); for dgexpv, lwsp=n*(m+2)+5*(m+2)^2+ideg+1
#' @param wsp workspace to store some results in; should be a double with lwsp elements
#' @param liwsp length of integer workspace; for dgexpv, liwsp=m+2
#' @param iwsp integer workspace
#' @param itrace option, set to 0
#' @param iflag option, set to 0
#' @param ia i indices of Qmat nonzero values
#' @param ja j indices of Qmat nonzero values
#' @param a nonzero values of Qmat (ia, ja, a are columns of a COO-formatted Q matrix)
#' @param nz number of non-zeros in Qmat
#'
#' EXPOKIT needs the input matrix to be transposed compared to normal.
#' COO format is required for EXPOKIT.
#' @return \code{w_output_probs} the output probabilities (= \code{myDGEXPV} variable \code{w}, or the fifth output
#' in the output from .Call("mydgexpv_", ...), given the (first) input t-value.
#' @seealso \code{\link{expokit_dgexpv_Qmat}}
#' @seealso \code{\link{expokit_wrapalldgexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' 
#' # dgexpv and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
#' # dgexpv is specifically designed for Markov chains and so may be slower, but more accurate.
#' 
#' # dgexpv, single t-value
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=2)
#' 
#' # This function runs a for-loop itself (sadly, we could not get mapply() to work
#' # on a function that calls dgexpv/dgexpv), returning a list of probability matrices.
#' 
#' # dgexpv functions
#' list_of_P_matrices_dgexpv = expokit_wrapalldgexpv_tvals(Qmat=Qmat, 
#' tvals=tvals, transpose_needed=TRUE)
#' list_of_P_matrices_dgexpv
#' 
expokit_mydgexpv_wrapper <- function(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz)
	{
	res2 = NULL
	
	# This must be mydgexpv_, not mydgexpv_ !!!!
	
	res2 <- .C("mydgexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz))
	
	w_output_probs = matrix(res2[[5]], ncol=n, byrow=TRUE)
	
	return(w_output_probs)
	}
	


#' Check if a row is all zeros
#'
#' Q matrices with all-zero rows will crash .Call(wrapalldmexpv_, ...) and .Call(wrapalldgexpv_, ...),
#' and therefore will crash expokit_wrapalldmexpv_tvals() and expokit_wrapalldgexpv_tvals() when 
#' these are set (the default) to return the full P matrix.  These functions work fine with
#' zero rows if \code{inputprobs_for_fast} is supplied, meaning that only the output probabilities
#' of each state are returned.
#'
#' @param tmprow A row of a Q transition matrix
#' @return \code{TRUE} if tmprow is all zeros, FALSE if not.
#' @seealso \code{\link{expokit_wrapalldmexpv_tvals}}
#' @seealso \code{\link{expokit_wrapalldgexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' 
#' # DMEXPV and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
#' # DMEXPV is specifically designed for Markov chains and so may be slower, but more accurate.
#' 
#' # DGEXPV, single t-value
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=2)
#' 
#' # This function runs the for-loop itself (sadly, we could not get mapply() to work
#' # on a function that calls dmexpv/dgexpv), returning a list of probability matrices.
#' 
#' # DGEXPV functions
#' list_of_P_matrices_dgexpv = expokit_wrapalldgexpv_tvals(Qmat=Qmat, 
#' tvals=tvals, transpose_needed=TRUE)
#' list_of_P_matrices_dgexpv
row_allzero_TF <- function(tmprow)
	{
	return(all(tmprow == 0))
	}

#' Check if a Q matrix has rows with all zeros
#'
#' Q matrices with all-zero rows will crash .Call(wrapalldmexpv_, ...) and .Call(wrapalldgexpv_, ...),
#' and therefore will crash expokit_wrapalldmexpv_tvals() and expokit_wrapalldgexpv_tvals() when 
#' these are set (the default) to return the full P matrix.  These functions work fine with
#' zero rows if \code{inputprobs_for_fast} is supplied, meaning that only the output probabilities
#' of each state are returned.
#'
#' @param matvec Q transition matrix
#' @return A list of TRUE/FALSE, as long as the number of rows. \code{TRUE}=the is all zeros, \code{FALSE}=the row has nonzero values.
#' @seealso \code{\link{expokit_wrapalldmexpv_tvals}}
#' @seealso \code{\link{expokit_wrapalldgexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' 
#' # DMEXPV and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
#' # DMEXPV is specifically designed for Markov chains and so may be slower, but more accurate.
#' 
#' # DGEXPV, single t-value
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=2)
#' 
#' # This function runs the for-loop itself (sadly, we could not get mapply() to work
#' # on a function that calls dmexpv/dgexpv), returning a list of probability matrices.
#' 
#' # DGEXPV functions
#' list_of_P_matrices_dgexpv = expokit_wrapalldgexpv_tvals(Qmat=Qmat, 
#' tvals=tvals, transpose_needed=TRUE)
#' list_of_P_matrices_dgexpv
findrows_w_all_zeros <- function(matvec)
	{
	return(apply(X=matvec, MARGIN=1, FUN=row_allzero_TF))
	}



#' Run EXPOKIT's dmexpv on one or more t-values
#'
#' The function runs EXPOKIT's \code{dmexpv} function on a Q matrix and \emph{one or more} time values.  If \code{Qmat} is NULL (default), a default matrix is input.
#'
#' @param Qmat an input Q transition matrix
#' @param tvals one or more time values to exponentiate by (doesn't have to literally be a time value, obviously)
#' @param inputprobs_for_fast If NULL (default), the full probability matrix (Pmat) is returned. However, the full
#' speed of EXPOKIT on sparse matrices will be exploited if inputprobs_for_fast=c(starting probabilities). In this case
#' these starting probabilities are input to \code{myDMEXPV} directly, as \code{v}, and \code{w}, the output probabilities,
#' are returned.
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' @param transform_to_coo_TF Should the matrix be tranposed to COO?  COO format is required
#' for EXPOKIT's sparse-matrix functions (like dmexpv and unlike the padm-related 
#' functions. Default TRUE; if FALSE, user must put a COO-formated matrix in \code{Qmat}. Supplying the
#' coo matrix is probably faster for repeated calculations on large matrices.
#' @param coo_n If a COO matrix is input, \code{coo_n} specified the order (# rows, equals # columns) of the matrix.
#' @param force_list_if_1_tval Default FALSE, but set to TRUE if you want a single matrix to be returned
#' inside a list
#' @param check_for_0_rows If TRUE or a numeric value, the input Qmat is checked for all-zero rows, since these will crash the FORTRAN wrapalldmexpv function. A small nonzero value set to check_for_0_rows or the default (0.0000000000001) is input to  off-diagonal cells in the row (and the diagonal value is normalized), which should fix the problem.
#' @return \code{tmpoutmat} the output matrix, if 1 t-value is input; \code{list_of_matrices_output},
#' if more than 1 t-value is input; to get a single output matrix in a list, set \code{force_list_if_1_tval=TRUE}
#' @seealso \code{\link{expokit_dmexpv_wrapper}}
#' @seealso \code{\link{expokit_dmexpv_Qmat}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' 
#' # DMEXPV and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
#' # DMEXPV is specifically designed for Markov chains and so may be slower, but more accurate.
#' 
#' # DGEXPV, single t-value
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=2)
#' 
#' # This function runs the for-loop itself (sadly, we could not get mapply() to work
#' # on a function that calls dmexpv/dgexpv), returning a list of probability matrices.
#' 
#' # DGEXPV functions
#' list_of_P_matrices_dgexpv = expokit_wrapalldgexpv_tvals(Qmat=Qmat, 
#' tvals=tvals, transpose_needed=TRUE)
#' list_of_P_matrices_dgexpv
#' 
expokit_wrapalldmexpv_tvals <- function(Qmat=NULL, tvals=c(2.1), inputprobs_for_fast=NULL, transpose_needed=TRUE, transform_to_coo_TF=TRUE, coo_n=NULL, force_list_if_1_tval=FALSE, check_for_0_rows=TRUE)
	{
	defaults = '
	Qmat=NULL
	t = 2.1
	inputprobs_for_fast=NULL
	transpose_needed=TRUE
	COO_needed=TRUE
	transform_to_coo_TF=TRUE
	coo_n=NULL
	force_list_if_1_tval=FALSE
	check_for_0_rows=FALSE
	check_for_0_rows=1e-15
	'
	
	# Check if Qmat is blank
	if (is.null(Qmat))
		{
		# Default Qmat
		warning("You supplied no matrix, so a default matrix is being used. Obviously you can't use this for anything real. YOU HAVE BEEN WARNED!!")
		
		Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
		}
	
	if (is.null(tvals))
		{
		warning("You supplied no time values, so default time values are being used. Obviously you can't use this for anything real. YOU HAVE BEEN WARNED!!")
		tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
		}

	matvec = Qmat
	
	# Zero rows will crash the FORTRAN wrapalldmexpv function, and
	# therefore crash R.  This is annoying.
	if (is.null(inputprobs_for_fast))
		{ # Return the full Pmat (slow)
		#######################################################
		# Return the Pmat
		#######################################################
		
		# Zero rows will crash the FORTRAN wrapalldmexpv function, and
		# therefore crash R.  This is annoying.
		if (check_for_0_rows != FALSE)
			{
			# If not false, check_for_0_rows is either TRUE or numeric

			# Get T/F for rows with all zeros
			rows_w_all_zeros_TF = findrows_w_all_zeros(matvec)
			
			# If all FALSE, do nothing
			if (all(rows_w_all_zeros_TF == FALSE))
				{
				# Do nothing
				pass = 1
				} else {
				# indices of TRUE
				rows_allzero_indices = seq(1, length(rows_w_all_zeros_TF), 1)[rows_w_all_zeros_TF]

				# Here, you have to input a small value for each zero
				if (is.numeric(check_for_0_rows))
					{
					check_for_0_rows = check_for_0_rows
					} else {
					# 1e-15 appears to be the precision limit of the FORTRAN code
					check_for_0_rows = 1e-15
					}
				# Input the given value into all zeros
				newrowvals = rep(check_for_0_rows, ncol(matvec))
				matvec[rows_allzero_indices, ] = newrowvals
				diagonal_val = -1 * sum(newrowvals[-1])
				matvec[rows_allzero_indices, rows_allzero_indices] = diagonal_val
				
				cat("\nWARNING: ", sum(rows_w_all_zeros_TF), " rows of the Q matrix Qmat had all zeros. This will crash .Call('wrapalldmexpv_', ...)\nand therefore expokit_wrapalldmexpv_tvals() run with the inputprobs_for_fast=NULL option (producing a full Pmat),\nand therefore R.  Replacement value for 0:  check_for_0_rows=", check_for_0_rows, ".\n", sep="")
				}
			}
		}	

		
	# Count the number of NON-zeros (nz)
	# and input the matrix size
	if (transform_to_coo_TF == TRUE)
		{
		# COO format
		# http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29

		# number of non-zeros
		nz  = sum(matvec != 0)
		
		# These vectors have that length
		ia  = integer(length=nz)
		ja  = integer(length=nz)
		a   = double(length=nz)	
		n=nrow(matvec)
		} else {
		n = coo_n
		# (And make a regular matrix from COO)

		# number of non-zeros
		# Assumes that coo-formatted matrix columns are
		# ia, ja, a
		nz  = sum(matvec[,"a"] != 0)
		}
	
	

	#######################################################
	ideg = as.integer(6)
	#######################################################
	#n=nrow(Qmat)	# don't do this here, you might have a coo matrix
	m=n-1
	# t=as.numeric(2.1)
	
	# v should have as many elements as n; first element = 1 (?)
	if (is.null(inputprobs_for_fast))
		{
		# Input space-fillers, these get written over by wrapall
		v=double(n)
		v[1] = 1
		# Input the input probabilities, these get used directly by myDMEXPV/myDGEXPV
		} else {
		v = double(n)
		v = inputprobs_for_fast
		}
	
	# w is the same length
	w = double(length=n)
	tol=as.numeric(0.01)
	
	# length of wsp
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 4*(m+2)^2+ideg+1)
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 5*(m+2)^2+ideg+1)
	lwsp = as.integer(n*(m+2)+5*(m+2)^2+ideg+1)
	
	#lwsp = 100
	wsp = double(length=lwsp)
	
	# length of iwsp
	liwsp = max(m+2, 7)
	iwsp = integer(length=liwsp)
	

	if ((transform_to_coo_TF == TRUE) && (transpose_needed == TRUE))
		{
		tmatvec = t(matvec)
		#rowSums(tmatvec)
		#colSums(tmatvec)
		}


	
	# type="O" is being used here, this is supposed to be the
	# default for norm(), although it throws an error if not
	# specified
	# 
	# From the help:
	# type - character string, specifying the type of matrix norm to be
	# computed. A character indicating the type of norm desired. 
	# 	"O", "o" or "1"
	# 		specifies the one norm, (maximum absolute column sum);
	if ((exists("anorm") == FALSE) || is.null(anorm))
		{
		# Use the 1-norm or one-norm
		if (transform_to_coo_TF==FALSE && transpose_needed==FALSE)
			{
			tmpQmat1 = coo2mat(matvec, n=coo_n)
			tmpQmat2 = t(tmpQmat1)
			anorm = as.numeric(norm(tmpQmat2, type="O"))
			} else {
			anorm = as.numeric(norm(matvec, type="O"))
			}
		}	
	
	itrace = 0
	iflag = 0	
	
	
	#a = as.numeric(tmatvec)
	#a = as.numeric(matvec)
	if (transform_to_coo_TF == TRUE)
		{		
		tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmatvec)
		} else {
		tmpmat_in_REXPOKIT_coo_fmt = matvec
		}
	# Either way, store the rows/columns in the input variables for FORTRAN
	ia = tmpmat_in_REXPOKIT_coo_fmt[,"ia"]
	ja = tmpmat_in_REXPOKIT_coo_fmt[,"ja"]
	a = tmpmat_in_REXPOKIT_coo_fmt[,"a"]

	num_tvals = length(tvals)


	# Run the wrapper function	
	if (is.null(inputprobs_for_fast))
		{ # Return the full Pmat (slow)
		#######################################################
		# Return the Pmat
		#######################################################
		
		# Create the space for res (the returned Pmat)
		res = double(length=n*n)

	
		# If there is more than 1 t-value, or if the user desires a list even for a single
		# t-value, return a list
		if ((num_tvals > 1) || (force_list_if_1_tval==TRUE))
			{

			# Loop through the list of tvals, get the prob. matrix for each
			# sadly, mapply() etc. crash when tried on expokit_dmexpv_wrapper
			
			# Set up empty matrix
			NA_matrix = matrix(NA, nrow=n, ncol=n)
			
			# Set up list of empty matrices
			list_of_matrices_output = replicate(NA_matrix, n=num_tvals, simplify=FALSE)
			
			for (i in 1:num_tvals)
				{
				t = tvals[i]
				list_of_matrices_output[[i]] = expokit_dmexpv_wrapper(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)
	
				} # end forloop
			
			return(list_of_matrices_output)
			
			} else {
			# If there is only 1 t value, just return 1 matrix
			#res2 <- .C("wrapalldmexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res))
			
			#tmpoutmat = matrix(res2[[18]], nrow=n, byrow=TRUE)
	
			t=tvals
			output_Pmat = expokit_dmexpv_wrapper(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)		
			
			#print(tmpoutmat)
			return(output_Pmat)
			}
		} else {
		#######################################################
		# Return the output probabilities
		#######################################################
				
		# Be sure to input the input probabilities
		v = inputprobs_for_fast

		# If there is more than 1 t-value, or if the user desires a list even for a single
		# t-value, return a list
		if ((num_tvals > 1) || (force_list_if_1_tval==TRUE))
			{
			# Loop through the list of tvals, get the prob. matrix for each
			# sadly, mapply() etc. crash when tried on expokit_dmexpv_wrapper
			
			# Set up empty matrix
			list_of_outprobs_output = matrix(NA, nrow=num_tvals, ncol=n)

			for (i in 1:num_tvals)
				{
				t = tvals[i]
			
				# This must be mydmexpv_, not myDMEXPV_ !!!!

				res2 <- .C("mydmexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz))
		
				# w, list item #5, contains the output probabilities
				w_output_probs = res2[[5]]

				list_of_outprobs_output[i,] = w_output_probs
	
				}
			} else {
			
			# If there is only 1 t value, just return 1 matrix
			#res2 <- .C("wrapalldmexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res))
			
			#tmpoutmat = matrix(res2[[18]], nrow=n, byrow=TRUE)
	
			t=tvals

			# This must be mydmexpv_, not myDMEXPV_ !!!!

			res2 <- .C("mydmexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz))
			
			# w, list item #5, contains the output probabilities
			w_output_probs = res2[[5]]

			list_of_outprobs_output[1,] = w_output_probs
	
			return(w_output_probs)
			}

		return(list_of_matrices_output)
		}
	}



#' Convert matrix to COO format using SparseM function
#'
#' Converts a matrix to COO format using the SparseM function, presumably this
#' is faster than using a for-loop.\cr
#'
#' \code{EXPOKIT}'s \code{dmexp}-type functions deal with sparse matrices.
#' These have a lot of zeros, and thus can be compressed
#' into COO (coordinated list) format, which is described here:\cr
#'
#' \url{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}\cr
#'
#' In \code{EXPOKIT} and its wrapper functions, a COO-formated matrix is input as
#' 3 vectors (first two integer, the third double):\cr
#'
#' ia = row number\cr
#' ja = column number\cr
#' a = value of that cell in the matrix (skipping 0 cells)\cr
#' 
#' @param tmpmat A square matrix
#' @return tmpmat_in_REXPOKIT_coo_fmt A \code{cbind} of \code{ia}, \code{ja}, and \code{a} 
#' @seealso \code{\link{mat2coo_forloop}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples # Example use:
#' @examples
#' # Make a Q matrix
#' tmpmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Convert to coo format
#' tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmpmat)
#' tmpmat_in_REXPOKIT_coo_fmt
#' 
mat2coo <- function(tmpmat)
	{
	defaults = '
	tmpmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
	'
	
	numrows = nrow(tmpmat)
	numcols = ncol(tmpmat)

	if (numrows != numcols)
		{
		stop("ERROR! mat2coo(tmpmat) says that in tmpmat, nrows != ncols")
		return(NA)
		}

	numcells = numrows ^2
	
	# require(sfsmisc)
	# xy.grid
	
	x = 1:numrows
	y = 1:numrows

	# Seems to be slow when numrow > 1000
	# as.vector(tmpmat) appends col1vals, then col2vals, etc., so ji = xy
	# cells_ij = expand.grid(x, y)
	# tmpa = as.vector(tmpmat)
	# 
	# # Remove 0s
	# TF = tmpa != 0
	# ia = cells_ij[,1][TF]	
	# ja = cells_ij[,2][TF]	
	# a = tmpa[TF]

	require(SparseM)	# required for the as.matrix.coo function
	
	# This produces a matrix in coo format
	# (this is an S4 object)
	tmpmat_in_SparseMcoo_fmt = as.matrix.coo(tmpmat)
	tmpmat_in_REXPOKIT_coo_fmt = SparseM_coo_to_REXPOKIT_coo(tmpmat_in_SparseMcoo_fmt)
	
	return(tmpmat_in_REXPOKIT_coo_fmt)
	}



#' Convert a SparseM COO matrix to a plain matrix
#'
#' Converts a SparseM COO-formatted matrix (an S4 object) to a plain matrix, with \cr
#' column #1 = ia = i index\cr
#' column #2 = ja = j index\cr
#' column #3 = a = nonzero values of the matrix\cr
#'
#' Background: COO (coordinated list) format, is described here:\cr
#'
#' \url{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}\cr
#'
#' In \code{EXPOKIT} and its wrapper functions, a COO-formated matrix is input as
#' 3 vectors (first two integer, the third double):\cr
#'
#' ia = row number\cr
#' ja = column number\cr
#' a = value of that cell in the matrix (skipping 0 cells)\cr
#' 
#' @param tmpmat_in_SparseMcoo_fmt A square matrix S4 object derived from SparseM's as.matrix.coo
#' @return tmpmat_in_REXPOKIT_coo_fmt A \code{cbind} of \code{ia}, \code{ja}, and \code{a} 
#' @seealso \code{\link{mat2coo_forloop}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples # Example use:
#' # Make a Q matrix
#' tmpmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Covert to SparseM coo format
#' tmpmat_in_SparseMcoo_fmt = as.matrix.coo(tmpmat)
#' 
#' # Convert to REXPOKIT coo format
#' tmpmat_in_REXPOKIT_coo_fmt = SparseM_coo_to_REXPOKIT_coo(tmpmat_in_SparseMcoo_fmt)
#' tmpmat_in_REXPOKIT_coo_fmt
#' 
SparseM_coo_to_REXPOKIT_coo <- function(tmpmat_in_SparseMcoo_fmt)
	{
	tmpcoo = tmpmat_in_SparseMcoo_fmt
	
	# We just need the 3 columns: i index, j index, and nonzero values
	tmpmat_in_REXPOKIT_coo_fmt = cbind(tmpcoo@ia, tmpcoo@ja, tmpcoo@ra)
	
	# Apply appropriate column names
	colnames(tmpmat_in_REXPOKIT_coo_fmt) = c("ia", "ja", "a")
	
	return(tmpmat_in_REXPOKIT_coo_fmt)
	}




#' Convert a COO-formated matrix to standard square format
#'
#' \code{EXPOKIT}'s \code{dmexp}-type functions deal with sparse matrices.
#' These have a lot of zeros, and thus can be compressed
#' into COO (coordinated list) format, which is described here:\cr
#'
#' \url{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}\cr
#'
#' In \code{EXPOKIT} and its wrapper functions, a COO-formated matrix is input as
#' 3 vectors (first two integer, the third double):\cr
#'
#' ia = row number\cr
#' ja = column number\cr
#' a = value of that cell in the matrix (skipping 0 cells)\cr
#'
#' This function takes a 3-column matrix or data.frame (basically \code{cbind(ia, ja, a)})
#' and the order of the matrix, \code{n} (n = the order of the matrix, i.e. number of
#' rows/columns) and converts back to standard square format.\cr
#'
#' @param coomat a 3-column matrix or data.frame (basically \code{cbind(ia, ja, a)})
#' @param n the order of the matrix
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' @return outmat
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples # Example use:
#' ia = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)
#' ja = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4)
#' a  = c(-1.218, 0.126, 0.168, 0.126, 0.504, -0.882, 0.504, 
#' 0.672, 0.336, 0.252, -1.050, 0.252, 0.378, 0.504, 0.378, -1.050)
#' coomat = cbind(ia, ja, a)
#' print(coomat)
#' n = 4
#' Qmat = coo2mat(coomat, n)
#' print(Qmat)
coo2mat <- function(coomat, n=max(max(coomat[,1]), max(coomat[,2])), transpose_needed=FALSE)
	{
	defaults='
	
	'
	
	# Make an empty matrix of 0s
	outmat = matrix(double(length=n*n), nrow=n)
	
	# go through each row of coomat
	ia = coomat[,1]
	ja = coomat[,2]
	a = coomat[,3]
	
	if (transpose_needed == FALSE)
		{
		for (k in 1:length(ia))
			{
			#cat(ia[k], ja[k], a[k], "\n")
			outmat[ia[k], ja[k]] = a[k]
			}
		} else {
		for (k in 1:length(ia))
			{
			#cat(ia[k], ja[k], a[k], "\n")
			outmat[ja[k], ia[k]] = a[k]
			}		
		}
	
	return(outmat)
	}



#' Convert matrix to COO format using nested for-loops
#'
#' Converts a matrix to COO format. This version of the function uses
#' for-loops, which is presumably less efficient than \code{\link{mat2coo}}.
#'
#' @param tmpmat A square matrix
#' @return tmpmat_in_REXPOKIT_coo_fmt A \code{cbind} of \code{ia}, \code{ja}, and \code{a} 
#' @seealso \code{\link{mat2coo}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples # Example use:
#' # Make a Q matrix
#' tmpmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Convert to REXPOKIT coo format
#' tmpmat_in_REXPOKIT_coo_fmt = mat2coo_forloop(tmpmat)
#' tmpmat_in_REXPOKIT_coo_fmt
#' 
mat2coo_forloop <- function(tmpmat)
	{
	# Number of non-zeros
	nz = sum(tmpmat != 0)
	
	# Blank columns for COO matrix
	ia = integer(length=nz)
	ja = integer(length=nz)
	a = integer(length=nz)
	
	count = 0
	for (i in 1:nrow(tmpmat))
		{
		for (j in 1:ncol(tmpmat))
			{
			if (tmpmat[i,j] != 0)
				{
				count = count+1
				ia[count] = i
				ja[count] = j
				a[count] = tmpmat[i,j]
				}
			}
		}
	tmpmat_in_REXPOKIT_coo_fmt = cbind(ia, ja, a)
	return(tmpmat_in_REXPOKIT_coo_fmt)
	}







#######################################################
# 
# NOTE: DGEXPV section.  Same code as dmexpv, but EXPOKIT's DGEXPV should be
# faster than DMEXPV, however DMEXPV runs an accuracy check appropriate for
# Markov chains, which is not done in DGEXPV.
#
#######################################################




#' EXPOKIT dgexpv matrix exponentiation on Q matrix
#'
#' This function converts a matrix to COO format and exponentiates
#' it via the EXPOKIT dgexpv function (designed for sparse matrices)
#' and wrapper functions \code{wrapalldgexpv_} around dgexpv.\cr
#'
#'
#' NOTE: DGEXPV vs. DMEXPV. According to the EXPOKIT documentation, DGEXPV should be
#' faster than DMEXPV, however DMEXPV runs an accuracy check appropriate for
#' Markov chains, which is not done in DGEXPV.\cr
#'
#' From EXPOKIT:\cr
#'
#' \code{*     The method used is based on Krylov subspace projection}\cr
#' \code{*     techniques and the matrix under consideration interacts only}\cr
#' \code{*     via the external routine 'matvec' performing the matrix-vector} \cr
#' \code{*     product (matrix-free method).}\cr
#' \code{*}\cr
#' \code{*     This [DMEXPV, not DGEXPV -- NJM] is a customised version for Markov Chains. This means that a}\cr
#' \code{*     check is done within this code to ensure that the resulting vector }\cr
#' \code{*     w is a probability vector, i.e., w must have all its components }\cr
#' \code{*     in [0,1], with sum equal to 1. This check is done at some expense}\cr
#' \code{*     and the user may try DGEXPV which is cheaper since it ignores }\cr
#' \code{*     probability constraints.}\cr
#'
#' I (NJM) have not noticed a difference between the outputs of these two functions, but it might
#' occur with large matrices.
#'
#' COO (coordinated list) format is a compressed format that is
#' required for EXPOKIT's sparse-matrix functions (like dgexpv and
#' unlike EXPOKIT's padm-related functions. COO format is described here:\cr
#'
#' \url{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}\cr
#'
#' If \code{Qmat} is NULL (default), a default matrix is input.\cr
#'
#' @param Qmat an input Q transition matrix
#' @param t a time value to exponentiate by
#' @param inputprobs_for_fast If NULL (default), the full probability matrix (Pmat) is returned. However, the full
#' speed of EXPOKIT on sparse matrices will be exploited if inputprobs_for_fast=c(starting probabilities). In this case
#' these starting probabilities are input to \code{myDMEXPV} directly, as \code{v}, and \code{w}, the output probabilities,
#' are returned.
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' @param transform_to_coo_TF Should the matrix be tranposed to COO?  COO format is required
#' for EXPOKIT's sparse-matrix functions (like dmexpv and unlike the padm-related 
#' functions. Default TRUE; if FALSE, user must put a COO-formated matrix in \code{Qmat}. Supplying the
#' coo matrix is probably faster for repeated calculations on large matrices.
#' @param coo_n If a COO matrix is input, \code{coo_n} specified the order (# rows, equals # columns) of the matrix.
#' @param anorm \code{dgexpv} requires an initial guess at the norm of the matrix. Using the
#' R function \code{\link{norm}} might get slow with large matrices. If so, the user
#' can input a guess manually (\code{Lagrange} seems to just use 1 or 0, if I
#' recall correctly).
#' @param check_for_0_rows If TRUE or a numeric value, the input Qmat is checked for all-zero rows, since these will crash the FORTRAN wrapalldmexpv function. A small nonzero value set to check_for_0_rows or the default (0.0000000000001) is input to  off-diagonal cells in the row (and the diagonal value is normalized), which should fix the problem.
#' @return \code{tmpoutmat} the output matrix. \code{wrapalldgexpv_} produces
#' additional output relating to accuracy of the output matrix etc.; these can be
#' by a direct call of dgexpv.
#' @seealso \code{\link{mat2coo}}
#' @seealso \code{\link{expokit_dgexpv_wrapper}}
#' @seealso \code{\link{expokit_wrapalldgexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples 	# Example:
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' 
#' # Exponentiate each with EXPOKIT's dgexpv (should be fast for large sparse matrices)
#' for (t in tvals)
#' 	{
#' 	Pmat = expokit_dgexpv_Qmat(Qmat=Qmat, t=t, transpose_needed=TRUE)
#' 	cat("\n\nTime=", t, "\n", sep="")
#' 	print(Pmat)
#' 	}
#'
#' # DMEXPV and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
#' # DMEXPV is specifically designed for Markov chains and so may be slower, but more accurate.
#' 
#' # DGEXPV, single t-value
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=2)
#' 
#' # This function runs the for-loop itself (sadly, we could not get mapply() to work
#' # on a function that calls dmexpv/dgexpv), returning a list of probability matrices.
#' 
#' # DGEXPV functions
#' list_of_P_matrices_dgexpv = expokit_wrapalldgexpv_tvals(Qmat=Qmat, 
#' tvals=tvals, transpose_needed=TRUE)
#' list_of_P_matrices_dgexpv
#' 
expokit_dgexpv_Qmat <- function(Qmat=NULL, t=2.1, inputprobs_for_fast=NULL, transpose_needed=TRUE, transform_to_coo_TF=TRUE, coo_n=NULL, anorm=NULL, check_for_0_rows=TRUE)
	{
	defaults = '
	Qmat=NULL
	t = 2.1
	inputprobs_for_fast=NULL
	transpose_needed=TRUE
	transform_to_coo_TF=TRUE
	coo_n=NULL
	anorm=NULL
	check_for_0_rows=FALSE
	check_for_0_rows=1e-15
	'
	
	matvec = Qmat
	
	# Check if Qmat is blank
	if (is.null(matvec))
		{
		# Default Qmat
		cat("\nWARNING: expokit_dgexpv_Qmat() was provided a Qmat with value NULL.  Example Qmat provided instead\n")
		matvec = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168,  0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
		}
	# Check if t is blank
	if (is.null(t))
		{
		# Default Qmat
		stop("\nSTOP ERROR: expokit_dgexpv_Qmat() was provided a t (time or times list) with value NULL.  \n")
		}


	# Zero rows will crash the FORTRAN wrapalldgexpv function, and
	# therefore crash R.  This is annoying.
	if (is.null(inputprobs_for_fast))
		{ # Return the full Pmat (slow)
		#######################################################
		# Return the Pmat
		#######################################################
		
		# Zero rows will crash the FORTRAN wrapalldgexpv function, and
		# therefore crash R.  This is annoying.
		if (check_for_0_rows != FALSE)
			{
			# If not false, check_for_0_rows is either TRUE or numeric

			# Get T/F for rows with all zeros
			rows_w_all_zeros_TF = findrows_w_all_zeros(matvec)
			
			# If all FALSE, do nothing
			if (all(rows_w_all_zeros_TF == FALSE))
				{
				# Do nothing
				pass = 1
				} else {
				# indices of TRUE
				rows_allzero_indices = seq(1, length(rows_w_all_zeros_TF), 1)[rows_w_all_zeros_TF]

				# Here, you have to input a small value for each zero
				if (is.numeric(check_for_0_rows))
					{
					check_for_0_rows = check_for_0_rows
					} else {
					# 1e-15 appears to be the precision limit of the FORTRAN code
					check_for_0_rows = 1e-15
					}
				# Input the given value into all zeros
				newrowvals = rep(check_for_0_rows, ncol(matvec))
				matvec[rows_allzero_indices, ] = newrowvals
				diagonal_val = -1 * sum(newrowvals[-1])
				matvec[rows_allzero_indices, rows_allzero_indices] = diagonal_val
				
				cat("\nWARNING: ", sum(rows_w_all_zeros_TF), " rows of the Q matrix Qmat had all zeros. This will crash .Call('wrapalldgexpv_', ...)\nand therefore expokit_wrapalldgexpv_tvals() run with the inputprobs_for_fast=NULL option (producing a full Pmat),\nand therefore R.  Replacement value for 0:  check_for_0_rows=", check_for_0_rows, ".\n", sep="")
				}
			}
		}
		
	
	# Count the number of NON-zeros (nz)
	# and input the matrix size
	if (transform_to_coo_TF == TRUE)
		{
		# COO format
		# http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29

		# number of non-zeros
		nz  = sum(matvec != 0)
		
		# These vectors have that length
		ia  = integer(length=nz)
		ja  = integer(length=nz)
		a   = double(length=nz)	
		n=nrow(matvec)
		
		} else {
		n = coo_n
		# (And make a regular matrix from COO)

		# number of non-zeros
		# Assumes that coo-formatted matrix columns are
		# ia, ja, a
		nz  = sum(matvec[,"a"] != 0)
		}

	# ideg = degree of polynomial, 6 is usually considered sufficient
	ideg = as.integer(6)
	#n=nrow(Qmat)	# don't do here, possibly coo-formatted
	m=n-1
	# t=as.numeric(2.1)
	
	# v should have as many elements as n; first element = 1 (?)
	if (is.null(inputprobs_for_fast))
		{
		# Input space-fillers, these get written over by wrapall
		v=double(n)
		v[1] = 1
		# Input the input probabilities, these get used directly by myDGEXPV/myDGEXPV
		} else {
		v = double(n)
		v = inputprobs_for_fast
		}
	
	# w is the same length
	w = double(length=n)
	tol=as.numeric(0.01)
	
	# lwsp = length of wsp
	# wsp = workspace to hold various variables, cells of the matrix, etc.
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 4*(m+2)^2+ideg+1)
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 5*(m+2)^2+ideg+1)
	lwsp = as.integer(n*(m+2)+5*(m+2)^2+ideg+1)
	
	#lwsp = 100
	wsp = double(length=lwsp)
	
	# length of iwsp
	liwsp = max(m+2, 7)
	iwsp = integer(length=liwsp)
	
	#matvec = matrix(data=Q, nrow=n, byrow=TRUE)
	matvec = matvec
	
	# Don't transform if already coo
	if ((transform_to_coo_TF == TRUE) && (transpose_needed == TRUE))
		{
		tmatvec = t(matvec)
		} else {
		tmatvec = matvec
		}
	#rowSums(tmatvec)
	#colSums(tmatvec)
	
	# This might (?) get slow with large matrices -- doesn't seem to
	if ((exists("anorm") == FALSE) || is.null(anorm))
		{
		# Use the 1-norm or one-norm
		if (transform_to_coo_TF==FALSE && transpose_needed==FALSE)
			{
			tmpQmat1 = coo2mat(matvec, n=coo_n)
			tmpQmat2 = t(tmpQmat1)
			anorm = as.numeric(norm(tmpQmat2, type="O"))
			} else {
			anorm = as.numeric(norm(matvec, type="O"))
			}
		}
	
	# The itrace flag, if set to 1, results in dgexpv printing some details of
	# the function's run to screen.
	itrace = 0
	iflag = 0	
	
	# Make the input COO matrix
	# COO = coordinate list format, useful for sparse matrices with lots of zeros:
	# http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29
	# ia = rownum in the matrix
	# ja = colnum in the matrix
	# a  = value of that cell
	
	if (transform_to_coo_TF == TRUE)
		{		
		tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmatvec)
		} else {
		tmpmat_in_REXPOKIT_coo_fmt = matvec
		}
	# Either way, store the rows/columns in the input variables for FORTRAN
	ia = tmpmat_in_REXPOKIT_coo_fmt[,"ia"]
	ja = tmpmat_in_REXPOKIT_coo_fmt[,"ja"]
	a = tmpmat_in_REXPOKIT_coo_fmt[,"a"]

	# Run the wrapper function	
	if (is.null(inputprobs_for_fast))
		{
		######################################
		# Return the full Pmat (slow)
		######################################
		
		# Create the space for res (the returned Pmat)
		res = double(length=n*n)
		
		res2 <- .C("wrapalldgexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res))
	
		# wrapalldgexpv_ returns all kinds of stuff, list item 18 is the P matrix
		# However, this may be an inefficient use of the dgexpv sparse matrix capabilities (Hansen)
		# Try mydgexpv_ to just get the ancestral probabilities (w, res2[[5]])
		output_Pmat = matrix(res2[[18]], nrow=n, byrow=TRUE)
		
		return(output_Pmat)
		} else {
		######################################
		# Instead of returning the full Pmat (slow), just return the output probabilities (fast)
		######################################
		
		# Be sure to input the input probabilities
		v = inputprobs_for_fast
		
		res2 <- .C("mydgexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz))
		
		# w, list item #5, contains the output probabilities
		w_output_probs = res2[[5]]
		
		return(w_output_probs)
		}
	}







#' EXPOKIT dgexpv wrapper function
#'
#' This function wraps the .C call to EXPOKIT for the dgexpv function.  Only the output probability
#' matrix is returned.
#'
#' NOTE: DGEXPV vs. DMEXPV. According to the EXPOKIT documentation, DGEXPV should be
#' faster than DMEXPV, however DMEXPV runs an accuracy check appropriate for
#' Markov chains, which is not done in DGEXPV.
#' 
#' @param n number of rows in Q matrix
#' @param m n-1
#' @param timeval the value to exponentiate the rate matrix by (often e.g. a time value)
#' @param v variable to store some results in; should have n elements (and perhaps start with 1)
#' @param w same length as v
#' @param tol tolerance for approximations; usually set to 0.01
#' @param anorm the norm of the Q matrix
#' @param lwsp length of workspace (wsp); for dgexpv, lwsp=n*(m+2)+5*(m+2)^2+ideg+1
#' @param wsp workspace to store some results in; should be a double with lwsp elements
#' @param liwsp length of integer workspace; for dgexpv, liwsp=m+2
#' @param iwsp integer workspace
#' @param itrace option, set to 0
#' @param iflag option, set to 0
#' @param ia i indices of Qmat nonzero values
#' @param ja j indices of Qmat nonzero values
#' @param a nonzero values of Qmat (ia, ja, a are columns of a COO-formatted Q matrix)
#' @param nz number of non-zeros in Qmat
#' @param res space for output probability matrix (n x n)
#'
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' COO format is required for EXPOKIT.
#' @return \code{tmpoutmat} the output matrix for the (first) input t-value
#' @seealso \code{\link{expokit_dgexpv_Qmat}}
#' @seealso \code{\link{expokit_wrapalldgexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples # Example building the inputs from scratch:
#'
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' timeval = tvals[2]
#'
#' 	ideg = as.integer(6)
#'	n=nrow(Qmat)
#'	m=n-1
#'	# t=as.numeric(2.1)
#'	
#'	# v should have as many elements as n; first element = 1 (?)
#'	v=double(n)
#'	v[1] = 1
#'	
#'	# w is the same length
#'	w = double(length=n)
#'	tol=as.numeric(0.01)
#'	
#'	# length of wsp
#'	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 4*(m+2)^2+ideg+1)
#'	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 5*(m+2)^2+ideg+1)
#'	lwsp = as.integer(n*(m+2)+5*(m+2)^2+ideg+1)
#'	
#'	#lwsp = 100
#'	wsp = double(length=lwsp)
#'	
#'	# length of iwsp
#'	liwsp = max(m+2, 7)
#'	iwsp = integer(length=liwsp)
#'	
#'	res = double(length=n*n)
#'	
#'	#matvec = matrix(data=Q, nrow=n, byrow=TRUE)
#'	matvec = Qmat
#'	tmatvec = t(matvec)
#'	rowSums(tmatvec)
#'	colSums(tmatvec)
#'	
#'	# type="O" is being used here, this is supposed to be the
#'	# default for norm(), although it throws an error if not
#'	# specified
#'	# 
#'	# From the help:
#'	# type - character string, specifying the type of matrix norm to be
#'	# computed. A character indicating the type of norm desired. 
#'	# 	"O", "o" or "1"
#'	# 		specifies the one norm, (maximum absolute column sum);
#'	anorm = as.numeric(norm(matvec, type="O"))
#'	#anorm = 1
#'	
#'	
#'	itrace = 0
#'	iflag = 0	
#'	
#'	
#'	#a = as.numeric(tmatvec)
#'	#a = as.numeric(matvec)
#'	tmpmat = tmatvec
#'	tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmpmat)
#'	ia = tmpmat_in_REXPOKIT_coo_fmt[,"ia"]
#'	ja = tmpmat_in_REXPOKIT_coo_fmt[,"ja"]
#'	a = tmpmat_in_REXPOKIT_coo_fmt[,"a"]
#' 
#'  # Number of non-zeros
#'  nz = nrow(Qmat) * ncol(Qmat)
#' 
#' # Run the wrapper function	
#' 
#' tmpoutmat = expokit_dgexpv_wrapper(n, m, timeval, v, w, tol, anorm, wsp, 
#' lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)
#' 
#' print(tmpoutmat)
#'
expokit_dgexpv_wrapper <- function(n, m, timeval, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)
	{
	res2 = NULL
	
	res2 <- .C("wrapalldgexpv_", as.integer(n), as.integer(m), as.double(timeval), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res))
	
	tmpoutmat = matrix(res2[[18]], nrow=n, byrow=TRUE)
	
	return(tmpoutmat)
	}

	


#' Run EXPOKIT's dgexpv on one or more t-values
#'
#' The function runs EXPOKIT's \code{dgexpv} function on a Q matrix and \emph{one or more} time values.  If \code{Qmat} is NULL (default), a default matrix is input.\cr
#'
#' NOTE: DGEXPV vs. DMEXPV. According to the EXPOKIT documentation, DGEXPV should be
#' faster than DMEXPV, however DMEXPV runs an accuracy check appropriate for
#' Markov chains, which is not done in DGEXPV.\cr
#' 
#' @param Qmat an input Q transition matrix
#' @param tvals one or more time values to exponentiate by (doesn't have to literally be a time value, obviously)
#' @param inputprobs_for_fast If NULL (default), the full probability matrix (Pmat) is returned. However, the full
#' speed of EXPOKIT on sparse matrices will be exploited if inputprobs_for_fast=c(starting probabilities). In this case
#' these starting probabilities are input to \code{myDMEXPV} directly, as \code{v}, and \code{w}, the output probabilities,
#' are returned.
#' @param transpose_needed If TRUE (default), matrix will be transposed (apparently
#' EXPOKIT needs the input matrix to be transposed compared to normal)
#' @param transform_to_coo_TF Should the matrix be tranposed to COO?  COO format is required
#' for EXPOKIT's sparse-matrix functions (like dmexpv and unlike the padm-related 
#' functions. Default TRUE; if FALSE, user must put a COO-formated matrix in \code{Qmat}. Supplying the
#' coo matrix is probably faster for repeated calculations on large matrices.
#' @param coo_n If a COO matrix is input, \code{coo_n} specified the order (# rows, equals # columns) of the matrix.
#' @param force_list_if_1_tval Default FALSE, but set to TRUE if you want a single matrix to be returned
#' inside a list
#' @param check_for_0_rows If TRUE or a numeric value, the input Qmat is checked for all-zero rows, since these will crash the FORTRAN wrapalldmexpv function. A small nonzero value set to check_for_0_rows or the default (0.0000000000001) is input to  off-diagonal cells in the row (and the diagonal value is normalized), which should fix the problem.
#' @return \code{tmpoutmat} the output matrix, if 1 t-value is input; \code{list_of_matrices_output},
#' if more than 1 t-value is input; to get a single output matrix in a list, set \code{force_list_if_1_tval=TRUE}
#' @seealso \code{\link{expokit_dgexpv_wrapper}}
#' @seealso \code{\link{expokit_dgexpv_Qmat}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples # Example:
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' 
#' # Exponentiate each with EXPOKIT's dgexpv (should be fast for large sparse matrices)
#' for (t in tvals)
#' 	{
#' 	Pmat = expokit_dgexpv_Qmat(Qmat=Qmat, t=t, transpose_needed=TRUE)
#' 	cat("\n\nTime=", t, "\n", sep="")
#' 	print(Pmat)
#' 	}
#'
#' # DMEXPV and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
#' # DMEXPV is specifically designed for Markov chains and so may be slower, but more accurate.
#' 
#' # DMEXPV, single t-value
#' 
#' # DGEXPV, single t-value
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=2)
#' 
#' # These functions runs the for-loop itself (sadly, we could not get mapply() to work
#' # on a function that calls dmexpv/dgexpv), returning a list of probability matrices.
#' 
#' # DGEXPV functions
#' list_of_P_matrices_dgexpv = expokit_wrapalldgexpv_tvals(Qmat=Qmat, 
#' tvals=tvals, transpose_needed=TRUE)
#' list_of_P_matrices_dgexpv
#'
expokit_wrapalldgexpv_tvals <- function(Qmat=NULL, tvals=c(2.1), inputprobs_for_fast=NULL, transpose_needed=TRUE, transform_to_coo_TF=TRUE, coo_n=NULL, force_list_if_1_tval=FALSE, check_for_0_rows=TRUE)
	{
	defaults = '
	Qmat=NULL
	t = 2.1
	inputprobs_for_fast=NULL
	transpose_needed=TRUE
	COO_needed=TRUE
	transform_to_coo_TF=TRUE
	coo_n=NULL
	force_list_if_1_tval=FALSE
	check_for_0_rows=FALSE
	check_for_0_rows=1e-15
	'
	
	# Check if Qmat is blank
	if (is.null(Qmat))
		{
		# Default Qmat
		warning("You supplied no matrix, so a default matrix is being used. Obviously you can't use this for anything real. YOU HAVE BEEN WARNED!!")
		
		Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
		}
	
	if (is.null(tvals))
		{
		warning("You supplied no time values, so default time values are being used. Obviously you can't use this for anything real. YOU HAVE BEEN WARNED!!")
		tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
		}

	matvec = Qmat
	
	# Zero rows will crash the FORTRAN wrapalldgexpv function, and
	# therefore crash R.  This is annoying.
	if (is.null(inputprobs_for_fast))
		{ # Return the full Pmat (slow)
		#######################################################
		# Return the Pmat
		#######################################################
		
		# Zero rows will crash the FORTRAN wrapalldgexpv function, and
		# therefore crash R.  This is annoying.
		if (check_for_0_rows != FALSE)
			{
			# If not false, check_for_0_rows is either TRUE or numeric

			# Get T/F for rows with all zeros
			rows_w_all_zeros_TF = findrows_w_all_zeros(matvec)
			
			# If all FALSE, do nothing
			if (all(rows_w_all_zeros_TF == FALSE))
				{
				# Do nothing
				pass = 1
				} else {
				# indices of TRUE
				rows_allzero_indices = seq(1, length(rows_w_all_zeros_TF), 1)[rows_w_all_zeros_TF]

				# Here, you have to input a small value for each zero
				if (is.numeric(check_for_0_rows))
					{
					check_for_0_rows = check_for_0_rows
					} else {
					# 1e-15 appears to be the precision limit of the FORTRAN code
					check_for_0_rows = 1e-15
					}
				# Input the given value into all zeros
				newrowvals = rep(check_for_0_rows, ncol(matvec))
				matvec[rows_allzero_indices, ] = newrowvals
				diagonal_val = -1 * sum(newrowvals[-1])
				matvec[rows_allzero_indices, rows_allzero_indices] = diagonal_val
				
				cat("\nWARNING: ", sum(rows_w_all_zeros_TF), " rows of the Q matrix Qmat had all zeros. This will crash .Call('wrapalldgexpv_', ...)\nand therefore expokit_wrapalldgexpv_tvals() run with the inputprobs_for_fast=NULL option (producing a full Pmat),\nand therefore R.  Replacement value for 0:  check_for_0_rows=", check_for_0_rows, ".\n", sep="")
				}
			}
		}	

		
	# Count the number of NON-zeros (nz)
	# and input the matrix size
	if (transform_to_coo_TF == TRUE)
		{
		# COO format
		# http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29

		# number of non-zeros
		nz  = sum(matvec != 0)
		
		# These vectors have that length
		ia  = integer(length=nz)
		ja  = integer(length=nz)
		a   = double(length=nz)	
		n=nrow(matvec)
		} else {
		n = coo_n
		# (And make a regular matrix from COO)

		# number of non-zeros
		# Assumes that coo-formatted matrix columns are
		# ia, ja, a
		nz  = sum(matvec[,"a"] != 0)
		}
	
	

	#######################################################
	ideg = as.integer(6)
	#######################################################
	#n=nrow(Qmat)	# don't do this here, you might have a coo matrix
	m=n-1
	# t=as.numeric(2.1)
	
	# v should have as many elements as n; first element = 1 (?)
	if (is.null(inputprobs_for_fast))
		{
		# Input space-fillers, these get written over by wrapall
		v=double(n)
		v[1] = 1
		# Input the input probabilities, these get used directly by myDMEXPV/myDGEXPV
		} else {
		v = double(n)
		v = inputprobs_for_fast
		}
	
	# w is the same length
	w = double(length=n)
	tol=as.numeric(0.01)
	
	# length of wsp
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 4*(m+2)^2+ideg+1)
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 5*(m+2)^2+ideg+1)
	lwsp = as.integer(n*(m+2)+5*(m+2)^2+ideg+1)
	
	#lwsp = 100
	wsp = double(length=lwsp)
	
	# length of iwsp
	liwsp = max(m+2, 7)
	iwsp = integer(length=liwsp)
	

	if ((transform_to_coo_TF == TRUE) && (transpose_needed == TRUE))
		{
		tmatvec = t(matvec)
		#rowSums(tmatvec)
		#colSums(tmatvec)
		}


	
	# type="O" is being used here, this is supposed to be the
	# default for norm(), although it throws an error if not
	# specified
	# 
	# From the help:
	# type - character string, specifying the type of matrix norm to be
	# computed. A character indicating the type of norm desired. 
	# 	"O", "o" or "1"
	# 		specifies the one norm, (maximum absolute column sum);
	if ((exists("anorm") == FALSE) || is.null(anorm))
		{
		# Use the 1-norm or one-norm
		if (transform_to_coo_TF==FALSE && transpose_needed==FALSE)
			{
			tmpQmat1 = coo2mat(matvec, n=coo_n)
			tmpQmat2 = t(tmpQmat1)
			anorm = as.numeric(norm(tmpQmat2, type="O"))
			} else {
			anorm = as.numeric(norm(matvec, type="O"))
			}
		}	
	
	itrace = 0
	iflag = 0	
	
	
	#a = as.numeric(tmatvec)
	#a = as.numeric(matvec)
	if (transform_to_coo_TF == TRUE)
		{		
		tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmatvec)
		} else {
		tmpmat_in_REXPOKIT_coo_fmt = matvec
		}
	# Either way, store the rows/columns in the input variables for FORTRAN
	ia = tmpmat_in_REXPOKIT_coo_fmt[,"ia"]
	ja = tmpmat_in_REXPOKIT_coo_fmt[,"ja"]
	a = tmpmat_in_REXPOKIT_coo_fmt[,"a"]

	num_tvals = length(tvals)


	# Run the wrapper function	
	if (is.null(inputprobs_for_fast))
		{ # Return the full Pmat (slow)
		#######################################################
		# Return the Pmat
		#######################################################
		
		# Create the space for res (the returned Pmat)
		res = double(length=n*n)

	
		# If there is more than 1 t-value, or if the user desires a list even for a single
		# t-value, return a list
		if ((num_tvals > 1) || (force_list_if_1_tval==TRUE))
			{

			# Loop through the list of tvals, get the prob. matrix for each
			# sadly, mapply() etc. crash when tried on expokit_dgexpv_wrapper
			
			# Set up empty matrix
			NA_matrix = matrix(NA, nrow=n, ncol=n)
			
			# Set up list of empty matrices
			list_of_matrices_output = replicate(NA_matrix, n=num_tvals, simplify=FALSE)
			
			for (i in 1:num_tvals)
				{
				timeval = tvals[i]
				list_of_matrices_output[[i]] = expokit_dgexpv_wrapper(n, m, timeval, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)
	
				} # end forloop
			
			return(list_of_matrices_output)
			
			} else {
			# If there is only 1 t value, just return 1 matrix
			#res2 <- .C("wrapalldgexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res))
			
			#tmpoutmat = matrix(res2[[18]], nrow=n, byrow=TRUE)
	
			timeval=tvals
			output_Pmat = expokit_dgexpv_wrapper(n, m, timeval, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)		
			
			#print(tmpoutmat)
			return(output_Pmat)
			}
		} else {
		#######################################################
		# Return the output probabilities
		#######################################################
				
		# Be sure to input the input probabilities
		v = inputprobs_for_fast

		# If there is more than 1 t-value, or if the user desires a list even for a single
		# t-value, return a list
		if ((num_tvals > 1) || (force_list_if_1_tval==TRUE))
			{
			# Loop through the list of tvals, get the prob. matrix for each
			# sadly, mapply() etc. crash when tried on expokit_dgexpv_wrapper
			
			# Set up empty matrix
			list_of_outprobs_output = matrix(NA, nrow=num_tvals, ncol=n)

			for (i in 1:num_tvals)
				{
				t = tvals[i]
			
				# This must be mydgexpv_, not myDGEXPV_ !!!!

				res2 <- .C("mydgexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz))
		
				# w, list item #5, contains the output probabilities
				w_output_probs = res2[[5]]

				list_of_outprobs_output[i,] = w_output_probs
	
				}
			} else {
			
			# If there is only 1 t value, just return 1 matrix
			#res2 <- .C("wrapalldgexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res))
			
			#tmpoutmat = matrix(res2[[18]], nrow=n, byrow=TRUE)
	
			t=tvals

			# This must be mydgexpv_, not myDGEXPV_ !!!!

			res2 <- .C("mydgexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz))
			
			# w, list item #5, contains the output probabilities
			w_output_probs = res2[[5]]

			list_of_outprobs_output[1,] = w_output_probs
	
			return(w_output_probs)
			}

		return(list_of_matrices_output)
		}
	}
