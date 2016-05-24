#' @include fermat.R
 
#roxygenize()





#######################################################
# Utilities
#######################################################

#' String splitting shortcut
#' 
#' \code{\link[base]{strsplit}} returns the results inside a list, which is annoying. \code{strsplit3} shortens the process.
#'
#' @param x A string to split
#' @param ... Other arguments to \code{\link[base]{strsplit}}.  The argument \code{split} is \emph{required}.
#' @return \code{out} The output from inside the list.
#' @export
#' @seealso \code{\link[base]{strsplit}}
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' 
#' # strsplit returns the results inside a list element
#' out = strsplit("ABC", split="")
#' out
#' # I.e....
#' out[[1]]
#' 
#' # If this is annoying/ugly in the code, use strsplit3:
#' out = strsplit3("ABC", split="")
#' out
#' 
strsplit3 <- function(x, ...)
	{
	out = strsplit(x, ...)[[1]]
	return(out)
	}

# Note: Rd files barely understand LaTeX, so just do basics. The PDF files will do
# nice LaTeX, however. 
#######################################################
# numstates_from_numareas:
#######################################################
#' Calculate the number of states, given a certain number of areas
#' 
#' This function calculates the number of discrete states that are needed to 
#' represent the possible combinations of presence and absence in a set of 
#' discrete areas.  The number of states is a function of the number of areas,
#' and the maximum allowed range size (in number of areas) of a species.
#' 
#' For example, with 3 areas (A, B, C), there are 8 possible states, if a null
#' range is allowed (null, A, B, C, AB, BC, AC, ABC).  If the maximum range size is
#' 2 areas, then there are only 7 possible states.
#' 
#' The formula for the number of geographic states, based on the number of areas (\emph{N}),
#' is the sum of \emph{N} choose \emph{k}, from \emph{k}=1 to \emph{m}
#' (maximum range size) \deqn{s = \sum_{k=1}^{m}{N\choose k}}{s = sum(k=1...m)(N choose k)}
#' 
#' This equation assumes that the null range (a species lives in no areas, i.e. is extinct)
#' is not allowed. In the LAGRANGE program of \cite{ReeSmith2008}), the null range is included
#' in the transition matrix, and thus this is one more state.  This situation is represented in 
#' \code{numstates_from_numareas} by setting \code{include_null_range=TRUE}.
#'
#' Users might manually remove states from the states list, if prior information indicates that
#' some configurations of presence/absence in areas are impossible as geographic ranges for
#' species.  If so, they should manually subtract from the number of states.
#'
#' @param numareas The number of areas in the analysis.
#' @param maxareas The maximum number of areas that any single species/lineage can occupy.
#' @param include_null_range If FALSE (default), the null range is not included in the count.
#' If TRUE, the null range is included, adding +1 to the count of the states.
#' @return \code{nstates} Number of states
#' @export
#' @seealso \code{\link[stats]{convolve}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' numstates_from_numareas(numareas=3, maxareas=3, 
#' include_null_range=FALSE)
#' numstates_from_numareas(numareas=3, maxareas=3, 
#' include_null_range=TRUE)
#' numstates_from_numareas(numareas=3, maxareas=2, 
#' include_null_range=TRUE)
#' numstates_from_numareas(numareas=3, maxareas=1, 
#' include_null_range=TRUE)
#' numstates_from_numareas(numareas=7, maxareas=7, 
#' include_null_range=TRUE)
#' numstates_from_numareas(numareas=7, maxareas=2, 
#' include_null_range=TRUE)
#' numstates_from_numareas(numareas=8, maxareas=8, 
#' include_null_range=TRUE)
#' numstates_from_numareas(numareas=8, maxareas=2, 
#' include_null_range=TRUE)
#' numstates_from_numareas(numareas=20, maxareas=20, 
#' include_null_range=TRUE)
#' numstates_from_numareas(numareas=20, maxareas=2, 
#' include_null_range=TRUE)
#' numstates_from_numareas(numareas=20, maxareas=3, 
#' include_null_range=TRUE)
numstates_from_numareas <- function(numareas=3, maxareas=numareas, include_null_range=FALSE)
	{
	# The formula for the number of geographic states, based on the number of areas,
	# is:
	#
	# sum_from_k=1_to_m (N choose k)
	# 
	nstates = 0
	
	for (i in 1:maxareas)
		{
		tmp_nstates = choose(n=numareas, k=i)
		nstates = nstates + tmp_nstates
		}
	
	if (include_null_range == TRUE)
		{
		nstates = nstates + 1
		}
	
	return(nstates)
	}




#######################################################
# areas_list_to_states_list_old
#######################################################
#' Convert a list of areas to a list of geographic ranges (states); original R version
#' 
#' This is the original R version of the function which converts a list of possible areas to
#' a list of all possible states (geographic ranges).  This gets slow for large numbers of areas.
#' 
#' The function is mostly replaced by \code{\link[cladoRcpp]{rcpp_areas_list_to_states_list}} in optimized code, but is still used in some places
#' for display purposes.
#' 
#' @param areas a list of areas (character or number; the function converts these to numbers, starting with 0)
#' @param maxareas maximum number of areas in this analyses
#' @param include_null_range \code{TRUE} or \code{FALSE}, should the \code{NULL} range be included in the possible states? (e.g., LAGRANGE default is yes)
#' @param split_ABC \code{TRUE} or \code{FALSE} If \code{TRUE} the output will consist of a list of lists (c("A","B","C"), c("A","B"), c("A","D"), etc.); 
#' if \code{FALSE}, the list of areas will be collapsed ("ABC", "AB", "AD", etc.).
#' @return \code{states_list} A list of the states.
#' @export
#' @seealso \code{\link{numstates_from_numareas}}, \code{\link{rcpp_areas_list_to_states_list}}
#' @note No notes.
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @examples
#' areas = c("A","B","C")
#' areas_list_to_states_list_old(areas=areas, maxareas=length(areas), 
#' include_null_range=TRUE, split_ABC=TRUE)
#' areas_list_to_states_list_old(areas=areas, maxareas=length(areas), 
#' include_null_range=TRUE, split_ABC=FALSE)
#' areas_list_to_states_list_old(areas=areas, maxareas=length(areas), 
#' include_null_range=FALSE, split_ABC=TRUE)
#' areas_list_to_states_list_old(areas=areas, maxareas=length(areas), 
#' include_null_range=FALSE, split_ABC=FALSE)
#' areas_list_to_states_list_old(areas=areas, maxareas=2, 
#' include_null_range=TRUE, split_ABC=TRUE)
#' areas_list_to_states_list_old(areas=areas, maxareas=2, 
#' include_null_range=TRUE, split_ABC=FALSE)
#' areas_list_to_states_list_old(areas=areas, maxareas=2, 
#' include_null_range=FALSE, split_ABC=TRUE)
#' areas_list_to_states_list_old(areas=areas, maxareas=2, 
#' include_null_range=FALSE, split_ABC=FALSE)
#' areas_list_to_states_list_old(areas=areas, maxareas=1, 
#' include_null_range=TRUE, split_ABC=TRUE)
#' areas_list_to_states_list_old(areas=areas, maxareas=1, 
#' include_null_range=TRUE, split_ABC=FALSE)
#' areas_list_to_states_list_old(areas=areas, maxareas=1, 
#' include_null_range=FALSE, split_ABC=TRUE)
#' areas_list_to_states_list_old(areas=areas, maxareas=1, 
#' include_null_range=FALSE, split_ABC=FALSE)
#' 
areas_list_to_states_list_old <- function(areas=c("A","B","C"), maxareas=length(areas), include_null_range=TRUE, split_ABC=TRUE)
	{
	
	# Error trap
	if (maxareas > length(areas))
		{
		maxareas = length(areas)
		}
	
	
	# Initialize the states_list to the correct size
	nstates = numstates_from_numareas(numareas=length(areas), maxareas=maxareas, include_null_range=include_null_range)
	
	states_list = rep(NA, times=nstates)

	# Add null range (globally extinct) to the states list
	if (include_null_range == TRUE)
		{
		# Start the index at 1 (first entry will be s=2)
		s = 1
		} else {
		s = 0
		}
	
	if (split_ABC == FALSE)
		{
		# Option #1: Don't split states
		# Add range combinations to the list
		for (m in 1:maxareas)
			{
			states_matrix = combn(x=areas, m=m)
			
			# Collapse the columns into txt
			tmp_states = apply(X=states_matrix, 2, paste, collapse="")
			
			# Add the states to the list
			for (i in 1:length(tmp_states))
				{
				states_list[[(s=s+1)]] = tmp_states[[i]]
				}
			}
		}
	
	# Unlist states_list
	#states_list = unlist(states_list)



	if (split_ABC == TRUE)
		{

		# Option #2: Do split states
		# Add range combinations to the list
		for (m in 1:maxareas)
			{
			states_matrix = combn(x=areas, m=m)
			
			# Collapse the columns into txt
			tmp_states = apply(X=states_matrix, 2, paste, collapse=",")
			
			# Add the states to the list
			for (i in 1:length(tmp_states))
				{
				states_list[[(s=s+1)]] = tmp_states[[i]]
				}
			}

		# Split the txt
		states_list = mapply(FUN=strsplit3, split=",", states_list)
		names(states_list) = NULL		
		}

			
	# Add null range (globally extinct) to the states list
	if (include_null_range == TRUE)
		{
		s = 1
		states_list[[s]] = c("_")
		} else {
		s = 1
		#states_list[[s]] = c(NULL)		
		}
	return(states_list)
	}



#######################################################
# rcpp_areas_list_to_states_list
#######################################################
#' Make a list of 0-based indices of possible combinations of input areas
#'
#' Given a list of areas (actually a list of anything; all that is important is the length of the list)
#' \code{rcpp_areas_list_to_states_list} calculates all possible combinations of these areas,
#' listing them by the 0-based indices that specify the position of each area in the list.
#'
#' Using 0-based indexing is convenient in the C++ code called by the other functions, rather than having
#' to keep track of the various people might label their areas (names, abbreviations, letters, numbers).
#' 
#' As in LAGRANGE (Ree & Smith 2008), the maximum range size (i.e. the maximum number of areas in a range) can be specified by the 
#' user.  Having a smaller maximum range size drastically reduces the number of states, and thus the size of
#' the transition matrix and the cladogenesis matrix.
#'
#' @param areas a list of areas (character or number; the function converts these to numbers, starting with 0)
#' @param maxareas maximum number of areas in this analyses
#' @param include_null_range \code{TRUE} or \code{FALSE}, should the \code{NULL} range be included in the possible states? (e.g., \code{LAGRANGE} default is yes)
#' @return \code{R_states_list} A list of the states, where each state is a list of areas in the form of 0-based indices
#' @export
#' @seealso \code{\link{numstates_from_numareas}}, \code{\link{areas_list_to_states_list_old}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Specify the areas
#' areas_list = c("A", "B", "C")
#' areas_list
#' 
#' # Let's try Rcpp_combn_zerostart, in case that is the source of a
#' # problem found via AddressSanitizer
#' Rcpp_combn_zerostart(n_to_choose_from=4, k_to_choose=2, maxlim=1e+07)
#' Rcpp_combn_zerostart(n_to_choose_from=4, k_to_choose=3, maxlim=1e+07)
#'
#'
#' \dontrun{
#' 
#' # Calculate the list of 0-based indices for each possible geographic range, i.e.
#' # each combination of areas
#' states_list = rcpp_areas_list_to_states_list(areas=areas_list, maxareas=3, 
#' include_null_range=FALSE)
#' states_list
#' states_list = rcpp_areas_list_to_states_list(areas=areas_list, maxareas=3, 
#' include_null_range=TRUE)
#' states_list
#' states_list = rcpp_areas_list_to_states_list(areas=areas_list, maxareas=2, 
#' include_null_range=TRUE)
#' states_list
#' states_list = rcpp_areas_list_to_states_list(areas=areas_list, maxareas=1, 
#' include_null_range=TRUE)
#' states_list
#' 
#' }
#' 
rcpp_areas_list_to_states_list <- function(areas=c("A","B","C"), maxareas=length(areas), include_null_range=TRUE)
	{
	# Error trap
	if (maxareas > length(areas))
		{
		maxareas = length(areas)
		}
	
	# Initialize the states_list to the correct size
	nstates = numstates_from_numareas(numareas=length(areas), maxareas=maxareas, include_null_range=include_null_range)
	
	R_areas_indices = seq(0, length(areas)-1, 1)
	R_states_list = as.list(rep(-1, times=nstates))
	R_maxareas = maxareas
	R_include_null_range = include_null_range


	# Call the fast C++ function
	#R_states_list = .Call( "cpp_areas_list_to_states_list", R_areas_indices=R_areas_indices, R_maxareas=R_maxareas, R_include_null_range=R_include_null_range, PACKAGE = "cladoRcpp" )
	# n = R_areas_indices
# 	m = 2
# 	nrows = m
# 	ncols = choose(n,m)
# 	combmat = matrix(data=0, nrow=nrows, ncol=ncols)
# 	
# 	length_output = nrows * ncols
# 	

	
	R_states_list = .Call("cpp_areas_list_to_states_list", as.integer(R_areas_indices), as.integer(R_maxareas), as.logical(R_include_null_range))
	
	if (include_null_range == TRUE)
		{
		# First item may be NULL
		if (R_states_list[[1]] == -1)
			{
			R_states_list[[1]] = NA
			}
		}

	return(R_states_list)
	}



#######################################################
# rcpp_states_list_to_DEmat
#######################################################
#' C++ conversion of a states list to a dispersal-extinction matrix (DEmat)
#'
#' This function takes a list of states/ranges, a matrix describing relative dispersal probability (dmat) for each pair of areas, and
#' a list describing the local extirpation probability for each area (elist), and calculates a transition matrix Qmat accordingly.
#'
#' The size of the matrix will expand dramatically with the number of areas.  See \code{\link{numstates_from_numareas}} for the calculation.
#' 
#' Above 7 or so areas, making \code{Qmat} a COO-formatted matrix (COO=Coordinate list, see wikipedia, \url{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29} ) which can then be used in \code{rexpokit}'s sparse-matrix algorithms,
#' should be more efficient. (Sparse matrices are matrices made of mostly 0s.)
#'
#' @param areas_list a list of lists of areas (numbers, starting with 0)
#' @param states_list a list of lists of areas (numbers, starting with 0)
#' @param dmat dispersal matrix from area to area
#' @param elist a list of extinction probabilities
#' @param amat A matrix specifying the probability of instantaneous transition from one area to another (as in standard character rate matrices).
#' @param include_null_range include the null () range (NA) in the matrix (LAGRANGE default=TRUE)
#' @param normalize_TF should the columns be -1 * rowsums?
#' @param makeCOO_TF should the returned matrix be COO or standard dense (the latter is default).
#' @param min_precision what is the effective minimum size for 0
#' @return dmat (a standard Q matrix)
#' @export
#' @seealso \code{\link{numstates_from_numareas}}, \code{\link[stats]{convolve}}
#' @references
#'   \url{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Specify the areas
#' areas_list_txt = c("A", "B", "C")
#' areas_list_txt
#' 
#' # rcpp_states_list_to_DEmat function requires a 0-based list of areas
#' areas_list = seq(0, length(areas_list_txt)-1, 1)
#' areas_list
#' 
#' \dontrun{
#' 
#' # Calculate the list of 0-based indices for each possible 
#' #geographic range, i.e. each combination of areas
#' states_list = rcpp_areas_list_to_states_list(areas=areas_list, maxareas=3, 
#' include_null_range=FALSE)
#' states_list
#' states_list = rcpp_areas_list_to_states_list(areas=areas_list, maxareas=3, 
#' include_null_range=TRUE)
#' states_list
#' states_list = rcpp_areas_list_to_states_list(areas=areas_list, maxareas=2, 
#' include_null_range=TRUE)
#' states_list
#' states_list = rcpp_areas_list_to_states_list(areas=areas_list, maxareas=1, 
#' include_null_range=TRUE)
#' states_list
#' 
#' # Hard-code the along-branch dispersal and extinction rates
#' d = 0.2
#' e = 0.1
#' 
#' # Calculate the dispersal weights matrix and the extinction weights matrix
#' # Equal dispersal in all directions (unconstrained)
#' areas = areas_list
#' distances_mat = matrix(1, nrow=length(areas), ncol=length(areas))
#' dmat = matrix(d, nrow=length(areas), ncol=length(areas))
#' dmat
#' 
#' # Equal extinction probability for all areas
#' elist = rep(e, length(areas))
#' elist
#' 
#' # Set up the instantaneous rate matrix (Q matrix, Qmat)
#' # DON'T force a sparse-style (COO-formatted) matrix here
#' force_sparse = FALSE
#' Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, 
#' include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
#' Qmat
#' 
#' # DO force a sparse-style (COO-formatted) matrix here
#' force_sparse = TRUE
#' Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, 
#' include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
#' Qmat
#' 
#' 
#' # Repeat with an amat
#' amat = dmat
#' amat[is.numeric(amat)] = 0.33
#' 
#' # Set up the instantaneous rate matrix (Q matrix, Qmat)
#' # DON'T force a sparse-style (COO-formatted) matrix here
#' force_sparse = FALSE
#' Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, amat, 
#' include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
#' Qmat
#' 
#' # DO force a sparse-style (COO-formatted) matrix here
#' force_sparse = TRUE
#' Qmat = rcpp_states_list_to_DEmat(areas_list, states_list, dmat, elist, amat, 
#' include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=force_sparse)
#' Qmat
#' }
#' 
#' 

rcpp_states_list_to_DEmat <- function(areas_list, states_list, dmat, elist, amat=NULL, include_null_range=TRUE, normalize_TF=TRUE, makeCOO_TF=FALSE, min_precision=1e-26)
	{
	
	# If amat is NULL, give it 0s
	if (is.null(amat))
		{
		amat = dmat
		amat[is.numeric(amat)] = 0.0
		amat
		}
	
	
	# Initialize the states_list to the correct size
	#nstates = numstates_from_numareas(numareas=length(areas), maxareas=maxareas, include_null_range=include_null_range)
	
	#R_areas_indices = seq(0, length(areas_list)-1, 1)
	R_areas_indices = areas_list
	#R_states_list = as.list(rep(NA, times=nstates))
	R_states_list = states_list
	
	# convert NAs to 
	
	R_include_null_range = include_null_range
	
	# convert NA ranges to -1
	R_include_null_range[is.na(R_include_null_range)] = -1

	# Call the fast C++ function
	#R_states_list = .Call( "cpp_areas_list_to_states_list", R_areas_indices=R_areas_indices, R_maxareas=R_maxareas, R_include_null_range=R_include_null_range, PACKAGE = "cladoRcpp" )
	
	
	# Convert TF to 1/0
	if (normalize_TF == TRUE)
		{
		R_normalize_TF = 1
		} else {
		R_normalize_TF = 0
		}
	# Convert TF to 1/0
	if (makeCOO_TF == TRUE)
		{
		DEmat_COO = .Call("cpp_states_list_to_DEmat_COO", as.list(R_areas_indices), as.list(R_states_list), as.matrix(dmat), as.numeric(elist), as.matrix(amat), as.integer(R_normalize_TF), as.numeric(min_precision))
		
	
		return(as.data.frame(DEmat_COO))
		} else {
		DEmat = .Call("cpp_states_list_to_DEmat", as.list(R_areas_indices), as.list(R_states_list), as.matrix(dmat), as.numeric(elist), as.matrix(amat), as.integer(R_normalize_TF))
		return(DEmat)
		}

	return()
	}



#######################################################
# rcpp_convolve:
#######################################################
#' Run C++ version of convolve(x,y, conj=TRUE, type="open")
#'
#' This function runs a C++ version of the R function \code{\link{convolve}}, 
#' specifically: \code{convolve(x,y, conj=TRUE, type="open")}
#' 
#' The R function \code{\link{convolve}} is an example of an R function that 
#' gets very slow when the input vectors are large. This C++ version, \code{rcpp_convolve}
#' can be dramatically faster for large vectors.
#'
#' \code{rcpp_convolve} produces the same output as: \code{convolve(ca, cb, conj=TRUE, type="open")}
#' 
#' Note: The C++ code is from the Rcpp examples in: Eddelbuettel & Francois (2011). Rcpp: Seamless R and C++ Integration. \emph{Journal of Statistical Software}, 40(8), 1-18.
#'
#' @param a a numeric vector
#' @param b a numeric vector
#' @return \code{convolve_result_vector} the vector which is the product of the convolution
#' @export
#' @seealso \code{\link[Rcpp]{Rcpp}}, \code{\link{convolve}}, \code{\link{rcpp_mult2probvect}}, \code{\link{Rcpp_combn_zerostart}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Eddelbuettel_Francois_2011
#' @author C++ code by: Dirk Eddelbuettel <edd at debian.org> & Romain Francois (2011); This R wrapper & documentation: Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Set up 2 vectors, then convolve them
#' ca = c(1,2,3,4,5)
#' cb = c(2,2,2,2,2)
#' rcpp_convolve(a=ca, b=cb)
#'
#' # Same as:
#' convolve(ca, cb, conj=TRUE, type="open")
#'
rcpp_convolve <- function(a, b){
	convolve_result_vector = .Call( "convolve3cpp", a=a, b=b, PACKAGE = "cladoRcpp" )
	return(convolve_result_vector)
}


#######################################################
# rcpp_mult2probvect
#######################################################
#' Get the product of multiplying each pair of values in a vector (cross-product)
#'
#' This function calls a C++ function which multiplies two vectors by each other elementwise,
#' such that the output is of length(a) * length(b).
#'
#' This is the cross-product operation, which exists in R (\code{\link[base]{\%o\%}} or 
#' \code{\link[base]{tcrossprod}}). However, it is handy to have is as a C++ function
#' for calculating the probability of pairs of descendant states, given the
#' probability of each state individually.
#'
#' @param a a numeric vector
#' @param b a numeric vector
#' @return \code{tcross_product_vector} the vector which is the product of the convolution
#' @export
#' @seealso \code{\link[base]{\%o\%}}, \code{\link[base]{tcrossprod}}, \code{\link{Rcpp_combn_zerostart}}, \code{\link{rcpp_convolve}}
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' ca = c(1,2,3,4,5)
#' cb = c(2,2,2,2,2)
#' rcpp_mult2probvect(a=ca, b=cb)
#'
#' # Same as:
#' c(ca %o% cb)
#' 
#' # Or:
#' c(outer(ca, cb))
#' 
#' # Or:
#' tcrossprod(ca, cb)
#'
rcpp_mult2probvect <- function(a, b){
	tcross_product_vector = .Call( "mult2probvect", a=a, b=b, PACKAGE = "cladoRcpp" )
	return(tcross_product_vector)
}




#######################################################
# Rcpp_combn_zerostart
#######################################################
#' Get all the combinations of descendent state pairs, in 0-based index form
#'
#' Given the number of states, this function returns all of the
#' pairs of indexes corresponding to those states.
#' 
#' The C++ version is MUCH faster than the plain-R version.
#'
#' @param n_to_choose_from N in N choose K
#' @param k_to_choose K in N choose K
#' @param maxlim To avoid memory overruns, the number of combinations can be no larger than \code{maxlim}
#' (default: 1e+07)
#' @return \code{outarray} an integer matrix with \code{outarray} rows; the number of columns is the
#' number of combinations.
#' @export
#' @seealso \code{\link{rcpp_calc_anclikes_sp}}, \code{\link{rcpp_mult2probvect}}, \code{\link{rcpp_convolve}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' Rcpp_combn_zerostart(n_to_choose_from=4, k_to_choose=2, maxlim=1e+07)
#' Rcpp_combn_zerostart(n_to_choose_from=4, k_to_choose=3, maxlim=1e+07)
#'
Rcpp_combn_zerostart <- function(n_to_choose_from, k_to_choose, maxlim=1e+07)
	{
	defaults ='
	n_to_choose_from = 5
	k_to_choose = 6
	maxlim=1e+07
	'
	
	n = n_to_choose_from
	m = k_to_choose
	
	# HEAD OFF MORE ERRORS
	if (n < 1)
		{
		txt = cat("\nERROR: n_to_choose_from must be > 0.  Your value is ", n_to_choose_from, "\n", sep="")
		cat(txt)
		return(stop(txt))
		}

	if (m < 1)
		{
		txt = cat("\nERROR: k_to_choose must be > 0.  Your value is ", k_to_choose, "\n", sep="")
		cat(txt)
		return(stop(txt))
		}

	if (n < m)
		{
		txt = cat("\nERROR: k_to_choose must be less than or equal to n_to_choose_from.\n",
		"Your values are:\n",
		"n_to_choose_from =		", n_to_choose_from, "\n",
		"k_to_choose =		", k_to_choose, "\n", sep="")
		cat(txt)
		return(stop(txt))
		}


	
	# HEAD OFF ERROR
	predicted_number_of_cells_to_fill = choose(n,m)
	
	if ((predicted_number_of_cells_to_fill > maxlim) || (predicted_number_of_cells_to_fill < 1))
		{
		txt = paste("ERROR: n=", n_to_choose_from, ", k=", k_to_choose, ", n choose k=", predicted_number_of_cells_to_fill, " > maxlim=", maxlim, "\nCalculating something this big may crash your computer!; (Note: (n choose k)<1 also raises this error.", sep="")
		cat(txt)
		return(stop(txt))
		}
	
	
	outarray = .Call("cpp_combn_zerostart", as.integer(n), as.integer(m), as.integer(maxlim))

	#R_states_list <- matrix(out$res,nrow=m,byrow=F)

	return(outarray)
	}



#######################################################
# rcpp_calc_anclikes_sp_prebyte:
#######################################################
#' Calculate probability of ancestral states below a speciation event, given probabilities of the states on each descendant branch
#'
#' This is the pre-byte compiled version of \code{\link{rcpp_calc_anclikes_sp}}.  \code{\link{rcpp_calc_anclikes_sp}} is 
#' byte-compiled, which (might) make it faster.  See \code{\link{rcpp_calc_anclikes_sp}} for full description and help.
#'
#' This function gets slow for large state spaces.
#' 
#' For information on byte-compiling, see \url{http://www.r-statistics.com/2012/04/speed-up-your-r-code-using-a-just-in-time-jit-compiler/} and \code{\link[compiler:cmpfun]{cmpfun}} in the \code{\link[compiler:cmpfun]{compiler}} package.
#'
#' @param Rcpp_leftprobs Probabilities of the states at the base of the left descendant branch
#' @param Rcpp_rightprobs Probabilities of the states at the base of the right descendant branch
#' @param l List of state indices (0-based)
#' @param s Relative weight of sympatric "subset" speciation. Default \code{s=1} mimics LAGRANGE model.
#' @param v Relative weight of vicariant speciation. Default \code{v=1} mimics LAGRANGE model.
#' @param j Relative weight of "founder event speciation"/jump speciation. Default \code{j=0} mimics LAGRANGE model.
#' @param y Relative weight of fully sympatric speciation (range-copying). Default \code{y=1} mimics LAGRANGE model.
#' @param dmat If given, a matrix of rank numareas giving multipliers for the probability
#' of each dispersal event between areas. Default NULL, which sets every cell of the 
#' \code{dmat} matrix to value 1.  Users may construct their own parameterized \code{dmat}
#' (for example, making \code{dmat} a function of distance) for inclusion in ML or
#' Bayesian analyses.
#' @param maxent01s Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a subset-sympatric speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent
#' always has range size 1).
#' @param maxent01v Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a vicariance speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param maxent01j Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a founder-event speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param maxent01y Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a full-sympatric (range-copying) 
#' speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param max_minsize_as_function_of_ancsize If given, any state with a range larger that this value 
#' will be given a probability of zero (for the branch with the smaller rangesize).  This means that 
#' not every possible combination of ranges has to be checked, which can get very slow for large 
#' state spaces.
#' @param Rsp_rowsums A vector of size (numstates)  giving the sum of the relative probabilites of 
#' each combination of descendant states, assuming the probabilities of the left- and right-states are 
#' all equal (set to 1). This is thus the sum of the weights, and dividing by this normalization 
#' vector means that the each row of the speciation probability matrix will sum to 1.  Default assumes 
#' the weights sum to 1 but this is not usually the case. Rsp_rowsums need only be calculated once 
#' per tree+model combination, stored, and then re-used for each node in the tree, yielding 
#' significant time savings.
#' @param printmat Should the probability matrix output be printed to screen? (useful for debugging, but 
#' can be dramatically slow in R.app for some reason for even moderate numbers of states; perhaps 
#' overrunning the line length...)
#' @return \code{prob_ancestral_states} The probabilities of the ancestral states.
#' @export
#' @seealso \code{\link{rcpp_calc_anclikes_sp}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # For the basic logic of a probablistic cladogenesis model, see
#' ?rcpp_calc_anclikes_sp
#' 
#' # For examples of running the functions, see the comparison of all functions at:
#' # ?cladoRcpp
#'
rcpp_calc_anclikes_sp_prebyte <- function(Rcpp_leftprobs, Rcpp_rightprobs, l, s=1, v=1, j=0, y=1, dmat=NULL, maxent01s=NULL, maxent01v=NULL, maxent01j=NULL, maxent01y=NULL, max_minsize_as_function_of_ancsize=NULL, Rsp_rowsums=rep(1.0,length(Rcpp_leftprobs)), printmat=FALSE)
	{
	# Error check to avoid segfault
	if (all(c(length(Rcpp_leftprobs)==length(Rcpp_rightprobs), length(Rcpp_leftprobs)==length(l), length(l)==length(Rcpp_rightprobs))) == FALSE)
		{
		error_msg = paste("ERROR: These are not the same length (Rcpp_leftprobs, Rcpp_rightprobs, l): ", length(Rcpp_leftprobs), ", ", length(Rcpp_rightprobs), ", ", length(l), sep="")
		return(stop(error_msg))
		}
	
	
	# Take the max of the 0-based indices of the possible areas, and add 1
	numareas = max(unlist(l), na.rm=TRUE) + 1
	
	# max_minsize_as_function_of_ancsize can just be all 1s as default
	if (is.null(max_minsize_as_function_of_ancsize))
		{
		max_minsize_as_function_of_ancsize = rep(1, numareas)
		}
	
	# If needed, fill in dmat
	if (is.null(dmat))
		{
		dmat = matrix(data=1, nrow=numareas, ncol=numareas)
		}
	
	##################################################################
	# If needed, put e.g. 
	# 1 0 0
	# 1 0 0
	# 1 0 0
	# ...into the maxent distributions
	##################################################################
	
	# subset sympatric speciation
	if (is.null(maxent01s))
		{
		maxent01s = matrix(0, nrow=numareas, ncol=numareas)
		maxent01s[,1] = 1
		}
	
	# vicariant speciation
	# the max range size of the smaller range is 
	# numareas/2, rounded down
	if (is.null(maxent01v))
		{
		numareas_smaller_vicariant_range = floor(numareas / 2)
		#maxent01v = matrix(0, nrow=numareas, ncol=numareas_smaller_vicariant_range)
		maxent01v = matrix(0, nrow=numareas, ncol=numareas)
		maxent01v[, 1] = 1
		#maxent01v = relative_probabilities_of_vicariants(max_numareas=numareas, maxent_constraint_01v=maxent01v_param, NA_val=-1)
		}
	
	# founder-event speciation ("jump dispersal")
	if (is.null(maxent01j))
		{
		maxent01j = matrix(0, nrow=numareas, ncol=numareas)
		maxent01j[, 1] = 1
		}
	
	# full sympatric speciation (range copying)
	# This is essentially subset speciation where the "smaller" range is the
	# same size as the "bigger" range
	if (is.null(maxent01y))
		{
		maxent01y = matrix(0, nrow=numareas, ncol=numareas)
		maxent01y[, 1] = 1
		}
	
	
	
#	tmp = unlist(maxent01s[(maxent01s > 0 && !is.na(maxent01s))])
# 	max_min_rangesize_s = 1:length(maxent01s[maxent01s > 0])
# 	max_min_rangesize_v = 1:length(maxent01v[maxent01v > 0])
# 	max_min_rangesize_j = 1:length(maxent01j[maxent01j > 0])
# 	max_min_rangesize_y = 1:length(maxent01y[maxent01y > 0])
# 	
# 	max_min_rangesize = max(unlist(c(max_min_rangesize_s, max_min_rangesize_v, max_min_rangesize_j)))
	
	# Don't use a, b, etc directly, they get screwed up
	tmpa = Rcpp_leftprobs
	tmpb = Rcpp_rightprobs
	
	# Print the matrix output to screen?	
	if (printmat == TRUE)
		{
		Rprintmat = 1
		} else {
		Rprintmat = 0
		}

	# Call the fast C++ function
	prob_ancestral_states = .Call( "cpp_calc_anclikes_sp", Rprintmat=as.integer(Rprintmat), cpp_leftprobs=as.numeric(tmpa), cpp_leftprobs=as.numeric(tmpb), l=l, s=as.numeric(s), v=as.numeric(v), j=as.numeric(j), y=as.numeric(y), dmat=as.matrix(dmat), maxent01s=as.matrix(maxent01s), maxent01v=as.matrix(maxent01v), maxent01j=as.matrix(maxent01j), maxent01y=as.matrix(maxent01y), max_minsize_as_function_of_ancsize=as.integer(max_minsize_as_function_of_ancsize), Rsp_rowsums=as.numeric(Rsp_rowsums), PACKAGE = "cladoRcpp" )
	
	return(prob_ancestral_states)
	}




#######################################################
# rcpp_calc_anclikes_sp:
#######################################################
#' Calculate probability of ancestral states below a speciation event, given probabilities of the states on each descendant branch
#'
#' This function, given parameters on the Relative weight of different geographic range inheritance
#' scenarios at cladogenesis (speciation) events, calculates the probability of each possible ancestral
#' state given the probabilities of each possible combination of tip states.
#'
#' The Python/C++ program LAGRANGE (Ree & Smith 2008) gives a fixed equal probability to each 
#' range-inheritance scenario it allows:
#'
#' (1) sympatric speciation with 1 area (e.g. A --> A,A);\cr
#' (2) sympatric speciation where one species inherits the ancestral range, and the other inherits
#' a 1-area subset of the ancestral range (e.g. ABC --> ABC,B);\cr
#' (3) vicariant speciation with one daughter occupying an area of size 1 (e.g. ABCD --> ACD,B)\cr
#'
#' For example, if the ancestral range is ABC, the possible daughters are:
#' 
#' (Left, Right)
#' 
#' Vicariance:
#' A,BC
#' AB,C
#' AC,B
#' BC,A
#' C,AB
#' B,AC
#' 
#' Sympatric subset:
#' A,ABC
#' B,ABC
#' C,ABC
#' ABC,A
#' ABC,B
#' ABC,C
#' 
#' There are 12 possibilities, so LAGRANGE would give each a probability of 1/12, conditional on the
#' ancestor having range ABC.  All other imaginable scenarios are given probability 0 -- e.g., sympatric
#' speciation of a widespread range (ABC --> ABC,ABC), or jump dispersal leading to founder-event speciation
#' (ABC --> ABC,D).
#'
#' In \code{BioGeoBEARS}, the relative probability (or weight) of these categories is set by the \code{s} 
#' (sympatric-subset), \code{v} (vicariance), \code{j} (jump/founder-event), and \code{y} 
#' (sympatric-range-copying) parameters. These parameters do not have to sum to 1, they just give the 
#' \emph{relative} weight of an event of each type.  E.g., if \code{s=1}, \code{v=1}, \code{j=0}, \code{y=1}, 
#' then each allowed sympatric-range-copying, sympatric-subset, and vicariance event is given equal probability 
#' (this is the LAGRANGE cladogenesis model) .
#' 
#' The \code{rcpp_calc_anclikes_sp} function gets slow for large state spaces, as every possible combination of states at Left and
#' Right branches is checked.  Even in C++ this will get slow, as the \code{(number of states) = 2^(number of areas)},
#' and as the number of possible combinations of (ancestor, left,right) states is 
#' \code{(number of states)*(number of states)*(number of states)}.
#' 
#' Note: the \code{maxent} parameters allow the user to specify the probability distribution for different range sizes of the 
#' smaller-ranged descendant lineage.  The defaults set these parameters so that the LAGRANGE model is implemented (the smaller
#' descendant always has range size 1).
#'
#' See \code{\link{rcpp_calc_anclikes_sp_COOprobs}} and \code{\link{rcpp_calc_anclikes_sp_COOweights_faster}}
#' for successively faster solutions to this problem.
#' 
#' This is the byte-compiled version of \code{\link{rcpp_calc_anclikes_sp_prebyte}}.
#' \code{\link{rcpp_calc_anclikes_sp}} is byte-compiled, which (might) make it faster.  
#'
#' For information on byte-compiling, see \url{http://www.r-statistics.com/2012/04/speed-up-your-r-code-using-a-just-in-time-jit-compiler/} and \code{\link[compiler:cmpfun]{cmpfun}} in the \code{\link[compiler:cmpfun]{compiler}} package.
#'
#' @param Rcpp_leftprobs Probabilities of the states at the base of the left descendant branch
#' @param Rcpp_rightprobs Probabilities of the states at the base of the right descendant branch
#' @param l List of state indices (0-based)
#' @param s Relative weight of sympatric "subset" speciation. Default \code{s=1} mimics LAGRANGE model.
#' @param v Relative weight of vicariant speciation. Default \code{v=1} mimics LAGRANGE model.
#' @param j Relative weight of "founder event speciation"/jump speciation. Default \code{j=0} mimics LAGRANGE model.
#' @param y Relative weight of fully sympatric speciation (range-copying). Default \code{y=1} mimics LAGRANGE model.
#' @param dmat If given, a matrix of rank numareas giving multipliers for the probability
#' of each dispersal event between areas. Default NULL, which sets every cell of the 
#' \code{dmat} matrix to value 1.  Users may construct their own parameterized \code{dmat}
#' (for example, making \code{dmat} a function of distance) for inclusion in ML or
#' Bayesian analyses.
#' @param maxent01s Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a subset-sympatric speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent
#' always has range size 1).
#' @param maxent01v Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a vicariance speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param maxent01j Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a founder-event speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param maxent01y Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a full-sympatric (range-copying) 
#' speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param max_minsize_as_function_of_ancsize If given, any state with a range larger that this value will 
#' be given a probability of zero (for the branch with the smaller rangesize).  This means that not every
#' possible combination of ranges has to be checked, which can get very slow for large state spaces.
#' @param Rsp_rowsums A vector of size (numstates)  giving the sum of the relative probabilites of 
#' each combination of descendant states, assuming the probabilities of the left- and right-states are 
#' all equal (set to 1). This is thus the sum of the weights, and dividing by this normalization vector 
#' means that the each row of the speciation probability matrix will sum to 1.  Default assumes the 
#' weights sum to 1 but this is not usually the case. Rsp_rowsums need only be calculated once per 
#' tree+model combination, stored, and then re-used for each node in the tree, yielding significant 
#' time savings.
#' @param printmat Should the probability matrix output be printed to screen? (useful for debugging, but 
#' can be dramatically slow in R.app for some reason for even moderate numbers of states; perhaps 
#' overrunning the line length...)
#' @return \code{prob_ancestral_states} The probabilities of the ancestral states.
#' @export
#' @seealso \code{\link{rcpp_calc_anclikes_sp}}, \code{\link{rcpp_calc_anclikes_sp_COOprobs}}, 
#' \code{\link{rcpp_calc_anclikes_sp_COOweights_faster}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # For the basic logic of a probablistic cladogenesis model, see
#' ?rcpp_calc_anclikes_sp
#' 
#' # For examples of running the functions, see the comparison of all functions at:
#' # ?cladoRcpp
#'
rcpp_calc_anclikes_sp = compiler::cmpfun(rcpp_calc_anclikes_sp_prebyte)






#######################################################
# rcpp_calc_anclikes_sp_rowsums:  
#######################################################
#' Calculate the number of cladogenesis events of nonzero probability for each ancestral state
#'
#' This function takes the list of possible states (\emph{l}), and the parameters of a cladogenesis model
#' (\emph{s}, \emph{v}, \emph{j}, \emph{y}) (which are the relative weights of each of type of cladogenic range inheritance event)
#' and, for each ancestral state, sums the weights of allowed descendant events.  Dividing the weights in each row, by the sum of 
#' the weights for that row, provides the absolute probabilities of each transition, conditional on the ancestral state for that row.
#'
#' The inputs \code{Rcpp_leftprobs} and \code{Rcpp_rightprobs} are basically irrelevant here, but 
#' retained for symmetry with the other functions.  In effect, this function is identical with
#' \code{\link{rcpp_calc_anclikes_sp}} except that \code{Rcpp_leftprobs} and \code{Rcpp_rightprobs}
#' are arrays of 1s of \code{length(l)}, i.e. \code{length(number_of_states)}.
#'
#' This function is no longer used in \code{BioGeoBEARS}, but has been retained to enable easy counting of
#' the number of events.  When all nonzero-probability events are of equal probability (e.g. as in LAGRANGE; Ree & Smith 2008)
#' the function could be used for normalization, but it is safer to use \code{\link{rcpp_calc_anclikes_sp}} or 
#' one of the faster COO-like equivalents.
#'
#' @param Rcpp_leftprobs Probabilities of the states at the base of the left descendant branch
#' @param Rcpp_rightprobs Probabilities of the states at the base of the right descendant branch
#' @param l List of state indices (0-based)
#' @param s Relative weight of sympatric "subset" speciation. Default \code{s=1} mimics LAGRANGE model.
#' @param v Relative weight of vicariant speciation. Default \code{v=1} mimics LAGRANGE model.
#' @param j Relative weight of "founder event speciation"/jump speciation. Default \code{j=0} mimics LAGRANGE model.
#' @param y Relative weight of fully sympatric speciation (range-copying). Default \code{y=1} mimics LAGRANGE model.
#' @param dmat If given, a matrix of rank numareas giving multipliers for the probability
#' of each dispersal event between areas. Default NULL, which sets every cell of the 
#' \code{dmat} matrix to value 1.  Users may construct their own parameterized \code{dmat}
#' (for example, making \code{dmat} a function of distance) for inclusion in ML or
#' Bayesian analyses.
#' @param maxent01s Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a subset-sympatric speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent
#' always has range size 1).
#' @param maxent01v Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a vicariance speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param maxent01j Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a founder-event speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param maxent01y Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a full-sympatric (range-copying) 
#' speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param max_minsize_as_function_of_ancsize If given, any state with a range larger that this value will 
#' be given a probability of zero (for the branch with the smaller rangesize).  This means that not every
#' possible combination of ranges has to be checked, which can get very slow for large state spaces.
#' @return \code{Rsp_rowsums} A vector of size (numstates)  giving the number of events of nonzero probability 
#' for each ancestral states.
#' @param printmat Should the probability matrix output be printed to screen? (useful for debugging, but 
#' can be dramatically slow in R.app for some reason for even moderate numbers of states; perhaps 
#' overrunning the line length...)
#' @export
#' @seealso \code{\link{rcpp_calc_anclikes_sp}}, \code{\link{rcpp_calc_anclikes_sp_COOprobs}}, 
#' \code{\link{rcpp_calc_anclikes_sp_COOweights_faster}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#'   @cite ReeSmith2008
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # For the basic logic of a probablistic cladogenesis model, see
#' ?rcpp_calc_anclikes_sp
#' 
#' # For examples of running the functions, see the comparison of all functions at:
#' # ?cladoRcpp
#'
rcpp_calc_anclikes_sp_rowsums <- function(Rcpp_leftprobs, Rcpp_rightprobs, l, s=1, v=1, j=0, y=1, dmat=NULL, maxent01s=NULL, maxent01v=NULL, maxent01j=NULL, maxent01y=NULL, max_minsize_as_function_of_ancsize=NULL, printmat=TRUE)
	{
	#print(Rcpp_leftprobs)
	#print(Rcpp_rightprobs)
	
	# Take the max of the indices of the possible areas, and add 1
	numareas = max(unlist(l), na.rm=TRUE) + 1

	# max_minsize_as_function_of_ancsize can just be all 1s as default
	if (is.null(max_minsize_as_function_of_ancsize))
		{
		max_minsize_as_function_of_ancsize = rep(1, numareas)
		}
	
	# If needed, fill in dmat
	if (is.null(dmat))
		{
		dmat = matrix(data=1, nrow=numareas, ncol=numareas)
		}
	
	##################################################################
	# If needed, put e.g. 
	# 1 0 0
	# 1 0 0
	# 1 0 0
	# ...into the maxent distributions
	##################################################################
	
	# subset sympatric speciation
	if (is.null(maxent01s))
		{
		maxent01s = matrix(0, nrow=numareas, ncol=numareas)
		maxent01s[,1] = 1
		}
	
	# vicariant speciation
	# the max range size of the smaller range is 
	# numareas/2, rounded down
	if (is.null(maxent01v))
		{
		numareas_smaller_vicariant_range = floor(numareas / 2)
		#maxent01v = matrix(0, nrow=numareas, ncol=numareas_smaller_vicariant_range)
		maxent01v = matrix(0, nrow=numareas, ncol=numareas)
		maxent01v[,1] = 1
		#maxent01v = relative_probabilities_of_vicariants(max_numareas=numareas, maxent_constraint_01v=maxent01v_param, NA_val=-1)
		}
	
	# founder-event speciation ("jump dispersal")
	if (is.null(maxent01j))
		{
		maxent01j = matrix(0, nrow=numareas, ncol=numareas)
		maxent01j[,1] = 1
		}
	
	# full sympatric speciation (range copying)
	# This is essentially subset speciation where the "smaller" range is the
	# same size as the "bigger" range
	if (is.null(maxent01y))
		{
		maxent01y = matrix(0, nrow=numareas, ncol=numareas)
		maxent01y[,1] = 1
		}
	
	
	
#	tmp = unlist(maxent01s[(maxent01s > 0 && !is.na(maxent01s))])
# 	max_min_rangesize_s = 1:length(maxent01s[maxent01s > 0])
# 	max_min_rangesize_v = 1:length(maxent01v[maxent01v > 0])
# 	max_min_rangesize_j = 1:length(maxent01j[maxent01j > 0])
# 	max_min_rangesize_y = 1:length(maxent01y[maxent01y > 0])
# 	
# 	max_min_rangesize = max(unlist(c(max_min_rangesize_s, max_min_rangesize_v, max_min_rangesize_j)))

	# Don't use a, b, etc directly, they get screwed up
	tmpa = Rcpp_leftprobs
	tmpb = Rcpp_rightprobs

	# Print the matrix output to screen?	
	if (printmat == TRUE)
		{
		Rprintmat = 1
		} else {
		Rprintmat = 0
		}

	
	# Call the fast C++ function
	sp_rowsums = .Call( "cpp_calc_anclikes_sp_rowsums", Rprintmat=as.integer(Rprintmat), cpp_leftprobs=as.numeric(tmpa), cpp_leftprobs=as.numeric(tmpb), l=l, s=as.numeric(s), v=as.numeric(v), j=as.numeric(j), y=as.numeric(y), dmat=as.matrix(dmat), maxent01s=as.matrix(maxent01s), maxent01v=as.matrix(maxent01v), maxent01j=as.matrix(maxent01j), maxent01y=as.matrix(maxent01y), max_minsize_as_function_of_ancsize=as.integer(max_minsize_as_function_of_ancsize), PACKAGE = "cladoRcpp" )

	#print(Rcpp_leftprobs)
	#print(Rcpp_rightprobs)
	
	return(sp_rowsums)
	}




#######################################################
# rcpp_calc_anclikes_sp_COOprobs:
#######################################################
#' Faster version of rcpp_calc_anclikes_sp
#'
#' This function is a faster version of \code{\link{rcpp_calc_anclikes_sp}}. Like
#' \code{\link{rcpp_calc_anclikes_sp}}, this function calculates the conditional
#' probability of every allowed combination of ancestral range, left descendent range,
#' and right descendent range.
#' 
#' This function improves upon \code{\link{rcpp_calc_anclikes_sp}} by
#' returning a COO-like list of the nonzero cells in the transition matrix
#' for the speciation event.
#' 
#' (COO = Coordinate list format for a matrix, see 
#' \url{http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29}
#' 
#' Whereas a COO-formatted square matrix stores, for each nonzero cell, the row #, column #, and 
#' cell value, \code{\link{rcpp_calc_anclikes_sp}} returns lists containing, for each nonzero cell:
#'
#' 1. 0-based index of the ancestral state\cr
#' 2. 0-based index of the left state\cr
#' 3. 0-based index of the right state\cr
#' 4. Value of the specified nonzero cell\cr
#'
#' Time savings over \code{\link{rcpp_calc_anclikes_sp}} are realized by skipping many 
#' ancestor/descendent combinations which are impossible transitions on the model, and 
#' neither recording, nor storing, nor passing them.  This becomes important with 
#' large state spaces.
#'
#' @param Rcpp_leftprobs Probabilities of the states at the base of the left descendant branch
#' @param Rcpp_rightprobs Probabilities of the states at the base of the right descendant branch
#' @param l List of state indices (0-based)
#' @param s Relative weight of sympatric "subset" speciation. Default \code{s=1} mimics LAGRANGE model.
#' @param v Relative weight of vicariant speciation. Default \code{v=1} mimics LAGRANGE model.
#' @param j Relative weight of "founder event speciation"/jump speciation. Default \code{j=0} mimics LAGRANGE model.
#' @param y Relative weight of fully sympatric speciation (range-copying). Default \code{y=1} mimics LAGRANGE model.
#' @param dmat If given, a matrix of rank numareas giving multipliers for the probability
#' of each dispersal event between areas. Default NULL, which sets every cell of the 
#' \code{dmat} matrix to value 1.  Users may construct their own parameterized \code{dmat}
#' (for example, making \code{dmat} a function of distance) for inclusion in ML or
#' Bayesian analyses.
#' @param maxent01s Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a subset-sympatric speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent
#' always has range size 1).
#' @param maxent01v Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a vicariance speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param maxent01j Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a founder-event speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param maxent01y Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a full-sympatric (range-copying) 
#' speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param max_minsize_as_function_of_ancsize If given, any state with a range larger that this value will 
#' be given a probability of zero (for the branch with the smaller rangesize).  This means that not every
#' possible combination of ranges has to be checked, which can get very slow for large state spaces.
#' @param printmat Should the probability matrix output be printed to screen? (useful for debugging, but 
#' can be dramatically slow in R.app for some reason for even moderate numbers of states; perhaps 
#' overrunning the line length...)
#' @return \code{list_weights_of_transitions} A list of 3 lists. Each list has (numstates) items, 
#' representing the ancestral states.  List #1 gives the 0-based state index for the nonzero left descendents
#' of each ancestral state. List #2 gives the 0-based state index for the nonzero right descendents
#' of each ancestral state. List #3 gives the weight of each nonzero transition from each ancestral state.
#' Summing these weights within each ancestral state for list #3 gives the total of the weights for
#' each ancestral state.  Dividing the weights by the sum of weights gives the conditional probability
#' of each descendent state, conditional on the ancestral state.  These conditional probabilities
#' need only be calculated once per 
#' tree+model combination, stored, and then re-used for each node in the tree, yielding significant 
#' time savings.
#' @export
#' @seealso \code{\link{rcpp_calc_anclikes_sp}}, \code{\link{rcpp_calc_anclikes_sp_COOprobs}}, 
#' \code{\link{rcpp_calc_anclikes_sp_COOweights_faster}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#'	 @cite ReeSmith2008
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # For the basic logic of a probablistic cladogenesis model, see
#' ?rcpp_calc_anclikes_sp
#' 
#' # For examples of running the functions, see the comparison of all functions at:
#' # ?cladoRcpp
#'
rcpp_calc_anclikes_sp_COOprobs <- function(Rcpp_leftprobs, Rcpp_rightprobs, l, s=1, v=1, j=0, y=1, dmat=NULL, maxent01s=NULL, maxent01v=NULL, maxent01j=NULL, maxent01y=NULL, max_minsize_as_function_of_ancsize=NULL, printmat=TRUE)
	{
	#print(Rcpp_leftprobs)
	#print(Rcpp_rightprobs)
	
	# Take the max of the indices of the possible areas, and add 1
	numareas = max(unlist(l), na.rm=TRUE) + 1

	# max_minsize_as_function_of_ancsize can just be all 1s as default
	if (is.null(max_minsize_as_function_of_ancsize))
		{
		max_minsize_as_function_of_ancsize = rep(1, numareas)
		}
	
	# If needed, fill in dmat
	if (is.null(dmat))
		{
		dmat = matrix(data=1, nrow=numareas, ncol=numareas)
		}
	
	##################################################################
	# If needed, put e.g. 
	# 1 0 0
	# 1 0 0
	# 1 0 0
	# ...into the maxent distributions
	##################################################################
	
	# subset sympatric speciation
	if (is.null(maxent01s))
		{
		maxent01s = matrix(0, nrow=numareas, ncol=numareas)
		maxent01s[,1] = 1
		}
	
	# vicariant speciation
	# the max range size of the smaller range is 
	# numareas/2, rounded down
	if (is.null(maxent01v))
		{
		numareas_smaller_vicariant_range = floor(numareas / 2)
		#maxent01v = matrix(0, nrow=numareas, ncol=numareas_smaller_vicariant_range)
		maxent01v = matrix(0, nrow=numareas, ncol=numareas)
		maxent01v[,1] = 1
		#maxent01v = relative_probabilities_of_vicariants(max_numareas=numareas, maxent_constraint_01v=maxent01v_param, NA_val=-1)
		}
	
	# founder-event speciation ("jump dispersal")
	if (is.null(maxent01j))
		{
		maxent01j = matrix(0, nrow=numareas, ncol=numareas)
		maxent01j[,1] = 1
		}
	
	# full sympatric speciation (range copying)
	# This is essentially subset speciation where the "smaller" range is the
	# same size as the "bigger" range
	if (is.null(maxent01y))
		{
		maxent01y = matrix(0, nrow=numareas, ncol=numareas)
		maxent01y[,1] = 1
		}
	
	
	
#	tmp = unlist(maxent01s[(maxent01s > 0 && !is.na(maxent01s))])
# 	max_min_rangesize_s = 1:length(maxent01s[maxent01s > 0])
# 	max_min_rangesize_v = 1:length(maxent01v[maxent01v > 0])
# 	max_min_rangesize_j = 1:length(maxent01j[maxent01j > 0])
# 	max_min_rangesize_y = 1:length(maxent01y[maxent01y > 0])
# 	
# 	max_min_rangesize = max(unlist(c(max_min_rangesize_s, max_min_rangesize_v, max_min_rangesize_j)))

	# Don't use a, b, etc directly, they get screwed up
	tmpa = Rcpp_leftprobs
	tmpb = Rcpp_rightprobs

	# Print the matrix output to screen?	
	if (printmat == TRUE)
		{
		Rprintmat = 1
		} else {
		Rprintmat = 0
		}

	
	# Call the fast C++ function
	list_weights_of_transitions = .Call( "cpp_calc_anclikes_sp_COOprobs", Rprintmat=as.integer(Rprintmat), cpp_leftprobs=as.numeric(tmpa), cpp_leftprobs=as.numeric(tmpb), l=l, s=as.numeric(s), v=as.numeric(v), j=as.numeric(j), y=as.numeric(y), dmat=as.matrix(dmat), maxent01s=as.matrix(maxent01s), maxent01v=as.matrix(maxent01v), maxent01j=as.matrix(maxent01j), maxent01y=as.matrix(maxent01y), max_minsize_as_function_of_ancsize=as.integer(max_minsize_as_function_of_ancsize), PACKAGE = "cladoRcpp" )

	#print(Rcpp_leftprobs)
	#print(Rcpp_rightprobs)
	
	return(list_weights_of_transitions)
	}







#######################################################
# rcpp_calc_anclikes_sp_COOweights_faster
#######################################################
#' Even faster version of rcpp_calc_anclikes_sp
#'
#' This function improves on \code{\link{rcpp_calc_anclikes_sp}} and
#' \code{\link{rcpp_calc_anclikes_sp_COOprobs}}.  In addition to the compressed
#' COO-like storage format, the internal C++ code here explicitly
#' enumerates the allowed transitions, rather than searching through
#' every possibility and testing whether or not it is allowed.  This 
#' appears to scale well to very large state spaces.
#'
#' This should be faster, i.e. by look for each type of event individually.
#' 
#' Returns results as 4 columns: ancestral index, left index, right index, conditional
#' probability given ancestral states (assuming likelihood of descendants is 1). Indexes
#' are 0-based.
#' 
#' Keep in mind that cladogenesis matrices exclude the null state
#' (a range of 0 areas), so if your states list starts with the 
#' null range (as is typical/default in DEC-style models)
#' then to get the R 1-based state indices requires e.g. 
#' COO_weights_columnar[[1]] + 2.
#'
#' When the calculation is run at each node in the tree, all that is required is one
#' pass through the COO-like array, with the downpassed probabilities of the
#' states on the left and right branches multiplied by the probability column.
#'
#' @param Rcpp_leftprobs Probabilities of the states at the base of the left descendant branch
#' @param Rcpp_rightprobs Probabilities of the states at the base of the right descendant branch
#' @param l List of state indices (0-based)
#' @param s Relative weight of sympatric "subset" speciation. Default \code{s=1} mimics LAGRANGE model.
#' @param v Relative weight of vicariant speciation. Default \code{v=1} mimics LAGRANGE model.
#' @param j Relative weight of "founder event speciation"/jump speciation. Default \code{j=0} mimics LAGRANGE model.
#' @param y Relative weight of fully sympatric speciation (range-copying). Default \code{y=1} mimics LAGRANGE model.
#' @param dmat If given, a matrix of rank numareas giving multipliers for the probability
#' of each dispersal event between areas. Default NULL, which sets every cell of the 
#' \code{dmat} matrix to value 1.  Users may construct their own parameterized \code{dmat}
#' (for example, making \code{dmat} a function of distance) for inclusion in ML or
#' Bayesian analyses.
#' @param maxent01s Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a subset-sympatric speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent
#' always has range size 1).
#' @param maxent01v Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a vicariance speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param maxent01j Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a founder-event speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param maxent01y Matrix giving the relative weight of each possible descendant rangesize for 
#' the smaller range, for a given ancestral rangesize, for a full-sympatric (range-copying) 
#' speciation event. 
#' Default is \code{NULL}, which means the script will set up the LAGRANGE model (one descendent 
#' always has range size 1).
#' @param max_minsize_as_function_of_ancsize If given, any state with a range larger that this value will 
#' be given a probability of zero (for the branch with the smaller rangesize).  This means that not every
#' possible combination of ranges has to be checked, which can get very slow for large state spaces.
#' @param printmat Should the probability matrix output be printed to screen? (useful for debugging, but 
#' can be dramatically slow in R.app for some reason for even moderate numbers of states; perhaps 
#' overrunning the line length...)
#' @return \code{COO_weights_columnar} Transition weights matrix in COO-like format as 4 columns: 
#' ancestral index, left index, right index, and weight of the specified transition. Indexes are
#' 0-based. 
#' Keep in mind that cladogenesis matrices exclude the null state
#' (a range of 0 areas), so if your states list starts with the 
#' null range (as is typical/default in DEC-style models)
#' then to get the R 1-based state indices requires e.g. 
#' COO_weights_columnar[[1]] + 2.
#' 
#' Dividing the
#' weights by the sum of the weights for a particular ancestral state yields the conditional
#' probabilities of each transition at the speciation event.
#' (assuming likelihood of descendants is 1).
#' @export
#' @seealso \code{\link{rcpp_calc_anclikes_sp}}, \code{\link{rcpp_calc_anclikes_sp_COOprobs}}, 
#' \code{\link{rcpp_calc_anclikes_sp_COOweights_faster}}
#' @seealso \code{\link{rcpp_calc_anclikes_sp}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#'   @cite ReeSmith2008
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # For the basic logic of a probablistic cladogenesis model, see
#' ?rcpp_calc_anclikes_sp
#' 
#' # For examples of running the functions, see the comparison of all functions at:
#' # ?cladoRcpp
#'
rcpp_calc_anclikes_sp_COOweights_faster <- function(Rcpp_leftprobs, Rcpp_rightprobs, l, s=1, v=1, j=0, y=1, dmat=NULL, maxent01s=NULL, maxent01v=NULL, maxent01j=NULL, maxent01y=NULL, max_minsize_as_function_of_ancsize=NULL, printmat=TRUE)
	{
	#print(Rcpp_leftprobs)
	#print(Rcpp_rightprobs)
	
	# Take the max of the indices of the possible areas, and add 1
	numareas = max(unlist(l), na.rm=TRUE) + 1

	# max_minsize_as_function_of_ancsize can just be all 1s as default
	if (is.null(max_minsize_as_function_of_ancsize))
		{
		max_minsize_as_function_of_ancsize = rep(1, numareas)
		}
	
	# If needed, fill in dmat
	if (is.null(dmat))
		{
		dmat = matrix(data=1, nrow=numareas, ncol=numareas)
		}
	
	##################################################################
	# If needed, put e.g. 
	# 1 0 0
	# 1 0 0
	# 1 0 0
	# ...into the maxent distributions
	##################################################################
	
	# subset sympatric speciation
	if (is.null(maxent01s))
		{
		maxent01s = matrix(0, nrow=numareas, ncol=numareas)
		maxent01s[,1] = 1
		}
	
	# vicariant speciation
	# the max range size of the smaller range is 
	# numareas/2, rounded down
	if (is.null(maxent01v))
		{
		numareas_smaller_vicariant_range = floor(numareas / 2)
		#maxent01v = matrix(0, nrow=numareas, ncol=numareas_smaller_vicariant_range)
		maxent01v = matrix(0, nrow=numareas, ncol=numareas)
		maxent01v[,1] = 1
		#maxent01v = relative_probabilities_of_vicariants(max_numareas=numareas, maxent_constraint_01v=maxent01v_param, NA_val=-1)
		}
	
	# founder-event speciation ("jump dispersal")
	if (is.null(maxent01j))
		{
		maxent01j = matrix(0, nrow=numareas, ncol=numareas)
		maxent01j[,1] = 1
		}
	
	# full sympatric speciation (range copying)
	# This is essentially subset speciation where the "smaller" range is the
	# same size as the "bigger" range
	if (is.null(maxent01y))
		{
		maxent01y = matrix(0, nrow=numareas, ncol=numareas)
		maxent01y[,1] = 1
		}
	
	
	
#	tmp = unlist(maxent01s[(maxent01s > 0 && !is.na(maxent01s))])
# 	max_min_rangesize_s = 1:length(maxent01s[maxent01s > 0])
# 	max_min_rangesize_v = 1:length(maxent01v[maxent01v > 0])
# 	max_min_rangesize_j = 1:length(maxent01j[maxent01j > 0])
# 	max_min_rangesize_y = 1:length(maxent01y[maxent01y > 0])
# 	
# 	max_min_rangesize = max(unlist(c(max_min_rangesize_s, max_min_rangesize_v, max_min_rangesize_j)))

	# Don't use a, b, etc directly, they get screwed up
	tmpa = Rcpp_leftprobs
	tmpb = Rcpp_rightprobs

	# Print the matrix output to screen?	
	if (printmat == TRUE)
		{
		Rprintmat = 1
		} else {
		Rprintmat = 0
		}

	
	# Call the fast C++ function
	COO_weights_columnar = .Call( "cpp_calc_anclikes_sp_COOweights_faster", Rprintmat=as.integer(Rprintmat), cpp_leftprobs=as.numeric(tmpa), cpp_leftprobs=as.numeric(tmpb), l=l, s=as.numeric(s), v=as.numeric(v), j=as.numeric(j), y=as.numeric(y), dmat=as.matrix(dmat), maxent01s=as.matrix(maxent01s), maxent01v=as.matrix(maxent01v), maxent01j=as.matrix(maxent01j), maxent01y=as.matrix(maxent01y), max_minsize_as_function_of_ancsize=as.integer(max_minsize_as_function_of_ancsize), PACKAGE = "cladoRcpp" )

	#print(Rcpp_leftprobs)
	#print(Rcpp_rightprobs)
	
	return(COO_weights_columnar)
	}










#######################################################
# rcpp_calc_anclikes_sp_using_COOprobs
#######################################################
#' Calculate ancestral likelihoods given a COO-like probability matrix
#'
#' This function does a pass through a COO-like transition probability matrix 
#' for a node, inputting the probabilities that have been passed down
#' from above for the left and right branch, and the sum of weights for
#' each ancestral state, and returns the ancestral relative probabilities.
#'
#' This C++ implementation should be slightly faster than the R version,
#' although for a simple pass through an array the difference may not
#' be great.
#'
#' @param Rcpp_leftprobs Probabilities of the states at the base of the left descendant branch
#' @param Rcpp_rightprobs Probabilities of the states at the base of the right descendant branch
#' @param RCOO_left_i_list 0-based index of the allowed left states\cr
#' @param RCOO_right_j_list 0-based index of the allowed right states\cr
#' @param RCOO_probs_list Value of the specified nonzero cells\cr
#' @param Rsp_rowsums A vector of size (numstates)  giving the sum of the relative probabilites of 
#' each combination of descendant states, assuming the probabilities of the left- and right-states are 
#' all equal (set to 1). This is thus the sum of the weights, and dividing by this normalization vector 
#' means that the each row of the speciation probability matrix will sum to 1.  Default assumes the 
#' weights sum to 1 but this is not usually the case. Rsp_rowsums need only be calculated once per 
#' tree+model combination, stored, and then re-used for each node in the tree, yielding significant 
#' time savings.
#' @param printmat Should the probability matrix output be printed to screen? (useful for debugging, but 
#' can be dramatically slow in R.app for some reason for even moderate numbers of states; perhaps 
#' overrunning the line length...)
#' @return \code{R_anc_relprobs} Vector of the probabilities of the ancestral states
#' @export
#' @seealso \code{\link{rcpp_calc_anclikes_sp}}
#' @seealso \code{\link{rcpp_calc_anclikes_sp}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # For the basic logic of a probablistic cladogenesis model, see
#' ?rcpp_calc_anclikes_sp
#' 
#' # For examples of running the functions, see the comparison of all functions at:
#' # ?cladoRcpp
#'
rcpp_calc_anclikes_sp_using_COOprobs <- function(Rcpp_leftprobs, Rcpp_rightprobs, RCOO_left_i_list, RCOO_right_j_list, RCOO_probs_list, Rsp_rowsums, printmat=TRUE)
	{
	# Don't use a, b, etc directly, they get screwed up
	tmpa = Rcpp_leftprobs
	tmpb = Rcpp_rightprobs

	# Print the matrix output to screen?	
	if (printmat == TRUE)
		{
		Rprintmat = 1
		} else {
		Rprintmat = 0
		}

	
	# Call the fast C++ function
	R_anc_relprobs = .Call( "cpp_calc_anclikes_sp_using_COOprobs", Rprintmat=as.integer(Rprintmat), leftprobs=as.numeric(tmpa), rightprobs=as.numeric(tmpb), RCOO_left_i_list=as.list(RCOO_left_i_list), RCOO_right_j_list=as.list(RCOO_right_j_list), RCOO_probs_list=as.list(RCOO_probs_list), Rsp_rowsums=as.numeric(Rsp_rowsums), PACKAGE = "cladoRcpp" )

	#print(Rcpp_leftprobs)
	#print(Rcpp_rightprobs)
	
	return(R_anc_relprobs)
	}









#######################################################
# rcpp_calc_rowsums_for_COOweights_columnar
#######################################################
#' Calculate sum of weights for each ancestral state
#'
#' This is a C++ implementation of \code{\link{rcpp_calc_anclikes_sp_rowsums}}.  It should
#' be substantially faster, as it requires only one pass through \code{COO_weights_columnar}.
#'
#'
#' @param COO_weights_columnar Transition probability matrix in COO-like format as 4 columns: 
#' ancestral index, left index, right index, conditional probability given ancestral states.
#' (assuming likelihood of descendants is 1). Indexes are 0-based.
#' Keep in mind that cladogenesis matrices exclude the null state
#' (a range of 0 areas), so if your states list starts with the 
#' null range (as is typical/default in DEC-style models)
#' then to get the R 1-based state indices requires e.g. 
#' COO_weights_columnar[[1]] + 2.
#' @param numstates The user should provide the number of states (WITHOUT counting the null range),
#' in case they are not all
#' present in \code{COO_weights_columnar}.  If empty, the function assumes that the highest 
#' index represents the last state, and adds 1 to get the number of states. This may be a 
#' hazardous assumption.
#' @param printmat Should the probability matrix output be printed to screen? (useful for debugging, but 
#' can be dramatically slow in R.app for some reason for even moderate numbers of states; perhaps 
#' overrunning the line length...)
#' @return \code{rowsums} A vector of size (numstates)  giving the sum of the relative probabilites of 
#' each combination of descendant states, assuming the probabilities of the left- and right-states are 
#' all equal (set to 1). This is thus the sum of the weights, and dividing by this normalization vector 
#' means that the each row of the speciation probability matrix will sum to 1.  Default assumes the 
#' weights sum to 1 but this is not usually the case. Rsp_rowsums need only be calculated once per 
#' tree+model combination, stored, and then re-used for each node in the tree, yielding significant 
#' time savings.
#' @export
#' @seealso \code{\link{rcpp_calc_anclikes_sp}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # For the basic logic of a probablistic cladogenesis model, see
#' ?rcpp_calc_anclikes_sp
#' 
#' # For examples of running the functions, see the comparison of all functions at:
#' # ?cladoRcpp
#'
rcpp_calc_rowsums_for_COOweights_columnar <- function(COO_weights_columnar, numstates=1+max(sapply(X=COO_weights_columnar, FUN=max)[1:3]), printmat=TRUE)
	{
	# Note: the default numstates max is NOT ideal
	

	# Print the matrix output to screen?	
	if (printmat == TRUE)
		{
		Rprintmat = 1
		} else {
		Rprintmat = 0
		}

	RCOO_weights_columnar_anc_i_list = COO_weights_columnar[[1]]
	RCOO_probs_list = COO_weights_columnar[[4]]
	
	# Call the fast C++ function
	rowsums = .Call( "cpp_calc_rowsums_for_COOweights_columnar", RCOO_weights_columnar_anc_i_list=as.integer(RCOO_weights_columnar_anc_i_list), RCOO_probs_list=as.numeric(RCOO_probs_list), Rnumstates=as.integer(numstates), PACKAGE = "cladoRcpp" )

	return(rowsums)
	}



#######################################################
# rcpp_calc_splitlikes_using_COOweights_columnar
#######################################################
#' Calculate the split likelihoods using \code{COO_weights_columnar}
#'
#' Calculates the split likelihoods using \code{COO_weights_columnar}, i.e. the weights as produced by
#' \code{\link{rcpp_calc_anclikes_sp_COOweights_faster}}.
#'
#' @param Rcpp_leftprobs Probabilities of the states at the base of the left descendant branch
#' @param Rcpp_rightprobs Probabilities of the states at the base of the right descendant branch
#' @param COO_weights_columnar Transition probability matrix in COO-like format as 4 columns: 
#' ancestral index, left index, right index, conditional probability given ancestral states.
#' (assuming likelihood of descendants is 1). Indexes are 0-based.
#' Keep in mind that cladogenesis matrices exclude the null state
#' (a range of 0 areas), so if your states list starts with the 
#' null range (as is typical/default in DEC-style models)
#' then to get the R 1-based state indices requires e.g. 
#' COO_weights_columnar[[1]] + 2.
#' @param Rsp_rowsums A vector of size (numstates)  giving the sum of the relative probabilites of 
#' each combination of descendant states, assuming the probabilities of the left- and right-states are 
#' all equal (set to 1). This is thus the sum of the weights, and dividing by this normalization vector 
#' means that the each row of the speciation probability matrix will sum to 1.  Default assumes the 
#' weights sum to 1 but this is not usually the case. Rsp_rowsums need only be calculated once per 
#' tree+model combination, stored, and then re-used for each node in the tree, yielding significant 
#' time savings.
#' @param printmat Should the probability matrix output be printed to screen? (useful for debugging, but 
#' can be dramatically slow in R.app for some reason for even moderate numbers of states; perhaps 
#' overrunning the line length...)
#' @return \code{splitlikes} Vector of the probabilities of each allowed split
#' @export
#' @seealso \code{\link{rcpp_calc_anclikes_sp}}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'   @cite Matzke_2012_IBS
#' @author Nicholas Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # For the basic logic of a probablistic cladogenesis model, see
#' ?rcpp_calc_anclikes_sp
#' 
#' # For examples of running the functions, see the comparison of all functions at:
#' # ?cladoRcpp
#' 
rcpp_calc_splitlikes_using_COOweights_columnar <- function(Rcpp_leftprobs, Rcpp_rightprobs, COO_weights_columnar, Rsp_rowsums, printmat=TRUE)
	{
	# Note: the default numstates max is NOT ideal
	

	# Print the matrix output to screen?	
	if (printmat == TRUE)
		{
		Rprintmat = 1
		} else {
		Rprintmat = 0
		}

	RCOO_weights_columnar_anc_i_list = COO_weights_columnar[[1]]
	RCOO_left_i_list = COO_weights_columnar[[2]]
	RCOO_right_j_list = COO_weights_columnar[[3]]
	RCOO_probs_list = COO_weights_columnar[[4]]
	
	# Call the fast C++ function
	splitlikes = .Call( "cpp_calc_splitlikes_using_COOweights_columnar", leftprobs=as.numeric(Rcpp_leftprobs), rightprobs=as.numeric(Rcpp_rightprobs), RCOO_weights_columnar_anc_i_list=as.integer(RCOO_weights_columnar_anc_i_list), RCOO_left_i_list=as.integer(RCOO_left_i_list), RCOO_right_j_list=as.integer(RCOO_right_j_list), RCOO_probs_list=as.numeric(RCOO_probs_list), Rsp_rowsums=as.numeric(Rsp_rowsums), PACKAGE = "cladoRcpp" )

	return(splitlikes)
	}







