#' Cladogenic probability calculations using Rcpp
#'
#' \tabular{ll}{
#' Package: \tab cladoRcpp\cr
#' Type: \tab Package\cr
#' Version: \tab 0.14.4\cr
#' Date: \tab 2014-05-17\cr
#' License: \tab GPL (>= 3)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' \bold{Summary:} This package implements in C++/Rcpp various cladogenesis-related calculations that are slow
#' in pure R. These include the calculation of the probability of various scenarios for the inheritance
#' of geographic range at the divergence events on a phylogenetic tree, and other calculations
#' necessary for models which are not 
#' continuous-time markov chains (CTMC), but where change instead occurs instantaneously at speciation 
#' events.  Typically these models must assess the probability of every possible combination of (ancestor 
#' state, left descendent state, right descendent state).  This means that there are up to (# of states)^3 
#' combinations to investigate, and in biogeographical models, there can easily be hundreds of states, so 
#' calculation time becomes an issue.  C++ implementation plus clever tricks (many combinations can be 
#' eliminated a priori) can greatly speed the computation time over naive R implementations.
#' 
#' CITATION INFO: This package is the result of my Ph.D. research, please cite the package if you use it!  
#' Type: \code{citation(package="cladoRcpp")} to get the citation information.
#' 
#' See also the citation information for the sister packages,
#' \code{citation(package="rexpokit")} and \code{citation(package="BioGeoBEARS")}
#' 
#' 
#' \bold{Further information:} In particular, in \code{cladoRcpp}, functions are implemented to calculate the probability, given a model, 
#' of various scenarios for the inheritance of geographic range at speciation
#' events, where the left and right branches may inherit ranges different from each other
#' and different from the ancestor.
#'
#' The documentation for \code{\link{rcpp_areas_list_to_states_list}}, and \code{\link{rcpp_calc_anclikes_sp}} contain the 
#' basic introduction to the logic of ancestral states and cladogenesis probabilities with historical-biogeography models.
#'
#' The widely-used historical biogeography program \code{LAGRANGE} (Ree & Smith 2008) 
#' has only one cladogenesis model, which is fixed and therefore not subject to inference. LAGRANGE's cladogenesis
#' model gives equal weight/equal probability to all allowed cladogenesis events.  LAGRANGE allows:
#'
#' \bold{1.} sympatric speciation (copying the ancestral range to descendant ranges), but only for ranges of size=1 area\cr
#' \bold{2.} vicariant speciation (the descendant range is divided between the 2 descendant species), but at least one of the descendants must have a ranges of size=1 area\cr
#' \bold{3.} sympatric "subset" speciation (one species starts inside the ancestral range, the other inherits the ancestral range); again, one of the descendants must have a ranges of size=1 area\cr
#'
#' But, another range inheritance scenario is imaginable: 
#'
#' \bold{4.} founder-event speciation, where one descendant species inherits the ancestral range, and the other species has a range completely outside of the ancestral range
#'
#' \code{cladoRcpp} allows specification of these different models, including allowing different weights for the different processes, if users would like to infer the
#' optimal model, rather than simply fixing it ahead of time.  The optimization and model choice is done with the help of the sister packages, \code{rexpokit} and \code{BioGeoBEARS}.
#' 
#' \emph{Note:} I began this package with a little bit of code from Rcpp and the various
#' examples that have been written with it, as well as from the following:
#' 
#' 1. phyloRcppExamples by Vladimir Minin 
#' (https://r-forge.r-project.org/scm/viewvc.php/pkg/phyloRcppExamples/?root=evolmod and
#' http://markovjumps.blogspot.com/2012/01/packaging-and-exposing-rcpp-functions.html) --
#' which shows how to do phylogenetic operations in C++, accessed with R
#' 
#' 2. rcppbugs by Whit Armstrong (https://github.com/armstrtw/rcppbugs and 
#' http://cran.r-project.org/web/packages/rcppbugs/index.html) -- which does BUGS-style 
#' MCMC via C++ functions wrapped in R.  It does this much faster than MCMCpack and
#' rjags.
#' 
#'
#' @name cladoRcpp-package
#' @aliases cladoRcpp
#' @docType package
#' @title Phylogenetic probability calculations using Rcpp
#' @note Some starting code borrowed from Rcpp examples, Whit Armstrong's rcppbugs, and Vladimir Minin's phyloRcppExamples.
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/biogeobears}
#' @bibliography /Dropbox/_njm/__packages/cladoRcpp_setup/cladoRcpp_refs.bib
#'  @cite Matzke_2012_IBS
#'  @cite Ree_etal_2005
#'  @cite ReeSmith2008
#' @keywords package Rcpp phyloRcppExamples rcppbugs RcppArmadillo
#' @seealso \code{\link{rcpp_calc_anclikes_sp}}, \code{\link{rcpp_areas_list_to_states_list}}, Rcpp, RcppArmadillo
#' @examples
#' 
#' # To get citation information for cladoRcpp, type:
#' citation(package="cladoRcpp")
#' 
#' # Please also cite the sister packages I created to utilize rexpokit:
#' # citation(package="rexpokit")		# Roger Sidje is a coauthor 
#'                                  	# of rexpokit and author of 
#'                                  	# the FORTRAN EXPOKIT
#' # citation(package="BioGeoBEARS")
#' 
#' library(cladoRcpp)
#' # Test this first as it causes problems for --gct or --use-valgrind
#' areas_list = c("A", "B", "C")
#' areas_list
#' 
#' # Calculate the list of 0-based indices for each possible 
#' # geographic range, i.e. each combination of areas
#' \dontrun{
#' 
#' states_list = rcpp_areas_list_to_states_list(areas=areas_list, maxareas=3, 
#' include_null_range=FALSE)
#' 
#' }
#' 
#' #################################################################################
#' # Examples using C++ to speed up the slow step of getting all possible combinations
#' # of areas (important when when number_of_areas >= 7, as this can mean 
#' # 2^number_of_areas states, and (2^number_of_areas)^2 imaginable descendant pairs
#' # from each ancestral state.
#' #################################################################################
#' 
#' #######################################################
#' # Set up 2 vectors, then convolve them
#' #######################################################
#' ca = c(1,2,3,4,5)
#' cb = c(2,2,2,2,2)
#' rcpp_convolve(a=ca, b=cb)
#' #'
#' # Same as:
#' convolve(ca, cb, conj=TRUE, type="open")
#' 
#' 
#' 
#' #######################################################
#' # Cross-products
#' #######################################################
#' ca = c(1,2,3,4,5)
#' cb = c(2,2,2,2,2)
#' rcpp_mult2probvect(a=ca, b=cb)
#' 
#' # Same as:
#' c(ca %o% cb)
#' 
#' # Or:
#' outer(ca, cb)
#' c(outer(ca, cb))
#' 
#' # Or:
#' tcrossprod(ca, cb)
#' c(tcrossprod(ca, cb))
#' 
#' 
#' 
#' 
#' 
#' #################################################################################
#' # Calculate the number of states (i.e., number of difference geographic ranges, 
#' # i.e. number of different combinations of presence/absence in areas) based on
#' # the number of areas
#' #################################################################################
#' numstates_from_numareas(numareas=3, maxareas=3, include_null_range=FALSE)
#' numstates_from_numareas(numareas=3, maxareas=3, include_null_range=TRUE)
#' numstates_from_numareas(numareas=3, maxareas=2, include_null_range=TRUE)
#' numstates_from_numareas(numareas=3, maxareas=1, include_null_range=TRUE)
#' numstates_from_numareas(numareas=7, maxareas=7, include_null_range=TRUE)
#' numstates_from_numareas(numareas=7, maxareas=2, include_null_range=TRUE)
#' numstates_from_numareas(numareas=8, maxareas=8, include_null_range=TRUE)
#' numstates_from_numareas(numareas=8, maxareas=2, include_null_range=TRUE)
#' numstates_from_numareas(numareas=20, maxareas=20, include_null_range=TRUE)
#' numstates_from_numareas(numareas=20, maxareas=2, include_null_range=TRUE)
#' numstates_from_numareas(numareas=20, maxareas=3, include_null_range=TRUE)
#' 
#' 
#' 
#' 
#' #################################################################################
#' # Generate the list of states based on the list of areas
#' # And then generate the continuous-time transition matrix (Q matrix) 
#' # for changes that happen along branches
#' # (the changes that happen at nodes are cladogenesis events)
#' #################################################################################
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
#' # Calculate the list of 0-based indices for each possible 
#' # geographic range, i.e. each combination of areas
#' \dontrun{
#' 
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
#' \dontrun{
#' 
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
#' }
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #################################################################################
#' # Calculate the probability of each (ancestral range) --> 
#' # (Left,Right descendant range pair) directly
#' #################################################################################
#' 
#' #######################################################
#' # Silly example, but which shows the math:
#' #######################################################
#' ca = c(1,2,3,4,5)
#' cb = c(2,2,2,2,2)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=0, y=1, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' # Another silly example, but which shows the normalization effect of specifying 
#' # Rsp_rowsums:
#' ca_1s = c(1,1,1,1,1)
#' cb_1s = c(1,1,1,1,1)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' 
#' # Get the Rsp_rowsums (sum across each row; each row=an ancestral state)
#' Rsp_rowsums = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca_1s, Rcpp_rightprobs=cb_1s, 
#' l=temp_states_indices,  s=0.33, v=0.33, j=0, y=0.33, printmat=TRUE)
#' Rsp_rowsums
#' 
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=0, y=1, 
#' Rsp_rowsums=Rsp_rowsums, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' #######################################################
#' # Silly example, but which shows the math -- redo with same cladogenesis model, 
#' # specified differently
#' #######################################################
#' ca = c(1,2,3,4,5)
#' cb = c(2,2,2,2,2)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=0.33, v=0.33, j=0, y=0.33, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' # Another silly example, but which shows the normalization effect of specifying 
#' # Rsp_rowsums:
#' ca_1s = c(1,1,1,1,1)
#' cb_1s = c(1,1,1,1,1)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' 
#' # Get the Rsp_rowsums (sum across each row; each row=an ancestral state)
#' Rsp_rowsums = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca_1s, Rcpp_rightprobs=cb_1s, 
#' l=temp_states_indices,  s=0.33, v=0.33, j=0, y=0.33, printmat=TRUE)
#' Rsp_rowsums
#' 
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=0.33, v=0.33, j=0, y=0.33, 
#' Rsp_rowsums=Rsp_rowsums, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' 
#' 
#' #######################################################
#' # Silly example, but which shows the math -- redo with different cladogenesis model
#' # (sympatric-copying only, maximum range size of 1)
#' #######################################################
#' ca = c(1,2,3,4,5)
#' cb = c(2,2,2,2,2)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=0, v=0, j=0, y=1, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' # Another silly example, but which shows the normalization effect of specifying 
#' # Rsp_rowsums:
#' ca_1s = c(1,1,1,1,1)
#' cb_1s = c(1,1,1,1,1)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' 
#' # Get the Rsp_rowsums (sum across each row; each row=an ancestral state)
#' Rsp_rowsums = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca_1s, Rcpp_rightprobs=cb_1s, 
#' l=temp_states_indices,  s=0, v=0, j=0, y=1, printmat=TRUE)
#' Rsp_rowsums
#' 
#' # Note that you get NaNs because some of your states (2 areas) are impossible on 
#' # this model
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=0, v=0, j=0, y=1, 
#' Rsp_rowsums=Rsp_rowsums, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' 
#' 
#' #######################################################
#' # Silly example, but which shows the math -- redo with different 
#' # cladogenesis model (BayArea, sympatric-copying only)
#' #######################################################
#' ca = c(1,2,3,4,5)
#' cb = c(2,2,2,2,2)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' numareas = 3
#' 	maxent01y = matrix(0, nrow = numareas, ncol = numareas)
#' 	maxent01y[, 1] = seq(1, numareas)
#' 	maxent01y[2:3, 2] = seq(2, numareas)
#' 	maxent01y[3, 3] = seq(3, numareas)
#' 
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=0, v=0, j=0, y=1, maxent01y=maxent01y, 
#' max_minsize_as_function_of_ancsize=rep(3,numareas), printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' # Another silly example, but which shows the normalization effect of specifying 
#' # Rsp_rowsums:
#' ca_1s = c(1,1,1,1,1)
#' cb_1s = c(1,1,1,1,1)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' 
#' # Get the Rsp_rowsums (sum across each row; each row=an ancestral state)
#' Rsp_rowsums = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca_1s, Rcpp_rightprobs=cb_1s, 
#' l=temp_states_indices, s=0, v=0, j=0, y=1, maxent01y=maxent01y, 
#' max_minsize_as_function_of_ancsize=rep(3,numareas), printmat=TRUE)
#' Rsp_rowsums
#' 
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=0, v=0, j=0, y=1, maxent01y=maxent01y, 
#' max_minsize_as_function_of_ancsize=rep(3,numareas), printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #######################################################
#' # Actual example
#' #######################################################
#' # When ca & cb are 1s, and s, v, j, y are 1s or 0s:
#' # ...this shows how many possible descendant pairs are possible from each 
#' # possible ancestor, under the model
#' # i.e., how many specific cladogenesis scenarios are possible from each 
#' # possible ancestor
#' # This is the LAGRANGE model
#' ca = c(1,1,1,1,1)
#' cb = c(1,1,1,1,1)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' Rsp_rowsums = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, 
#' l=temp_states_indices, s=1, v=1, j=0, y=1, printmat=TRUE)
#' Rsp_rowsums
#' 
#' # Calculate likelihoods of ancestral states if probabilities of each state
#' # at the base of the left and right branches ARE equal
#' # WITHOUT the weights correction
#' ca = c(0.2,0.2,0.2,0.2,0.2)
#' cb = c(0.2,0.2,0.2,0.2,0.2)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=0, y=1, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' # WITH the weights correction
#' ca = c(0.2,0.2,0.2,0.2,0.2)
#' cb = c(0.2,0.2,0.2,0.2,0.2)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=0, y=1, 
#' Rsp_rowsums=Rsp_rowsums, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' 
#' 
#' # Calculate likelihoods of ancestral states if probabilities of each state
#' # at the base of the left and right branches are NOT equal
#' ca = c(0.05,0.1,0.15,0.2,0.5)
#' cb = c(0.05,0.1,0.15,0.2,0.5)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=0, y=1, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' # WITH the weights correction
#' # Calculate likelihoods of ancestral states if probabilities of each state
#' # at the base of the left and right branches ARE equal
#' ca = c(0.05,0.1,0.15,0.2,0.5)
#' cb = c(0.05,0.1,0.15,0.2,0.5)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=0, y=1, 
#' Rsp_rowsums=Rsp_rowsums, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' # Calculate likelihoods of ancestral states if probabilities of each state
#' # at the base of the left and right branches are NOT equal
#' ca = c(0.05,0.1,0.15,0.2,0.5)
#' cb = c(0.05,0.1,0.15,0.2,0.5)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=0, y=1, 
#' Rsp_rowsums=Rsp_rowsums, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #######################################################
#' # Actual example -- for another model (allowing jump dispersal equal probability 
#' # with the rest)
#' #######################################################
#' # When ca & cb are 1s, and s, v, j, y are 1s or 0s:
#' # ...this shows how many possible descendant pairs are possible from each possible 
#' # ancestor, under the model
#' # i.e., how many specific cladogenesis scenarios are possible from each possible 
#' # ancestor
#' # This is the LAGRANGE model
#' ca = c(1,1,1,1,1)
#' cb = c(1,1,1,1,1)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' Rsp_rowsums = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, 
#' l=temp_states_indices, s=1, v=1, j=1, y=1, printmat=TRUE)
#' Rsp_rowsums
#' 
#' # Calculate likelihoods of ancestral states if probabilities of each state
#' # at the base of the left and right branches ARE equal
#' # WITHOUT the weights correction
#' ca = c(0.2,0.2,0.2,0.2,0.2)
#' cb = c(0.2,0.2,0.2,0.2,0.2)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=1, y=1, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' # WITH the weights correction
#' ca = c(0.2,0.2,0.2,0.2,0.2)
#' cb = c(0.2,0.2,0.2,0.2,0.2)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=1, y=1, 
#' Rsp_rowsums=Rsp_rowsums, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' 
#' 
#' # Calculate likelihoods of ancestral states if probabilities of each state
#' # at the base of the left and right branches are NOT equal
#' ca = c(0.05,0.1,0.15,0.2,0.5)
#' cb = c(0.05,0.1,0.15,0.2,0.5)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=1, y=1, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' # WITH the weights correction
#' # Calculate likelihoods of ancestral states if probabilities of each state
#' # at the base of the left and right branches ARE equal
#' ca = c(0.05,0.1,0.15,0.2,0.5)
#' cb = c(0.05,0.1,0.15,0.2,0.5)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=1, y=1, 
#' Rsp_rowsums=Rsp_rowsums, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' # Calculate likelihoods of ancestral states if probabilities of each state
#' # at the base of the left and right branches are NOT equal
#' ca = c(0.05,0.1,0.15,0.2,0.5)
#' cb = c(0.05,0.1,0.15,0.2,0.5)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=1, y=1, Rsp_rowsums=Rsp_rowsums, 
#' printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #######################################################
#' # Calculate the sums of each row (i.e. for each ancestral state) -- 
#' # changes only based on the model
#' #######################################################
#' # Standard LAGRANGE model
#' # Rcpp_leftprobs=ca, Rcpp_rightprobs=cb are irrelevant except for length,
#' # rcpp_calc_anclikes_sp_rowsums() actually treats them as arrays of 1s
#' # if s, v, j, y are 1s or 0s, then Rsp_rowsums = counts of the events
#' ca = c(0.05,0.1,0.15,0.2,0.5)
#' cb = c(0.05,0.1,0.15,0.2,0.5)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' Rsp_rowsums = rcpp_calc_anclikes_sp_rowsums(Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, 
#' l=temp_states_indices, s=1, v=1, j=0, y=1, printmat=TRUE)
#' Rsp_rowsums
#' 
#' # Standard LAGRANGE model, adding jump dispersal
#' ca = c(0.05,0.1,0.15,0.2,0.5)
#' cb = c(0.05,0.1,0.15,0.2,0.5)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' Rsp_rowsums = rcpp_calc_anclikes_sp_rowsums(Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, 
#' l=temp_states_indices, s=1, v=1, j=1, y=1, printmat=TRUE)
#' Rsp_rowsums
#' 
#' # Same models, parameterized differently
#' # Allowing jump dispersal to areas outside of the ancestral range
#' Rsp_rowsums = rcpp_calc_anclikes_sp_rowsums(Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, 
#' l=temp_states_indices, s=0.5, v=0.5, j=0, y=0.5, printmat=TRUE)
#' Rsp_rowsums
#' 
#' # Allowing jump dispersal to areas outside of the ancestral range
#' Rsp_rowsums = rcpp_calc_anclikes_sp_rowsums(Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, 
#' l=temp_states_indices, s=0.5, v=0.5, j=0.5, y=0.5, printmat=TRUE)
#' Rsp_rowsums
#' 
#' 
#' 
#' 
#' 
#' 
#' #######################################################
#' # The relative weights of the different types of cladogenesis events doesn't matter, 
#' # if the correction factor is included
#' #######################################################
#' 
#' # LAGRANGE+founder-event speciation
#' # Calculate likelihoods of ancestral states if probabilities of each state
#' # at the base of the left and right branches are NOT equal
#' # s,v,j,y weights set to 1
#' ca = c(0.05,0.1,0.15,0.2,0.5)
#' cb = c(0.05,0.1,0.15,0.2,0.5)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' uncorrected_condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(
#' Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=1, y=1, 
#' printmat=TRUE)
#' uncorrected_condlike_of_data_for_each_ancestral_state
#' 
#' Rsp_rowsums = rcpp_calc_anclikes_sp_rowsums(Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, 
#' l=temp_states_indices, s=1, v=1, j=1, y=1, printmat=TRUE)
#' Rsp_rowsums
#' 
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=1, y=1, 
#' Rsp_rowsums=Rsp_rowsums, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' 
#' # s,v,j,y weights set to 0.5
#' ca = c(0.05,0.1,0.15,0.2,0.5)
#' cb = c(0.05,0.1,0.15,0.2,0.5)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' uncorrected_condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(
#' Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, l=temp_states_indices, s=0.5, v=0.5, j=0.5, 
#' y=0.5, printmat=TRUE)
#' uncorrected_condlike_of_data_for_each_ancestral_state
#' 
#' Rsp_rowsums = rcpp_calc_anclikes_sp_rowsums(Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, 
#' l=temp_states_indices, s=0.5, v=0.5, j=0.5, y=0.5, printmat=TRUE)
#' Rsp_rowsums
#' 
#' condlike_of_data_for_each_ancestral_state = rcpp_calc_anclikes_sp(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=0.5, v=0.5, j=0.5, y=0.5, 
#' Rsp_rowsums=Rsp_rowsums, printmat=TRUE)
#' condlike_of_data_for_each_ancestral_state
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #######################################################
#' # For large state spaces (many areas, a great many possible geographic ranges i.e. 
#' # combinations of areas),
#' # rcpp_calc_anclikes_sp() gets slow, even with C++ implementation, as it loops 
#' # through every possible combination
#' # of ancestral and descendant states.  rcpp_calc_anclikes_sp_COOprobs() is a partial 
#' # speedup which takes various shortcuts.
#' #
#' # Instead of having the weights/probabilites represented internally, and producing 
#' # the conditional likelihoods as output, 
#' # rcpp_calc_anclikes_sp_COOprobs() produces 3 lists, giving the coordinates of 
#' # nonzero cells in the transition matrix.
#' # 
#' # List #1: 0-based index of states on the Left branch, for each of the ancestral 
#' # states
#' # List #2: 0-based index of states on the Right branch, for each of the ancestral 
#' # states
#' # List #3: Weight of each transition, for each of the ancestral states
#' #######################################################
#' 
#' 
#' # Silly example, but which shows the math:
#' ca = c(1,2,3,4,5)
#' cb = c(2,2,2,2,2)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' list_weights_of_transitions = rcpp_calc_anclikes_sp_COOprobs(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=0, y=1, printmat=TRUE)
#' list_weights_of_transitions
#' 
#' # List #1: 0-based index of states on the Left branch, for each of the ancestral states
#' # List #2: 0-based index of states on the Right branch, for each of the ancestral 
#' # states
#' # List #3: Weight of each transition, for each of the ancestral states
#' 
#' 
#' 
#' # Get the Rsp_rowsums (sums of the rows of the cladogenesis P matrix)
#' # Set the weights to 1
#' ca = c(1,1,1,1,1)
#' cb = c(1,1,1,1,1)
#' COO_probs_list_for_rowsums = rcpp_calc_anclikes_sp_COOprobs(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=0, y=1, printmat=TRUE)
#' COO_probs_list_for_rowsums
#' Rsp_rowsums = sapply(X=COO_probs_list_for_rowsums[[3]], FUN=sum)
#' Rsp_rowsums
#' 
#' 
#' # Calculate likelihoods of ancestral states if probabilities of each state
#' # at the base of the left and right branches ARE equal
#' ca = c(0.2,0.2,0.2,0.2,0.2)
#' cb = c(0.2,0.2,0.2,0.2,0.2)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' COO_weights_list = rcpp_calc_anclikes_sp_COOprobs(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=0, y=1, printmat=TRUE)
#' COO_weights_list
#' 
#' COO_weights_list_rowsums = sapply(X=COO_probs_list_for_rowsums[[3]], FUN=sum)
#' COO_weights_list_rowsums
#' 
#' # To see the transitional probabilities for each ancestral state, under the model:
#' COO_format_transition_probability_matrix = COO_probs_list_for_rowsums
#' 
#' for (i in 1:length(COO_format_transition_probability_matrix[[3]]))
#' 	{
#' 	COO_format_transition_probability_matrix[[3]][[i]] = 
#'  COO_format_transition_probability_matrix[[3]][[i]] / COO_weights_list_rowsums[[i]]
#' 	}
#' COO_format_transition_probability_matrix
#' 
#' # And you can see that the probabilities now sum to 1 for each row
#' sapply(X=COO_format_transition_probability_matrix[[3]], FUN=sum)
#' 
#' 
#' 
#' 
#' # The following is what you would get for the conditional likelihoods of the 
#' # data given each ancestral state, WITHOUT making each row
#' # of the transition matrix sum to 1
#' uncorrected_COO_condlikes_list = rcpp_calc_anclikes_sp_using_COOprobs(
#' Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, RCOO_left_i_list=COO_weights_list[[1]], 
#' RCOO_right_j_list=COO_weights_list[[2]], RCOO_probs_list=COO_weights_list[[3]], 
#' Rsp_rowsums=rep(1,length(ca)), printmat=TRUE)
#' uncorrected_COO_condlikes_list
#' 
#' # This is what you get if you correct, so that each row sums to 1, 
#' # using the sums of the rows to normalize
#' corrected_COO_condlikes_list = rcpp_calc_anclikes_sp_using_COOprobs(
#' Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, RCOO_left_i_list=COO_weights_list[[1]], 
#' RCOO_right_j_list=COO_weights_list[[2]], RCOO_probs_list=COO_weights_list[[3]], 
#' Rsp_rowsums=COO_weights_list_rowsums, printmat=TRUE)
#' corrected_COO_condlikes_list
#' 
#' 
#' 
#' 
#' 
#' #################################################################################
#' # rcpp_calc_anclikes_sp_COOweights_faster():
#' # An even faster method "intelligently" looks for allowed transitions with 
#' # nonzero weights
#' # The output is stored in 4 lists / columns in COO_weights_columnar:
#' # List #1. 0-based index of ancestral states / geographic ranges
#' # List #2. 0-based index of Left descendant states / geographic ranges
#' # List #3. 0-based index of Right descendant states / geographic ranges
#' # List #4. Weight (or probability, if each weight has been divided by the sum 
#' # of the weights for the row) of the
#' #    transition specified by that cell.
#' #################################################################################
#' 
#' 
#' 
#' 
#' # Get the Rsp_rowsums (sums of the rows of the cladogenesis P matrix)
#' # Set the weights to 1
#' ca = c(1,1,1,1,1) # ca and cb don't matter here, since we are just calculating 
#' # the weights
#' cb = c(1,1,1,1,1)
#' temp_states_indices = list(c(0), c(1), c(2), c(0,1), c(1,2))
#' COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(
#' Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, 
#' j=0, y=1, printmat=TRUE)
#' COO_weights_columnar
#' 
#' # List #1. 0-based index of ancestral states / geographic ranges
#' # List #2. 0-based index of Left descendant states / geographic ranges
#' # List #3. 0-based index of Right descendant states / geographic ranges
#' # List #4. Weight (or probability, if each weight has been divided by the 
#' # sum of the weights for the row) of the
#' #    transition specified by that cell.
#' 
#' # Calculate the sums of the weights for each row/ancestral state
#' numstates = 1+max(sapply(X=COO_weights_columnar, FUN=max)[1:3])
#' Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(
#' COO_weights_columnar=COO_weights_columnar, numstates=numstates)
#' Rsp_rowsums
#' 
#' 
#' 
#' # Silly example, but which shows the math:
#' ca = c(1,2,3,4,5)
#' cb = c(2,2,2,2,2)
#' 
#' # WITHOUT using appropriate correction (correction = dividing by the 
#' # sum of the weights for each row)
#' uncorrected_condlikes_list = rcpp_calc_splitlikes_using_COOweights_columnar(
#' Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, COO_weights_columnar=COO_weights_columnar, 
#' Rsp_rowsums=rep(1,numstates), printmat=TRUE)
#' uncorrected_condlikes_list
#' 
#' # WITH using appropriate correction (correction = dividing by the sum
#' # of the weights for each row)
#' corrected_condlikes_list = rcpp_calc_splitlikes_using_COOweights_columnar(
#' Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, COO_weights_columnar=COO_weights_columnar, 
#' Rsp_rowsums, printmat=TRUE)
#' corrected_condlikes_list
#' 
#' 
#' 
#' 
#' # Calculate likelihoods of ancestral states if probabilities of each state
#' # at the base of the left and right branches ARE equal
#' # Get the Rsp_rowsums (sums of the rows of the cladogenesis P matrix)
#' # Set the weights to 1
#' ca = c(0.2,0.2,0.2,0.2,0.2) # ca and cb don't matter here, since we are just 
#' # calculating the weights
#' cb = c(0.2,0.2,0.2,0.2,0.2)
#' COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs=ca, 
#' Rcpp_rightprobs=cb, l=temp_states_indices, s=1, v=1, j=0, y=1, printmat=TRUE)
#' COO_weights_columnar
#' 
#' # List #1. 0-based index of ancestral states / geographic ranges
#' # List #2. 0-based index of Left descendant states / geographic ranges
#' # List #3. 0-based index of Right descendant states / geographic ranges
#' # List #4. Weight (or probability, if each weight has been divided by the 
#' # sum of the weights for the row) of the
#' #    transition specified by that cell.
#' 
#' # Calculate the sums of the weights for each row/ancestral state
#' numstates = 1+max(sapply(X=COO_weights_columnar, FUN=max)[1:3])
#' Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(
#' COO_weights_columnar=COO_weights_columnar, numstates=numstates)
#' Rsp_rowsums
#' 
#' 
#' # WITHOUT using appropriate correction (correction = dividing by 
#' # the sum of the weights for each row)
#' uncorrected_condlikes_list = rcpp_calc_splitlikes_using_COOweights_columnar(
#' Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, COO_weights_columnar=COO_weights_columnar, 
#' Rsp_rowsums=rep(1,numstates), printmat=TRUE)
#' uncorrected_condlikes_list
#' 
#' # WITH using appropriate correction (correction = dividing by the sum of the 
#' # weights for each row)
#' corrected_condlikes_list = rcpp_calc_splitlikes_using_COOweights_columnar(
#' Rcpp_leftprobs=ca, Rcpp_rightprobs=cb, COO_weights_columnar=COO_weights_columnar, 
#' Rsp_rowsums, printmat=TRUE)
#' corrected_condlikes_list
#' 
 
