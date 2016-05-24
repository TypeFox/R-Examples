#' @title Probability of split of the tree prior.
#' @description Probability of split of the tree prior.
#' @param alpha base parameter of the tree prior, \eqn{\alpha \in [0,1)}.
#' @param beta power parameter of the tree prior, \eqn{\beta \geq 0}.
#' @param depth depth of the current node, \eqn{depth \geq 0}.
#' @return Returns the probability of split.
#' @examples
#' p_split(.95,.5)
#' p_split(.95,.5,1)
#' p_split(.95,.5,2)
#' @references Chipman, H. A., George, E. I. et McCulloch, R. E. (1998). Bayesian CART model search. Journal of the American Statistical Association, 93(443), 935-948.
#' @export
"p_split"

#' @title Expected value of the number of bottom nodes in the unrealistic case where we assume that the number of variables and possible splits are infinite (therefore P(T) is not dependent on the design matrix X) and \eqn{\beta=0} (Case #1).
#' @description Expected value of the number of bottom nodes in the unrealistic case where we assume that the number of variables and possible splits are infinite (therefore P(T) is not dependent on the design matrix X) and \eqn{\beta=0} (Case #1).
#' @param alpha base parameter of the tree prior, \eqn{alpha \in [0,1)}.
#' @return Returns the expected value of the number of bottom nodes.
#' @examples
#' E_alpha(.30)
#' E_alpha(.45)
#' E_alpha(.499)
#' E_alpha(.75)
#' @seealso \code{\link{Var_alpha}}
#' @references Jolicoeur-Martineau, A. (Currently in revision, expected 2016) \emph{Etude d'une loi a priori pour les arbres binaires de regression} (\emph{Study on the prior distribution of binary regression trees}) (Master thesis). UQAM university.
#' @export
"E_alpha"

#' @title Variance of the number of bottom nodes in the unrealistic case where we assume that the number of variables and possible splits are infinite (therefore P(T) is not dependent on the design matrix X) and \eqn{\beta=0} (Case #1).
#' @description Variance of the number of bottom nodes in the unrealistic case where we assume that the number of variables and possible splits are infinite (therefore P(T) is not dependent on the design matrix X) and \eqn{\beta=0} (Case #1).
#' @param alpha base parameter of the tree prior, \eqn{\alpha \in [0,1)}.
#' @return Returns the variance of the number of bottom nodes.
#' @examples
#' Var_alpha(.30)
#' Var_alpha(.45)
#' Var_alpha(.499)
#' Var_alpha(.75)
#' @seealso \code{\link{E_alpha}}
#' @references Jolicoeur-Martineau, A. (Currently in revision, expected 2016) \emph{Etude d'une loi a priori pour les arbres binaires de regression} (\emph{Study on the prior distribution of binary regression trees}) (Master thesis). UQAM university.
#' @export
"Var_alpha"

#' @title Number of bottom nodes and depth in the unrealistic case where we assume that the number of variables and possible splits are infinite (therefore P(T) is not dependent on the design matrix X) (Case #2).
#' @description Generate a tree and returns the number of bottom nodes and depth in the unrealistic case where we assume that the number of variables and possible splits are infinite (therefore P(T) is not dependent on the design matrix X) (Case #2).
#' @param alpha base parameter of the tree prior, \eqn{\alpha \in [0,1)}.
#' @param beta power parameter of the tree prior, \eqn{beta \geq 0}.
#' @param depth depth of the current node, \eqn{depth \geq 0}.
#' @return Returns a vector containing the number of bottom nodes and depth.
#' @examples
#' NumBotMaxDepth_inf(.95,.5)
#' @seealso \code{\link{NumBotMaxDepth}}, \code{\link{NumBotMaxDepthX}}
#' @export
"NumBotMaxDepth_inf"

#' @title Simulation of the tree prior in the unrealistic case where we assume that the number of variables and possible splits are infinite (therefore P(T) is not dependent on the design matrix X) (Case #2).
#' @description Generate \eqn{n_{iter}} trees from the prior distribution in the unrealistic case where we assume that the number of variables and possible splits are infinite (therefore P(T) is not dependent on the design matrix X) (Case #2).
#' @param alpha base parameter of the tree prior, \eqn{\alpha \in [0,1)}.
#' @param beta power parameter of the tree prior, \eqn{beta \geq 0}.
#' @param n_iter number of trees to generate, \eqn{n_{iter}>0}.
#' @return Returns a list containing, in the following order: the mean number of bottom nodes, the standard deviation of the number of bottom nodes, the mean of the depth, the standard deviation of the depth and a data.frame of vectors \eqn{(b_i,d_i)}, where \eqn{b_i} is the number of bottom nodes and \eqn{d_i} is the depth of the \eqn{i}th generated tree (\eqn{i=1, \ldots ,n_{iter}}).
#' @examples
#' results = BayesTreePriorOrthogonalInf(.95,.5)
#' @seealso \code{\link{BayesTreePriorOrthogonal}}, \code{\link{BayesTreePriorNotOrthogonal}}
#' @export
"BayesTreePriorOrthogonalInf"

#' @title Number of bottom nodes and depth in the case where the design matrix X is orthogonal (Case #3).
#' @description Generate a tree and returns the number of bottom nodes and depth in the case where the design matrix X is orthogonal (Case #3).
#' @param alpha base parameter of the tree prior, \eqn{\alpha \in [0,1)}.
#' @param beta power parameter of the tree prior, \eqn{beta \geq 0}.
#' @param x_size number of possible splits, \eqn{x_{size}>0}.
#' @param depth depth of the current node, \eqn{depth \geq 0}.
#' @return Returns a vector containing the number of bottom nodes and depth
#' @examples
#' NumBotMaxDepth(.95,.5,500)
#' @seealso \code{\link{NumBotMaxDepth_inf}}, \code{\link{NumBotMaxDepthX}}
#' @export
"NumBotMaxDepth"

#' @title Simulation of the tree prior in the case where the design matrix X is orthogonal (Case #3).
#' @description Generate \eqn{n_{iter}} trees from the prior distribution in the case where the design matrix X is orthogonal (Case #3).
#' @param alpha base parameter of the tree prior, \eqn{\alpha \in [0,1)}.
#' @param beta power parameter of the tree prior, \eqn{beta \geq 0}.
#' @param n_obs number of unique observations, \eqn{n_{obs}>1}.
#' @param n_iter number of trees to generate, \eqn{n_{iter}>0}.
#' @return Returns a list containing, in the following order: the mean number of bottom nodes, the standard deviation of the number of bottom nodes, the mean of the depth, the standard deviation of the depth and a data.frame of vectors \eqn{(b_i,d_i)}, where \eqn{b_i} is the number of bottom nodes and \eqn{d_i} is the depth of the \eqn{i}th generated tree (\eqn{i=1, \ldots ,n_{iter}}).
#' @examples
#' results1 = BayesTreePriorOrthogonal(.95,.5, 100)
#' results2 = BayesTreePriorOrthogonal(.95,.5, 250)
#' @seealso \code{\link{BayesTreePriorOrthogonalInf}}, \code{\link{BayesTreePriorNotOrthogonal}}
#' @export
"BayesTreePriorOrthogonal"

#' @title Unique splits that leads to children with more than \eqn{minpart} nodes.
#' @description Unique splits that leads to children with more than \eqn{minpart} nodes.
#' @param x vector containing the observations of a variable.
#' @param minpart minimum number of observations in the children nodes.
#' @param MIA set to TRUE if you want Missing Incorporated in Attributes (MIA) imputation to be used.
#' @return If \eqn{MIA} is TRUE and \eqn{minpart>1}, the possible splits could be different depending on whether we transfer the NAs to the left child or the right child; if this is the case then the function returns a list \eqn{(v1,v2)}, where \eqn{v1} is the vector containing the unique splits that leads to \eqn{minpart} nodes when transferring the NAs to the left child and \eqn{v2} is the vector containing the unique splits that leads to children with more than \eqn{minpart} nodes when transferring the NAs to the left child. Otherwise, it returns the vector containing the unique splits that leads to children with more than \eqn{minpart} nodes.
#' @examples
#' GetListUniqueSplits(c(1,4,7,3,0,2,2,3,4,7,7,7),minpart=1)
#' GetListUniqueSplits(c(1,4,7,3,0,2,2,3,4,7,7,7),minpart=3)
#' GetListUniqueSplits(c(1,4,7,3,0,2,2,3,4,7,7,7,NA,NA,NA),minpart=1, MIA=TRUE)
#' GetListUniqueSplits(c(1,4,7,3,0,2,2,3,4,7,7,7,NA,NA,NA),minpart=3, MIA=TRUE)
#' @export
"GetListUniqueSplits"

#' @title Number of bottom nodes and depth in the general case (Case #4).
#' @description Generate a tree and returns the number of bottom nodes and depth in the general case (Case #4).
#' @param alpha base parameter of the tree prior, \eqn{\alpha \in [0,1)}.
#' @param beta power parameter of the tree prior, \eqn{beta \geq 0}.
#' @param X data.frame of the design matrix.
#' @param depth depth of the current node, \eqn{depth \geq 0}.
#' @param minpart the minimum number of bottom nodes required in one of the child to be able to split, \eqn{minpart>0}.
#' @param pvars vector of probabilities for the choices of variables to split (Will automatically be normalized so that the sum equal to 1). It must be twice as large as the number of variables when \eqn{missingdummy} is TRUE.
#' @param MIA set to TRUE if you want Missing Incorporated in Attributes (MIA) imputation to be used.
#' @param missingdummy set to TRUE if you have dummy coded the NAs.
#' @return Returns a vector containing the number of bottom nodes and depth
#' @examples
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'     x1 = MASS::mcycle$times
#'     x1[sample(1:length(x1), 20)] <- NA
#'     x2= MASS::mcycle$accel
#'     x2[sample(1:length(x2), 20)] <- NA
#'     X = cbind(x1, x2)
#'     results1 = NumBotMaxDepthX(.95,.5, data.frame(X), minpart=5)
#'     X_dummies = is.na(X) + 0
#'     results2 = NumBotMaxDepthX(.95,.5, data.frame(cbind(X,X_dummies)), minpart=5, MIA=TRUE, 
#'     missingdummy=TRUE)
#' }
#' @seealso \code{\link{NumBotMaxDepth_inf}}, \code{\link{NumBotMaxDepth}}
#' @references Twala, B. E. T. H., Jones, M. C., & Hand, D. J. (2008). \emph{Good methods for coping with missing data in decision trees.} Pattern Recognition Letters, 29(7), 950-956.
#' @export
"NumBotMaxDepthX"

#' @title Simulation of the tree prior in the general case (Case #4).
#' @description Generate \eqn{n_{iter}} trees from the prior distribution in the general case (Case #4).
#' @param alpha base parameter of the tree prior, \eqn{\alpha \in [0,1)}.
#' @param beta power parameter of the tree prior, \eqn{\beta \geq 0}.
#' @param X data.frame of the design matrix.
#' @param n_iter number of trees to generate, \eqn{n_{iter}>0}.
#' @param minpart the minimum number of bottom nodes required in one of the child to be able to split, \eqn{minpart>0}.
#' @param pvars vector of probabilities for the choices of variables to split (Will automatically be normalized so that the sum equal to 1). It must be twice as large as the number of variables when \eqn{missingdummy} is TRUE.
#' @param MIA set to TRUE if you want Missing Incorporated in Attributes (MIA) imputation to be used.
#' @param missingdummy set to TRUE if you have dummy coded the NAs.
#' @return Returns a list containing, in the following order: the mean number of bottom nodes, the standard deviation of the number of bottom nodes, the mean of the depth, the standard deviation of the depth and a data.frame of vectors \eqn{(b_i,d_i)}, where \eqn{b_i} is the number of bottom nodes and \eqn{d_i} is the depth of the \eqn{i}th generated tree (\eqn{i=1, \ldots ,n_{iter}}).
#' @examples
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'     x1 = MASS::mcycle$times
#'     x1[sample(1:length(x1), 20)] <- NA
#'     x2= MASS::mcycle$accel
#'     x2[sample(1:length(x2), 20)] <- NA
#'     X = cbind(x1, x2)
#'     results1 = BayesTreePriorNotOrthogonal(.95,.5, data.frame(X), minpart=5)
#'     X_dummies = is.na(X) + 0
#'     results2 = BayesTreePriorNotOrthogonal(.95,.5, data.frame(cbind(X,X_dummies)), minpart=5, 
#'     MIA=TRUE, missingdummy=TRUE)
#' }
#' @seealso \code{\link{BayesTreePriorOrthogonalInf}}, \code{\link{BayesTreePriorOrthogonal}}
#' @export
"BayesTreePriorNotOrthogonal"

#' @title Simulation of the tree prior.
#' @description This is the main function to use for simulating from the prior. There are 4 cases : 
#' \itemize{
#'  \item{Case #1: }{Unrealistic case where we assume that the number of variables and possible splits are infinite (therefore \eqn{P(T)} is not dependent on the design matrix X) and \eqn{\beta=0}}
#'  \item{Case #2: }{Unrealistic case where we assume that the number of variables and possible splits are infinite (therefore \eqn{P(T)} is not dependent on the design matrix X)}
#'  \item{Case #3: }{Design matrix X is orthogonal}
#'  \item{Case #4: }{General case}
#' }
#' Case #1 will be used if no design matrix X or number of observations is given and \eqn{\beta = 0}. Case #2 will be used if no design matrix X or number of observations is given and \eqn{\beta \neq 0}. Case #3 will be used if no design matrix X is given but the number of observations is given. Case #4 will be used if the design matrix X is given. Note that case #4 is always slower, so if your design matrix is orthogonal, it would be advisable to enter the number of uniques observations rather than the design matrix X, to be able to use case #3.
#' 
#' @param alpha base parameter of the tree prior, \eqn{\alpha \in [0,1)}.
#' @param beta power parameter of the tree prior, \eqn{\beta \geq 0}.
#' @param X data.frame of the design matrix (Required for case #4).
#' @param n_iter number of trees to generate, \eqn{n_{iter}>0} (Used for case #2, #3 or #4).
#' @param n_obs number of unique observations, \eqn{n_{obs}>1} (Required for case #3).
#' @param minpart the minimum number of bottom nodes required in one of the child to be able to split, \eqn{minpart>0}.
#' @param package a optional string that can take the following values : "BayesTree", "tgp" or "bartMachine". It forces the algorithm to use the default paramameters used by the package specified (\eqn{minpart=5} for BayesTree, \eqn{minpart = max(c(10,dim(X)[2]+1))} for tgp and \eqn{minpart=1} for bartMachine). 
#' @param pvars vector of probabilities for the choices of variables to split (Will automatically be normalized so that the sum equal to 1). It must be twice as large as the number of variables when \eqn{missingdummy} is TRUE.
#' @param MIA set to TRUE if you want Missing Incorporated in Attributes (MIA) imputation to be used.
#' @param missingdummy set to TRUE if you want the NAs to be dummy coded.
#' @return In case #1, it returns a list containing, in the following order: the expectation and the variance of the number of bottom nodes. In cases #2, #3 or #4, it returns a list containing, in the following order: the mean number of bottom nodes, the standard deviation of the number of bottom nodes, the mean of the depth, the standard deviation of the depth and a data.frame of vectors \eqn{(b_i,d_i)}, where \eqn{b_i} is the number of bottom nodes and \eqn{d_i} is the depth of the \eqn{i}th generated tree (\eqn{i=1, \ldots ,n_{iter}}).
#' @examples
#' #Case 1 : Unrealistic case where we assume that the number of var/obs is infinite and beta=0
#' results1 = BayesTreePrior(.45,0)
#' #Case 2 : Unrealistic case where we assume that the number of var/obs is infinite
#' results2 = BayesTreePrior(.95,.5)
#' #Case 3 : Design matrix X is orthogonal
#' results3 = BayesTreePrior(.95,.5,n_obs=150)
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'     #Case 4 : General case, without missing data
#'     x1 = MASS::mcycle$times
#'     x2= MASS::mcycle$accel
#'     X = cbind(x1, x2)
#'     results4_nomiss = BayesTreePrior(.95,.5, data.frame(X), minpart=5, package="tgp")
#'     #Case 4 : General case, with missing data
#'     x1[sample(1:length(x1), 20)] <- NA
#'     x2[sample(1:length(x2), 20)] <- NA
#'     X = cbind(x1, x2)
#'     results4_miss = BayesTreePrior(.95,.5, data.frame(X), minpart=5, package="tgp", 
#'     MIA=TRUE, missingdummy=TRUE)
#' }
#' @references
#' Chipman, H. A., George, E. I., & McCulloch, R. E. (1998). \emph{Bayesian CART model search.} Journal of the American Statistical Association, 93(443), 935-948.
#'
#' Gramacy, R. B. (2007). tgp: \emph{an \R package for Bayesian nonstationary, semiparametric nonlinear regression and design by treed Gaussian process models.} Journal of Statistical Software, 19(9), 6.
#'
#' Chipman, H. A., George, E. I., & McCulloch, R. E. (2010). \emph{BART: Bayesian additive regression trees.} The Annals of Applied Statistics, 266-298.
#'
#' Kapelner, A., & Bleich, J. (2013). \emph{bartMachine: A powerful tool for machine learning.} stat, 1050, 8.
#'
#' Twala, B. E. T. H., Jones, M. C., & Hand, D. J. (2008). \emph{Good methods for coping with missing data in decision trees.} Pattern Recognition Letters, 29(7), 950-956.
#'
#' Jolicoeur-Martineau, A. (Currently in revision, expected 2016) \emph{Etude d'une loi a priori pour les arbres binaires de regression} (\emph{Study on the prior distribution of binary regression trees}) (Master thesis). UQAM university.
#' @import stats
#' @export
"BayesTreePrior"

p_split = function(alpha,beta,depth=0){
	return(alpha/((1+depth)^beta))
}

E_alpha = function(alpha){
	if (alpha <.5) return ((1-alpha)/(1-2*alpha))
	else return (Inf)
}

Var_alpha = function(alpha){
	if (alpha <.5) return (((1-alpha)*alpha)/((1-2*alpha)^3))
	else return (Inf)
}

NumBotMaxDepth_inf = function(alpha, beta, depth = 0){
	if (runif(1) <= p_split(alpha, beta, depth)){
		results_left = NumBotMaxDepth_inf(alpha, beta, depth+1)
		results_right = NumBotMaxDepth_inf(alpha, beta, depth+1)
		return(c(results_left[1] + results_right[1], max(results_left[2], results_right[2])))
	} 
	else return(c(1, depth))
}

BayesTreePriorOrthogonalInf = function(alpha, beta, n_iter=500){
	nodes = rep(0,n_iter)
	depth = rep(0,n_iter)
	for (i in 1:n_iter){
		results = NumBotMaxDepth_inf(alpha, beta)
		nodes[i] = results[1]
		depth[i] = results[2]
	}
	return(list(mean_nodes = mean(nodes), sd_nodes = sd(nodes), mean_depth = mean(depth), sd_depth = sd(depth), draws = cbind(nodes,depth)))
}

NumBotMaxDepth = function(alpha, beta, x_size, depth = 0){
	index = sample(1:x_size, 1)
	if (runif(1) <= p_split(alpha, beta, depth)){
		if (x_size == 1) return(c(2, depth+1))	
		else if (index == 1 || index == x_size){
			results = NumBotMaxDepth(alpha, beta, (x_size - 1), depth+1)
			return(c(1 + results[1], results[2]))
		} 
		else{
			results_left = NumBotMaxDepth(alpha, beta, (index - 1), depth+1)
			results_right = NumBotMaxDepth(alpha, beta, (x_size - index), depth+1)
			return(c(results_left[1] + results_right[1], max(results_left[2], results_right[2])))
		} 
	}
	else return(c(1, depth))
}

BayesTreePriorOrthogonal = function(alpha, beta, n_obs, n_iter=500)
{
	nodes = rep(0,n_iter)
	depth = rep(0,n_iter)
	for (i in 1:n_iter){
		results = NumBotMaxDepth(alpha, beta, n_obs-1)
		nodes[i] = results[1]
		depth[i] = results[2]
	}

	return(list(mean_nodes = mean(nodes), sd_nodes = sd(nodes), mean_depth = mean(depth), sd_depth = sd(depth), draws = cbind(nodes,depth)))
}

GetListUniqueSplits = function(x, minpart=1, MIA=FALSE)
{
	if (!MIA){
		if (minpart == 1){
			# Note that sort() automatically remove NAs
			x = unique(sort(x))
			x = x[-length(x)]
			return(x)
		}
		else{
			x_sorted = sort(x)
			# NAs only, return empty vector
			if (length(x_sorted) == 0) return(numeric(0))
			rle = rle(x_sorted)
			uniques = unlist(rle[2])
			if (length(uniques) <= 1) return(numeric(0))
			n_left = cumsum(unlist(rle[1]))
			n_right = length(x) - n_left
			reject = rep(FALSE, length(n_left))
			for (i in 1:length(n_left)){
				if (min(n_left[i],n_right[i]) < minpart) reject[i]=TRUE
				else reject[i]=FALSE
			}
			return(uniques[!reject])
		}
	}
	else{
	# If we are using MIA, we must return the two types of splits [(X<=x or X missing) v.s. (X>x)] and [(X<=x) v.s. (X>x or X missing)] as a list
		if (minpart == 1){
			# Note that sort() automatically remove NAs
			x = unique(sort(x))
			x = x[-length(x)]
			return(x)
		}
		else{
			NA_counts = sum(is.na(x))
			x_sorted = sort(x)
			# NAs only, return empty vector
			if (length(x_sorted) == 0) return(numeric(0))
			rle = rle(x_sorted)
			uniques = unlist(rle[2])
			if (length(uniques) <= 1) return(numeric(0))
			n_left = cumsum(unlist(rle[1]))
			if (NA_counts == 0){
				n_right = length(x) - n_left
				reject = rep(FALSE, length(n_left))
				for (i in 1:length(n_left)){
					if (min(n_left[i],n_right[i]) < minpart) reject[i]=TRUE
					else reject[i]=FALSE
				}
				uniques = uniques[!reject]
				return(uniques)
			}
			else{
				# We need to remove the last unique observation. 
				# Normally this is automatically dealt with when we search for reject=TRUE but here if we split on the last unique observation
				# it could lead to only NAs in one of the right children but the size would be bigger than minpart and we would have reject=FALSE.
				# We already have included dummy variables of NAs for that purposes in the variable so we can't have them as splits in here too.
				n_left = n_left[-length(n_left)]
				uniques = uniques[-length(uniques)]
				n_left_1 = n_left + NA_counts
				n_right_1 = length(x) - n_left_1
				n_left_2 = n_left
				n_right_2 = length(x) - n_left_2

				reject_1 = rep(FALSE, length(n_left_1))
				reject_2 = rep(FALSE, length(n_left_2))
				for (i in 1:length(n_left)){
					if (min(n_left_1[i],n_right_1[i]) < minpart) reject_1[i]=TRUE
					else reject_1[i]=FALSE
					if (min(n_left_2[i],n_right_2[i]) < minpart) reject_2[i]=TRUE
					else reject_2[i]=FALSE
				}
				uniques1 = uniques[!reject_1]
				uniques2 = uniques[!reject_2]
				if (length(uniques1) == 0) return(uniques2)
				else if (length(uniques2) == 0) return(uniques1)
				else return(list(v1 = uniques[!reject_1], v2 = uniques[!reject_2]))
			}
		}
	}
}

NumBotMaxDepthX = function(alpha, beta, X, depth = 0, minpart=1, pvars=NULL, MIA=FALSE, missingdummy=FALSE){
	# If we splits with size < minpart*2, this will necessarely leads to children with less than minpart nodes
	if (dim(X)[1] < minpart*2) return(c(1, depth))
	else if (runif(1) <= p_split(alpha, beta, depth)){
		# List of lists (v1,v2) is (MIA && minpart >1), otherwise list of vectors
		X_list = lapply(as.list(X),GetListUniqueSplits, minpart=minpart, MIA=MIA)
		# Remove variables with no available splits but remember their indexes
		empty_splits = sapply(X_list, function(x) length(x) == 0)
		X_list = X_list[!empty_splits]
		nvars = length(X_list)
		# Can't split. This can happen if you have some observations with the exact same values (collinear matrix) or if no variables leads to children with n > minpart.
		if (nvars == 0) return(c(1, depth))
		index_var = sample(1:nvars, 1, prob=pvars[!empty_splits])
		true_index_var = match(c(index_var),cumsum(!empty_splits))
		if (MIA){
			if (runif(1) <= 0.5){
				NA_to_left = TRUE
				if (class(X_list[[index_var]]) == "list") X_list[[index_var]] = X_list[[index_var]][[1]]
			} 
			else{
				NA_to_left = FALSE
				if (class(X_list[[index_var]]) == "list") X_list[[index_var]] = X_list[[index_var]][[2]]
			}
		}
		nobs = length(X_list[[index_var]])
		index_obs = sample(1:nobs, 1)
		# Assign true or false depending on if the variable is <= the index of observation
		if (MIA){
			# If the variable we split on is NOT a missing dummy variable
			if (!missingdummy || true_index_var <= dim(X)[2]/2){
				# bartMachine package source code shows that they implement MIA so that there is 1/2 probability of NAs going left and 1/2 probability of going right
				# Therefore, we will implement use a coin flip here and we will return the chosen list of splits
				# Case 1 : (X<=x or X missing) v.s. (X>x)
				if (NA_to_left) X_new_index <- ifelse(X[, true_index_var] <= X_list[[index_var]][index_obs] | is.na(X[, true_index_var]), TRUE, FALSE)
				# Case 2 : (X<=x) v.s. (X>x or X missing)
				else X_new_index <- ifelse(X[, true_index_var] <= X_list[[index_var]][index_obs] & !is.na(X[, true_index_var]), TRUE, FALSE)
			}
			#Case 3 : X missing vs X not missing
			else X_new_index <- ifelse(X[, true_index_var] <= X_list[[index_var]][index_obs], TRUE, FALSE)
		}
		else X_new_index <- ifelse(X[, true_index_var] <= X_list[[index_var]][index_obs], TRUE, FALSE)
		X_left = X[X_new_index,]
		X_right = X[!X_new_index,]
		dim_left = dim(X_left)[1]
		dim_right = dim(X_right)[1]
		if (dim_left == 1 && dim_right == 1) c(2, depth+1)
		else if (dim_right == 1){
			results = NumBotMaxDepthX(alpha, beta, X_left, depth+1, minpart=minpart, pvars=pvars, MIA=MIA, missingdummy=missingdummy)
			return(c(1 + results[1], results[2]))
		} 
		else if(dim_left == 1){
			results = NumBotMaxDepthX(alpha, beta, X_right, depth+1, minpart=minpart, pvars=pvars, MIA=MIA, missingdummy=missingdummy)
			return(c(1 + results[1], results[2]))
		} 
		else{
			results_left = NumBotMaxDepthX(alpha, beta, X_left, depth+1, minpart=minpart, pvars=pvars, MIA=MIA, missingdummy=missingdummy)
			results_right = NumBotMaxDepthX(alpha, beta, X_right, depth+1, minpart=minpart, pvars=pvars, MIA=MIA, missingdummy=missingdummy)
			return(c(results_left[1] + results_right[1], max(results_left[2], results_right[2])))
		} 
	}
	else return(c(1, depth))
}

BayesTreePriorNotOrthogonal = function(alpha, beta, X, n_iter=500, minpart=1, pvars=NULL, MIA=FALSE, missingdummy=FALSE)
{
	nodes = rep(0,n_iter)
	depth = rep(0,n_iter)
	for (i in 1:n_iter){
		results = NumBotMaxDepthX(alpha, beta, X, minpart=minpart, pvars=pvars, MIA=MIA, missingdummy=missingdummy)
		nodes[i] = results[1]
		depth[i] = results[2]
	}
	return(list(mean_nodes = mean(nodes), sd_nodes = sd(nodes), mean_depth = mean(depth), sd_depth = sd(depth), draws = cbind(nodes,depth)))
}


BayesTreePrior = function(alpha, beta, X=NULL, n_obs=NULL, n_iter=500, minpart=1, package=NULL, pvars=NULL, MIA=FALSE, missingdummy=FALSE)
{
	if (!is.null(package)){
		if (package=="BayesTree"){
			message("Using default parameters of BayesTree package : minpart = 5")
			minpart = 5
		}
		else if (package=="tgp"){
			message("Using default parameters of tgp package : minpart = max(c(10,dim(X)[2]+1))")
			if (!is.null(X)) minpart = max(c(10,dim(X)[2]+1))
		}
		else if (package=="bartMachine"){
			message("Using default parameters of bartMachine package : minpart = 1")
			minpart = 1
		}
		else if (package!="") warning("Argument package is not recognized, ignoring. Set package=NULL to remove warning.")
	} 
	if (mode(alpha) !="NULL" && (mode(alpha) != "numeric" || alpha < 0)) stop("alpha must be a number bigger or equal to 0")
	if (alpha > 1) warning("alpha should probably be <= 1, alpha > 1 will force the use of the alternative definition of the probability of split : min[1, alpha/((1+depth)^beta)] (Jolicoeur-Martineau, A. (Master thesis in revision, expected 2016))")
	if (mode(beta) !="NULL" && (mode(beta) != "numeric" || beta < 0)) stop("beta must be a number bigger than 0")
	if (mode(X) !="NULL" && class(X) != "data.frame" && class(X) != "matrix" && !is.vector(X)) stop("X must be a data.frame")
	if (mode(n_obs) !="NULL" && (mode(n_obs) != "numeric" || n_obs%%1!=0 || n_obs < 2)) stop("n_obs must be a integer bigger than 1")
	if (mode(n_iter) != "numeric" || n_iter%%1!=0 || n_iter < 1) stop("n_iter must be a integer bigger than 0")
	if (mode(minpart) != "numeric"|| minpart%%1!=0 || minpart < 1) stop("minpart must be a integer bigger than 0")
	if (mode(pvars) !="NULL" && !is.vector(pvars)) stop("pvars must be a vector")
	if(typeof(MIA) != "logical") stop("MIA must be a logical variable")
	if(typeof(missingdummy) != "logical") stop("missingdummy must be a logical variable")

	if (is.null(X)){
		if (is.null(n_obs)){
			if (beta == 0){
				message("Case 1 -> Unrealistic case where we assume that the number of variables and possible splits are infinite (therefore P(T) is not dependent on the design matrix X) and beta=0")
				message("Ignored arguments : n_iter, min_part, pvars, MIA, missingdummy")
				return(list(mean_nodes = E_alpha(alpha), sd_nodes = sqrt(Var_alpha(alpha))))
			}
			else{
				message("Case 2 -> Unrealistic case where we assume that the number of variables and possible splits are infinite (therefore P(T) is not dependent on the design matrix X)")
				message("Ignored arguments : min_part, pvars, MIA, missingdummy")
				return(BayesTreePriorOrthogonalInf(alpha, beta, n_iter))
			}
		}
		else{
			# Add check for orthogonality on X for this case;
			message("Case 3 -> The design matrix is orthogonal")
			message("Ignored arguments : min_part, pvars, MIA, missingdummy")
			return(BayesTreePriorOrthogonal(alpha, beta, n_obs, n_iter))
		}
	}
	else{
   		if(is.vector(X) && minpart==1 && is.null(pvars) && MIA==FALSE){
			message("Case 3 -> The design matrix is orthogonal")
			message("Ignored arguments : min_part, pvars, MIA, missingdummy")
   			return(BayesTreePriorOrthogonal(alpha, beta, n_obs=length(GetListUniqueSplits(X))))
   		}
   		if(dim(X)[2] == 1 && minpart==1 && is.null(pvars) && MIA==FALSE){
			message("Case 3 -> The design matrix is orthogonal")
			message("Ignored arguments : min_part, pvars, MIA, missingdummy")
   			return(BayesTreePriorOrthogonal(alpha, beta, n_obs=length(GetListUniqueSplits(as.vector(as.matrix(X))))))
   		}
		message("Case 4 -> General case")
		X = data.frame(X)
		if (min(sapply(X, is.numeric)) == 0) stop("X contains non-numeric variables. Please dummy code the categorical variables and make sure that every variable is numeric.")
		if (missingdummy == 1 && MIA==1){
				# We could remove dummies that are all zeros, but this won't make much of a difference computational-wise.
				X_dummies = is.na(X) + 0
				X = cbind(X,X_dummies)		
		}
		else if(missingdummy == 0 && MIA==1){
				# Do nothing
		}
		else{
			if (missingdummy == 1 && MIA==0) message("missingdummy=1 is ignored, please set MIA to 1.")
			X = X[complete.cases(X),]
		} 
		return(BayesTreePriorNotOrthogonal(alpha, beta, X, n_iter, minpart, pvars, MIA, missingdummy))
	}
}
