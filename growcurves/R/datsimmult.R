#' Repeated measures for two groups of subjects with two multiple membership (MM) terms
#' 
#' A simulation dataset containing repeated subject measures for 2 treatment arms,
#' (control = 0, treatment = 1), constructed from a four- MM term multiple membership model.
#' MM terms 1 and 3 employ an 'mmi' prior over the set of treatment groups, which are then linked
#' back to subjects with a weight matrix.  MM terms 2 and 4 employ an 'mmcar' prior 
#' over the number of treatment sessions attended, with correlations
#' between adjacent number of sessions equal to 0.25. Subject effects were randomly drawn from 10 clusters
#' with assignment weights/probabilities drawn from a Dirichlet distribution.   Cluster location
#' values were generated from a Gaussian base distribution.
#' 
#' \itemize{
#'   \item y. response.  there are \code{N = 897} total measures for \code{P = 299} subjects
#'   \item subject. subject identifier \code{(1,2,...,299}
#'   \item trt. treatment group identifier of length \code{N} (e.g. \code{(0,0,0,...,1,1,1,...)} , either \code{{0,1}} for control and treatment. 
#'   \item time. times in months for each repeated subject measure of length \code{N}.  There are 3 distinct time points.  e.g. \code{(0,3,6,0,3,6,0,0,3,,,,)}
#'   \item n.random. number of random effects per subject.  Set = 3.
#'   \item n.fix_degree. order of fixed effects.  Set = 2, for quadratic, meaning 3 effects (intercept, slope, quadratic) each, for treatment and control groups. 
#'   \item X The resultant \code{N x n.fixed} fixed effects design matrix.  \code{X = c(1, time, time^2, trt_1,trt_1*time, trt_1*time^2, age, income)}.
#'   \item coefs. true coefficient values for the time-based quadratic fixed effects generated from the \code{trt}, and \code{time} inputs.
#'   \item Z The resultant \code{N x n.random} by-subject random effects design matrix.  \code{Z = c(1, time, time^2}.
#'   \item Z.b The \code{N x 1} matrix  [Z_{1} * b_{1},...,Z_{n}*b_{n}]. where \code{n} equals the number of (unique) subjects.
#'   \item subj.aff. A list object with each term a vector that indexes \code{n[g]} subjects linked to each of \code{g = 1, ..., (G = 4)} CBT therapy groups.  Each group is specialized
#'			to its own multiple membership (MM) term for employment of engine function \code{dpgrowmmult}.  The number of treatment arm subjects for each group are 
#'			\code{(n[1] = 17, n[2] = 18, n[3] = 19, n[4] = 78)} for a total of \code{Paff = 132} subjects assigned to the treatment arm.
#'   \item W.subj.aff.  A list object containing \code{G, n[g] x S[g]} multiple membership weight matrices that together map the \code{P_aff = 132} affected subjects (in each element of the 	
#'			\code{subj.aff} list) to their particular sessions attended within their assigned group.  There are a total of 245 sessions allocated to the \code{G} groups 
#'			as \code{S[1] = 36, S[2] = 40, S[3] = 40, S[4] = 129}.
#'   \item group.  A list object of length equal to the number of MM terms under the 'mmcar' prior - In this case, terms 2 and 4.  The elements in this list
#'			object contain a vector of labels that align the MM effects in a given term to disjoint, non-communicating groups.  For these data, there are no 
#'			groupings, so the single entry in 'group' associated to the 'mmcar' prior on number of sessions attended is filled with 1's equal to the number of
#'			of MM effects in each term.
#'   \item Omega. A list object containing an \code{S[g] x S[g]} CAR adjacency matrix used to model prior dependence among effects for each MM term under the 'mmcar'	
#'			prior, where \code{S[g]} are the number of effects for term \code{i}.  Here, there are two terms (2 and 4) under "mmcar", so \code{Omega} contains two elements.
#'   \item us. A list object where each element is a vector of true MM session effect values for a given MM term.  There are four MM terms in these data, so there are 4 elements, each of length \code{S[g]}.
#'   \item s. A vector of true cluster memberships for each of the \code{n} subjects.
#'   \item b.star. A list object of true cluster location values for each of \code{M = 10} clusters.  Each entry contains the \code{n.random = 3} location values for that 
#'		cluster.
#'   \item b. A list object true random effect coefficient values for each of \code{n} subjects.  Each entry contains the \code{n.random = 3} effect values for that subject.
#'   \item coefs. A numeric vector of length {n.fixed + 1} that contains the true fixed effects coefficient values for the intercept plus the columns in \code{X} (in order).
#'   \item tau.b. true values for the prior precisions of the base Gaussian distribution for each of \code{n.random = 3} subject effects.
#'   \item tau.u. An \code{4 x 1} vector of precision parameters associated to the Gaussian prior formulations for each MM term (under either 'mmcar' or 'mmi' covariance
#'			constructions).
#'   \item tau.e. true value for overall model error.
#'   \item option. A character vector equal to \code{c("mmi","mmcar","mmi","mmcar")} that presents the prior formulations chosen for the 4 MM terms; 
#'	 The "mmi" employs an independent Gaussian prior while the "mmcar" allows for adjacency based dependence.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name datsimmult
#' @usage datsimmult
#' @format A list object of 22 variables for 897 total observations on 299 subjects
NULL

 