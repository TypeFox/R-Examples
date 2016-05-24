#' @title
#' Light Cone Reconstruction of States - Predictive State Estimation From Spatio-Temporal Data
#' @aliases LICORS
#' @docType package
#' @name LICORS-package
#' @author
#' Georg M. Goerg \email{gmg@@stat.cmu.edu}
#' @description
#' A package for predictive state estimation from spatio-temporal data. The main 
#' function is \code{\link{mixed_LICORS}}, which implements an EM algorithm
#' for predictive state recovery (see References).
#' 
#' This is an early release: some function names and arguments might/will
#' (slightly) change in the future, so regularly check with new package updates.
#' @keywords package
#' @references 
#' Goerg and Shalizi (2013), JMLR W\&CP 31:289-297. Also available at 
#' \url{arxiv.org/abs/1211.3760}.
#' 
#' Goerg and Shalizi (2012). Available at \url{arxiv.org/abs/1206.2398}.
#' @seealso The main function in this package: \code{\link{mixed_LICORS}}
#' @import zoo RColorBrewer FNN mvtnorm locfit Matrix
#' 
#' @section Details on Methodology - Predictive State Model for Spatio-temporal Processes:
#' 
#' \emph{For details and additional references please consult 
#' Goerg and Shalizi (2012, 2013).}
#' 
#' Let \eqn{\mathcal{D} = \lbrace X(\mathbf{r}, t) \mid \mathbf{r} \in \mathbf{S}, t = 1, \ldots, T \rbrace = (X_1, \ldots, X_{\tilde{N}})} 
#' be a sample from a spatio-temporal process, observed over an 
#' \eqn{N}-dimensional spatial grid \eqn{\mathbf{S}} and for \eqn{T} time steps.
#' We want to find a model that is optimal for forecasting a new 
#' \eqn{X(\mathbf{s}, u)} given the data \eqn{\mathcal{D}}.  To do this we need to
#' know
#' \deqn{
#' P(X(\mathbf{s}, u) \mid \mathcal{D})
#' }
#' 
#' In general this is too complicated/time-intensive since \eqn{\mathcal{D}} is
#' very high-dimensional.  But we know that in any physical system, information
#' can only propagate at a finite speed, and thus we can restrict the search for
#' optimal predictors to a subset 
#' \eqn{\ell^{-}(\mathbf{r}, t) \subset \mathcal{D}}; this is the 
#' \strong{past light cone (PLC)} at \eqn{(\mathbf{r}, t)}. 
#' 
#' There exists a mapping \eqn{\epsilon: \ell^{-} \rightarrow \mathcal{S}}, where
#' \eqn{\mathcal{S} = \lbrace s_1, \ldots, s_K \rbrace} is the predictive state 
#' space.  This mapping is such that
#' \deqn{
#' P(X_i \mid \ell^{-}_i) = P(X_i \mid  s_j),
#' }
#' where \eqn{s_j = \epsilon(\ell^{-}_i)} is the predictive state of PLC \eqn{i}.
#' Furthermore, the future is independent of the past given the predictive state:
#' \deqn{
#' P(X_i \mid \ell^{-}_i, s_j) = P(X_i \mid  s_j) .
#' }
#'
#' 
#' The likelihood of the joint process factorizes as a product of
#' predictive conditional distributions
#' \deqn{
#' P(X_1, \ldots, X_N ) \propto \prod_{i=1}^{N} P(X_i \mid \ell^{-}_i) 
#' =  \prod_{i=1}^{N} P(X_i \mid \epsilon(\ell^{-}_i)).
#' }
#' 
#' Since \eqn{s_j} is unknown this can be seen as the complete data
#' likelihood of a nonparametric finite mixture model over predictive states:
#' \deqn{
#' P(X_1, \ldots, X_N ) \propto \prod_{i=1}^{N} \sum_{j=1}^{K} 
#'              \mathbf{1}(\epsilon(\ell^{-}_i) = s_j) \times P(X_i \mid s_j).
#' }
#' 
#' This predictive state model is a provably optimal finite 
#' mixture model, where the ``parameter'' \eqn{\epsilon} is chosen to provide
#' optimal forecasts.
#' 
#' The \pkg{LICORS} R package implements methods to estimate this optimal mapping 
#' \eqn{\epsilon}.
#' 
#' @section Acronyms and common function arguments:
#' 
#' The R package uses a lot of acronyms and terminology from the References, 
#' which are provided here for the sake of clarity/easier function 
#' navigation:
#' \describe{
#'  \item{LCs}{light cones}
#'  \item{PLC}{past light cone; notation: \eqn{\ell^{-}}}
#'  \item{FLC}{future light cone; notation: \eqn{\ell^{+}}}
#'  \item{LICORS}{LIght COne Reconstruction of States}
#' }
#' 
#' Many functions use these acryonyms as part of their name.  Function arguments
#' that repeat over and over again are:
#' \describe{
#'  \item{\code{weight.matrix}}{an \eqn{N \times K} matrix, where \eqn{N} are 
#'                              the samples and \eqn{K} are the states. That is, 
#'                              each row contains a vector of length \eqn{K} that
#'                              adds up to one (the mixture weights).}
#'  \item{\code{states}}{a vector of length \eqn{N} with entry \eqn{i} being
#'                             the label \eqn{k = 1, \ldots, K} of PLC \eqn{i}}                    
#' } 
#' 
#' @examples
#' \dontrun{
#' # setup the light cone geometry
#' LC_geom = setup_LC_geometry(speed=1, horizon=list(PLC = 2, FLC = 0), shape ="cone")
#' # load the field
#' data(contCA00)
#' # get LC configurations from field
#' contCA_LCs = data2LCs(contCA00$observed, LC.coordinates = LC_geom$coordinates)
#' # run mixed LICORS
#' 
#' mod = mixed_LICORS(contCA_LCs, num.states_start = 10, initialization = "KmeansPLC", max_iter = 20)
#' 
#' plot(mod)
#' }
#' 
#' 
NULL
