###################################
## generic for non-session models
###################################
#' @include gpdpgrow.R
NULL

#' Bayesian instrinsic Gaussian Markov Random Field model for dependent time-indexed functions
#'
#' Estimates a collection of time-indexed functions under intrinsic Gaussian Markov random
#' field prior formulations where a Dirichlet process mixture allows sub-groupings of the functions to share the same 
#' iGMRF precision parameter.  The iGMRF formulation supports any number of additive precision terms,
#' expressing either or both of multiple trend and seasonality.
#'
#' @param y A multivariate continuous response, specified as an \emph{N x T} matrix, where \code{N}
#'     denotes the number of functions and \code{T}, the number of time points per function. Intermittent 
#'     missing-at-random values are allowed and will be estimated from the posterior
#'     predictive distribution.  Missing cells should be denoted with \code{NA}.
#' @param ipr An optional input vector of inclusion probabilities for each observation unit in the case
#'        the observed data were acquired through an informative sampling design, so that unbiased
#'        inference about the population requires adjustments to the observed sample.  Defaults to
#'        \code{ipr = rep(1,nrow(y))} indicating an iid sample.
#' @param q_order An integer vector of length \code{K} to select the order for each iGMRF precision term.
#'   e.g. If the first term is a RW2 and there is a second is a 3-month seasonality term, where the time
#'   points are indexed by month, then \code{q_order = c(2,3)}.  Defaults to \code{q_order = 2}
#' @param q_type A vector of length \code{K}, the number of iGMRF precision terms, with each entry 
#'    indicating whether the associated term is a trend (\code{"tr"}) or seasonality (\code{"sn"}) term.
#'    So all entries must be one of \code{c("tr","sn")}. Defaults to \code{q_type = "tr"}.
#' @param q_shape The value (in (0,infty)) for the shape hyperparameter for the Gamma base distribution for
#'    the iGMRF scale parameters, \code{kappa_star(k,m)}, where \code{k} denotes the term
#'    and \code{m}, the cluster.  Defaults to \code{q_shape = 0.3}.
#' @param q_rate The rate parameter of the Gamma base distribution on \code{kappa_star}. 
#'        Defaults to \code{q_rate = 0.0005}.
#' @param tau_shape The value (in (0,infty)) for the shape hyperparameter for the Gamma prior on the error
#'        precision parameter. Defaults to \code{tau_shape = 1.0}.
#' @param tau_rate The rate parameter of the Gamma prior distribution on \code{tau_e}. 
#'        Defaults to \code{tau_rate = 1}.
#' @param dp_shape The shape parameter for the Gamma prior on the DP concentration parameter, 
#'   \code{conc}. Defaults to \code{dp_shape = 1}.
#' @param dp_rate The rate parameter for the Gamma prior on the DP concentration parameter, 
#'   \code{conc}. Defaults to \code{dp_rate = 1}.
#' @param M_init Starting number of clusters of \code{nrow(y)} units to initialize sampler.
#'   Defaults to \code{M_init = nrow(y)}.
#' @param n.iter Total number of MCMC iterations.
#' @param n.burn Number of MCMC iterations to discard.  
#'   \code{gmrfdpgrow} will return \code{(n.iter - n.burn)} posterior samples.
#' @param n.thin Gap between successive sampling iterations to save.
#' @param progress A boolean value denoting whether to display a progress bar during model execution.
#'        Defaults to \code{progress = true}.
#' @param jitter A scalar double indicating amount of jitter to subract from the posterior 
#'        rate and shape hyperparameters of \code{tau_e} to stabilize computation.  
#'        Defaults to \code{jitter = 0.0}.
#' @param kappa_fast Boolean for whether to generate rate hyperparameter from full conditionals
#'                  versus joint Gaussian (on random effects, \code{bb}, given \code{kappa}.  The
#'                  former is faster, but numerically less stable. 
#'                  Defaults to \code{kappa_fast = FALSE}.
#' @return S3 \code{gmrfdpgrow} object, for which many methods are available to return and view results.  
#'        Generic functions applied to an object, \code{res} of class \code{gmrfdpgrow}, includes:
#'   \item{plot(res)}{ returns results plots, including fit functions versus data and allocation
#'                  of fitted functions into clusters}
#'  	\item{samples(res)}{ contains (\code{n.iter - n.burn}) posterior sampling iterations 
#'             for every model parameter}
#'	\item{resid(res)}{ contains the model residuals.}
#' @note The intended focus for this package are data composed of observed noisy functions (each of 
#'        length \code{T}) for a set of experimental units where the functions may express dependence
#'        among the experimental units
#' @keywords model
#' @seealso \code{\link{gmrfdpgrow}}
#' @examples 
#' {
#' library(growfunctions)
#' 
#' ## load the monthly employment count data for a collection of 
#' ## U.S. states from the Current 
#' ## Population Survey (cps)
#' data(cps)
#' ## subselect the columns of N x T, y, associated 
#' ## with the years 2008 - 2013
#' ## to examine the state level employment levels 
#' ## during the "great recession"
#' y_short   <- cps$y[,(cps$yr_label %in% c(2008:2013))]
#' 
#' ## Run the DP mixture of iGMRF's to estimate posterior 
#' ## distributions for model parameters
#' ## Under default RW2(kappa) = order 2 trend 
#' ## precision term
#' ## Run for 1500 iterations, with half as burn-in for a
#' ## more useful (converged) result.
#' res_gmrf            <- gmrfdpgrow(y = y_short, 
#'                                      n.iter = 40, 
#'                                      n.burn = 20, 
#'                                      n.thin = 1) 
#'                                      
#' ## 2 plots of estimated functions: 1. faceted by cluster and fit;
#' ## 2.  data for experimental units.
#' ## for a group of randomly-selected functions
#' fit_plots_gmrf      <- cluster_plot( object = res_gmrf, 
#'                                      units_name = "state", 
#'                                      units_label = cps$st, 
#'                                      single_unit = FALSE, 
#'                                      credible = TRUE )    
#' }
#' @aliases gmrfdpgrow
#' @aliases gmrfdpgrow.default
#' @author Terrance Savitsky \email{tds151@@gmail.com} Daniell toth \email{danielltoth@@yahoo.com}
#' @references 
#'	T. D. Savitsky and D. Toth (2014) Bayesian Non-parametric Models for Collections of Time-
#'     indexed Functions. submitted to: JRSS Series A (Statistics in Society).
#' @references
#'     T. D. Savitsky (2014) Bayesian Non-parametric Functional Mixture
#'     Estimation for Time-indexed data. submitted to: Annals of Applied Statistics.
#' @references
#'     T. D. Savitsky (2014) Bayesian Non-Parametric Mixture Estimation for Time-Indexed Functional
#'     Data for \code{R}. Submitted to: Journal of Statistical Software.     
#' @import reshape2 scales ggplot2 Rcpp RcppArmadillo Matrix spam mvtnorm     
#' @export 
gmrfdpgrow		<- function(y, ipr, q_order, q_type, q_shape, q_rate, tau_shape, tau_rate,
                        dp_shape, dp_rate, M_init, n.iter, n.burn, n.thin, progress, jitter, kappa_fast)
     UseMethod("gmrfdpgrow")

################################################
## default dispatch method for mm-session models
################################################
#' @export
gmrfdpgrow.default		<- function(y, ipr = NULL, q_order = 2, q_type = "tr", q_shape = 1.0, 
                            q_rate = 0.1, tau_shape = 1.0, tau_rate = 1.0, dp_shape = 1, 
                            dp_rate = 1, M_init = 2, 
                            n.iter = 10000, n.burn = 5000, n.thin = 5, progress = TRUE,
                            jitter = 0.0, kappa_fast = TRUE)
{ ## start function gmrfdpgrow.default
     ########################################################################
     ## check user inputs for correctness
     ########################################################################
     if( !(class(y) %in% c("data.frame","matrix")) )
     {
          stop("\nData, y, must be an N x T numeric data matrix.\n")
     } ## end check of inputted data.
     
     if( !all(q_type %in% c("tr","sn")) )
     {
          stop("\nPrecision matrix type for each term must be one of 'tr' for
               a trend matrix or 'sn' for a seasonal matrix.\n")
     } ## end check of allowed entries for q_type
     
     if( !all(q_order == floor(q_order)) ) ## ensure each input is an integer
     {
          stop("\nOrder of each precision matrix term must be specified as an integer.\n")
     } ## end check of allowable entries for q_order
     
     if( length(q_order) != length(q_type) )
     {
          stop("\nThe length of inputs 'q_type' and 'q_order' must be the same
               and equal to the number of desired iGMRF terms.\n")
     } ## end check to make sure q_order and q_type are of the same length equal to number of terms
     
     if( n.iter < n.burn )
     {
          stop("\nThe number of sampling iterations, n.iter, must be greater than the 
               numer of such iterations discarded as burn-in, n.burn.\n")
     } ## end check to make sure number of sampling iterations greater than number of burn-in
     #########################################################################
     ## set inclusion probability vector, ipr, equal to 1's if not input by user
     #########################################################################
     ## assumes data acquired in an iid sampling design
     if(is.null(ipr)){ipr <- rep(1,nrow(y))}
     
     #########################################################################
     ## compose trend and seasonal data objects
     #########################################################################
     ## render q_type to lower case
     q_type         <- tolower(q_type)
     q_order        <- as.numeric(q_order)
     ## if Q is seasonal, rank is (T - (seasonal period-1))
     ## if Q is trend, then rank is T - order
     q_order_adj    <- vector("numeric",length(q_order))
     T              <- ncol(y) # number of time observations per experimental unit
     N              <- nrow(y) # number of experimental units
     K_t            <- length(grep("tr",q_type)) ## number of trend precision terms
     K_s            <- length(grep("sn",q_type)) ## number of seasonal precision terms
     K              <- K_t + K_s
     ## each trend and seasonality component will have a unique order
     Q_k  <- matrix(0,T,T) ## T x T precision matrices for the K terms; e.g. Q_k = D_k - Omega_k
     D    <- matrix(0,K,T) ## each row will hold the T diagonal elements of Q_k
     C    <- vector("list",K) ## T x T C_k = D_k^-1 * Omega_k. (memo: diagonal element of Omega_k are 0)
     Omega_k <- matrix(0,T,T) ## T x T adjacency matrix, Omega_k = D_k - Q_k
     for( k in 1:K )
     { ## determine if GMRF precision matrix is trend or seasonal
          if(q_type[k] == "tr")
          {
               Q_k         <- as.matrix(precmat.RWn(T,order=q_order[k]))
               ## populate q_order_adj[k] based on whether Q is trend or seasonal
               q_order_adj[k]   <- q_order[k]
               
          }else{ ## q_type[k] == "sn"
               if(q_type[k] == "sn")
               {
                    Q_k         <- as.matrix(precmat.season(T, season=q_order[k]))
                    ## populate q_order_adj[k] based on whether Q is trend or seasonal
                    q_order_adj[k]      <- (q_order[k] - 1) ## rank of Q for drawing of kappa_star
               }
               
          }
          ## compute derived (sparse) matrices for C++
          D[k,]          <- diag(Q_k)
          Omega_k        <- diag(D[k,]) - Q_k  ## Q_k = diag(D[k,]) - Omega_k
          C[[k]]         <- Omega_k / D[k,]
          
          ## convert to dcgmatrices (under Matrix package) required for sparse representation
          ## as sp_mats in RcppArmadillo
          ## C[[k]]               <- as(C[[k]], "dgCMatrix") 
     } ## end loop over constructing K iGMRF precision matrices
     
     #########################################################################
     ## convert input R NA's in y to an atomic integer value (-9) for sampling in C++
     #########################################################################
     y[is.na(y)] <- -9 ## this value will trigger an update step in c++ for missing y values in each row.
     
     ## set progress (boolean on whether to display sampling progress bar) to lower for c++ input
     progress       <- as.integer(progress)
     kappa_fast     <- as.integer(kappa_fast)
     
     ################################################################
     ## conduct posterior sampling and capture results
     ################################################################
     ## set q_order_2 to 0 to stabilize computation
     res       <- gmrfdpPost(y, ipr, C, D, q_order_adj, q_type, n.iter, n.burn, n.thin, 
                             M_init, q_shape, q_rate, tau_shape, tau_rate,
                             dp_shape, dp_rate, progress, 
                             jitter, kappa_fast)
     cat(paste("Your additive set of iGMRF precision terms includes type = ", 
               c(q_type)," and order = ", c(q_order), "\n", sep=""))
     ## add in sn_order to construct test set covariance matrix used for prediction
     res$q_order                   <- q_order
     res$optpartition$y[y == -9]   <- NA ## convert missing values back from -9 to NA
     ## Add column names
     ## kappa, k is fast moving
     colnames(res$Kappa)           <- paste("kappa[",rep(1:K,times=N),",",rep(1:N,each=K),"]",sep="")
     ## bb
     colnames(res$bb)              <- paste("bb[",rep(1:N,times=T),",",rep(1:T,each=N),"]",sep="")
     ## Tau_e
     colnames(res$Tau_e)           <- "tau_e"
     ## assign column labels to each function, f_1, ... , f_K, where bb = f_1 + ... + f_K
     for( k in 1:K )  
     {
          colnames(res$f[[k]])     <- paste("f_",k,"[",rep(1:N,times=T),",",rep(1:T,each=N),"]",sep="")
     } ## end loop over K additive function terms
     
     ## need this object in cluster_plot() for plotting data vs. fitted functions
     res$residuals       <- colMeans(res$Residuals)
     ## save type and order identifiers for each term to res for use in predict_functions()
     res$q_order         <- q_order
     res$q_type          <- q_type
     class(res)          <- c("gmrfdpgrow")
     return(res)
     
} ## end function gpdpgrow()

#####################################################
## .Call statements to C++ functions
#####################################################
#' Run a Bayesian functional data model under an instrinsic GMRF prior whose precision parameters
#' employ a DP prior
#'
#' An internal function to \code{\link{gmrfdpgrow}}
#'
#' @export
#' @aliases gmrfdpPost
#' @param y An \emph{N x T} matrix of N observations of \emph{T x 1} functions
#' @param ipr An optional input vector of inclusion probabilities for each observation unit in the case
#'        the observed data were acquired through an informative sampling design, so that unbiased
#'        inference about the population requires adjustments to the observed sample.  Defaults to
#'        \code{ipr = rep(1,nrow(y))} indicating an iid sample.
#' @param C A list object of length, \code{K}, the number of iGMRF precision terms.
#'          Each entry contains a \emph{T x T} normalized adjacency matrix.  The diagonal entries are
#'          \code{0} and row \code{i} contains the weight for each entry {!=i} divided by the sum
#'          of the weights.
#' @param D A \emph{K x T} matrix, where \code{K} denotes the number of iGMRF terms.
#'          Row \code{k} contains the \code{T} elements of the diagonal of the term-\code{k}
#'          precision matrix, \code{Q_k}.
#'          Will increase with order and be equal, except for boundary corrections.     
#' @param q_order An integer vector where each entry contains the order of the associated \code{K}
#'        iGMRF precision terms
#'        matrix of Euclidean distances associated to each seasonal covariance term.
#' @param q_type A vector of length \code{K}, the number of iGMRF precision terms, with each entry 
#'        indicating whether the associated term is a trend (\code{"tr"}) or 
#'        seasonality (\code{"sn"}) term.
#' @param tau_shape The value (in (0,infty)) for the shape hyperparameter for the Gamma prior on the error
#'        precision parameter. Defaults to \code{tau_shape = 1.0}.
#' @param tau_rate The rate parameter of the Gamma prior distribution on \code{tau_e}. 
#'        Defaults to \code{tau_rate = 1}.
#' @param n.iter The number of MCMC sampling iterations
#' @param n.burn The number of warm-up iterations to discard
#' @param n.thin The interval or step size of post-burn-in samples to return
#' @param M_init Starting value of number of clusters for sampling cluster assignments. 
#' @param q_shape The shape parameter of the Gamma base distribution for the \code{kappa_star}
#'        locations used to sample the DP prior on the \code{P} GP covariance parameters, 
#'        \code{kappa}, for each experimental unit.
#' @param q_rate The rate parameter of the Gamma base distribution for the \code{kappa_star}
#'        locations used to sample the DP prior on the \code{P} GP covariance parameters, 
#'        \code{kappa}, for each experimental unit.  
#' @param dp_shape The shape parameter for the \eqn{\Gamma} prior on the DP concentration parameter.  
#'     The rate parameter is set of \code{1}.
#' @param dp_rate The rate parameter for the \eqn{\Gamma} prior on the DP concentration parameter. 
#'        Default value is \code{1}.
#' @param progress An indicator in \code{{0,1}} denoting whether to display a progress bar during model execution.
#'        \code{progress = 1} displays a progress bar. Defaults to \code{progress = 1}.
#' @param jitter A scalar double indicating amount of jitter to subract from the posterior 
#'        rate and shape hyperparameters of \code{tau_e} to stabilize computation.  
#'        Defaults to \code{jitter = 0.0}.
#' @param kappa_fast Boolean for whether to generate rate hyperparameter from full conditionals
#'                  versus joint Gaussian (on random effects, \code{bb}, given \code{kappa}.  The
#'                  former is faster, but numerically less stable. 
#'                  Defaults to \code{kappa_fast = FALSE}.
#' @return res A list object containing MCMC runs for all model parameters.
#' @seealso \code{\link{gpdpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{gmrfdpgrow}}
gmrfdpPost = function (y, ipr, C, D, q_order, q_type,
                       n.iter, n.burn, n.thin, M_init, 
                       q_shape, q_rate, tau_shape, tau_rate, dp_shape, dp_rate, progress, 
                       jitter, kappa_fast) 
{ 
     y    <- as.matrix(y)
     stopifnot(nrow(D) == length(C))
     stopifnot(ncol(y) == nrow(C[[1]]))
     
     ## adjust modeling order to rank of precision matrices.
     ## for RW terms, that is the order.  for season terms that is the (order+1).
     K          <- length(q_type) ## number of GMRF terms
     mod_order  <- rep(0,K)
     for( k in 1:K )
     {
          if(q_type[k] == "sn")
          {
               mod_order[k] <- (q_order[k]-1)
          }else{ ## q_type[k] == "tr"
               mod_order[k] <- q_order[k]
          }
     } ## end loop over K GMRF terms to set order vector for modeling
     res <- .Call("IGMRFDPMIX", y, C, D, mod_order, 
                  n.iter, n.burn, n.thin, M_init, 
                  q_shape, q_rate, tau_shape, tau_rate, dp_shape, dp_rate, progress, 
                  jitter, kappa_fast, ipr,
                  PACKAGE = "growfunctions")
} ## end function gmrfdpPost

#####################################################
## function to predict iGMRF at future time values
#####################################################
#' Use the model-estimated iGMRF precision parameters from gmrfdpgrow() to predict the iGMRF function at
#' future time points.  Inputs the \code{gmrfdpgrow} object of estimated parameters.
#'
#' A companion function to \code{\link{gmrfdpgrow}}
#'
#' @export 
#' @param object Object of class \code{gmrfdpgrow} returned from model run of \code{gmrfdpgrow()}
#' @param J Scalar denoting number of draws to take from posterior predictive for each unit.
#'          Defaults to \code{J = 500}.
#' @param T_test The number of equally-spaced time points to predict the iGMRF functions ahead of 
#'               of the functions estimated at \code{T_train} time points.       
#' @param ... further arguments passed to or from other methods.
#' @return out  A list object containing containing two matrices; the first is a P x (N*T)
#'                  matrix of predicted function values for each of P sampled iterations.  N is 
#'                  slow index and denotes the number of experimental units.  The second matrix is
#'                  an N x T average over the P sampled draws, composed in Rao-Blackwellized fashion.
#' @export predict_functions gmrfdpgrow
#' @aliases predict_functions.gmrfdpgrow
#' @seealso \code{\link{gmrfdpgrow}}
#' @examples 
#' \dontrun{
#' library(growfunctions)
#' data(cps)
#' y_short   <- cps$y[,(cps$yr_label %in% c(2010:2013))]
#' t_train   <- ncol(y_short)
#' N         <- nrow(y_short)
#' t_test    <- 4
#'  
#' ## Model Runs
#'
#' res_gmrf            = gmrfdpgrow(y = y_short, 
#'                                 q_order = c(2,4), 
#'                                 q_type = c("tr","sn"), 
#'                                 n.iter = 100, 
#'                                 n.burn = 50, 
#'                                 n.thin = 1) 
#' ## Prediction Model Runs
#' T_test             <- 4
#'
#' pred_gmrf          <- predict_functions( object = res_gmrf,
#'                                      J = 1000, 
#'                                      T_test = T_test )
#'
#' ## plot estimated and predicted functions
#' plot_gmrf       <- predict_plot(object = pred_gmrf, 
#'                                units_label = cps$st, 
#'                                single_unit = TRUE, 
#'                                credible = FALSE)
#' }
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as a companion function for \code{\link{gmrfdpgrow}} for prediction
predict_functions.gmrfdpgrow = function (object, J=500, T_test,...) 
{
     ## Extract dimensions 
     y                   <- object$optpartition$y
     T_train             <- ncol(y)
     T_tot               <- T_train + T_test
     ## q_type and q_order tie the matrix iGMRF T_tot x T_tot precision matrices, R[[k]], used for
     ## prediction to those used for estimation of functions and to the estimated precision parameters
     ## kappa[k,i]
     q_type              <- object$q_type ## "sn" or "tr" identifiers for each of K GMRF terms
     q_order             <- object$q_order ## order identifier for each of K GMRF terms
     K_t                 <- length(grep("tr",q_type)) ## number of trend precision terms
     K_s                 <- length(grep("sn",q_type)) ## number of seasonal precision terms
     K                   <- K_t + K_s ## total number of iGMRF terms
     
     ## Compose list of T_tot x T_tot precision matrices, R_k, where Q_ki = kappa_ki * R_k
     ## for precision parameter, kappa_ki, for term k and unit i.
     R_k  <- matrix(0,T_tot,T_tot) ## T_tot x T_tot precision matrices for the K terms; 
                               ## e.g. R_k = D_k - Omega_k
     R    <- vector("list",K) ## list of K, R_k matrices
     for( k in 1:K )
     { ## determine if GMRF precision matrix is trend or seasonal
          if(q_type[k] == "tr")
          {
               R_k         <- as.matrix(precmat.RWn(T_tot,order=q_order[k]))
               
          }else{ ## q_type[k] == "sn"
               if(q_type[k] == "sn")
               {
                    R_k         <- as.matrix(precmat.season(T_tot, season=q_order[k]))  
               }
               
          }
          ## using dense matrices, for now
          R[[k]]         <- R_k
          
     } ## end loop over constructing K iGMRF precision matrices 
     
     out <- .Call("predict_gmrf_bb", object, R, J, PACKAGE = "growfunctions")
     
     class(out)                    <- c("predict_functions.gmrfdpgrow")
     return(out)
} ## end function predict_functions()


####################################
## accessor methods
####################################

#' Produce samples of MCMC output
#'
#' provides posterior sampled values for every model parameter of a
#' \code{gmrfdpgrow} object
#'
#' @param object A \code{gmrfdpgrow} object
#' @param ... Ignored
#' @export 
#' @aliases samples.gmrfdpgrow
samples.gmrfdpgrow <- function(object,...)
{
     
     res            <- list(Kappa = object$Kappa, bb = object$bb, f = object$f,
                            Tau_e = object$Tau_e,
                            M = object$M, Conc = object$Conc)
     class(res) 	<- "samples.gmrfdpgrow"
     return(res)
}
