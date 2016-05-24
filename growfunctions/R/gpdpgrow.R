###################################
## generic for non-session models
###################################
#' @include gmrfdpgrow.R
NULL

#' Bayesian non-parametric dependent Gaussian process model for time-indexed functional data
#'
#' Estimates a collection of time-indexed functions with Gaussian process (GP) formulations 
#' where a Dirichlet process mixture allows sub-groupings of the functions to share the same 
#' GP covariance parameters.  The GP formulation supports any number of additive GP covariance terms,
#' expressing either or both of multiple trend and seasonality.
#'
#' @param y A multivariate continuous response, specified as an \emph{N x T} matrix, where \code{N}
#'     denotes the number of functions and \code{T}, the number of time points per function. Intermittent 
#'     missing-at-random values are allowed and will be estimated from the posterior
#'     predictive distribution.  Missing cells should be denoted with \code{NA}.  The sampling of
#'     missing values requires co-sampling the functions, \code{bb}, while these functions are 
#'     marginalized out for sampling under the case of no missing data.  So the sampler will run
#'     more slowly in the case of intermittent missingness.
#' @param ipr An optional input vector of inclusion probabilities for each observation unit in the case
#'        the observed data were acquired through an informative sampling design, so that unbiased
#'        inference about the population requires adjustments to the observed sample.  Defaults to
#'        \code{ipr = rep(1,nrow(y))} indicating an iid sample.
#' @param time_points Inputs a vector of common time points at which the collections of functions were
#'        observed (with the possibility of intermittent missingness).  The length of \code{time_points}
#'        should be equal to the number of columns in the data matrix, \code{y}.  Defaults to 
#'        \code{time_points = 1:ncol(y)}.
#' @param gp_cov A vector of length \code{L} to select the covariance function for each of
#'   \code{L} terms. Allowed inputs are \code{c("rq","se","sn")}, where "rq" denotes the 
#'   rational quadratic covariance function, "se", the square exponential, and "sn", seasonality.
#'   The seasonality covariance is a quasi-periodic multiplicative combination of a fixed length-
#'   scale periodic covariance kernel (where the period length is specified in 'sn_order') and
#'   a squared exponential kernel (of varying length-scale) to allow the periodicity in the seasonal
#'   term to evolve.  e.g. If 4 terms are wanted - 2 trend terms with "se" and 2 seasonality terms, the input would be
#'   \code{gp_cov = c("se","se","sn","sn")}. Defaults to \code{gp_cov = "rq"}.
#' @param sn_order A vector of length \code{L_s}, the number of terms in \code{gp_cov} where \code{"sn"}
#'   is selected, that denotes the seasonality order for each term; e.g. if the two "sn" terms above 
#'   are for 3 and 12 month seasonality, respectively, for monthly data, then 
#'   \code{sn_order = c(3,12)}. Defaults to \code{sn_order = NULL}.
#' @param jitter A scalar numerical value added to the diagonal elements of the T x T GP covariance 
#' matrix to stabilize computation.  Defaults to \code{jitter = 0.01}.
#' @param gp_shape The shape parameter of the Gamma base distribution for the DP prior on
#'   the P x N matrix of GP covariance parameters (where P 
#'   denotes the number of parameters for each of the N experimental units).
#'   Defaults to \code{gp_shape = 1}
#' @param gp_rate The rate parameter of the Gamma base distribution on GP covariance parameters. 
#'   Defaults to \code{gp_rate = 1}.
#' @param noise_shape The shape parameter of the Gamma base distribution on \code{tau_e}, the
#'   model noise precision parameter. Defaults to \code{noise_shape = 3}.
#' @param noise_rate The rate parameter of the Gamma base distribution on \code{tau_e}, the model
#'   noise precision parameter.  Defaults to \code{noise_rate = 1}.
#' @param dp_shape The shape parameter for the Gamma prior on the DP concentration parameter, 
#'   \code{conc}. Defaults to \code{dp_shape = 1}.
#' @param dp_rate The rate parameter for the Gamma prior on the DP concentration parameter, 
#'   \code{conc}. Defaults to \code{dp_rate = 1}.
#' @param M_init Starting number of clusters of \code{nrow(y)} units to initialize sampler.
#'   Defaults to \code{M_init = nrow(y)}.
#' @param lower The lower end of the range to be used in conditionally sampling the GP covariance 
#'   parameters (\code{kappa,tau_e}) in the slice sampler.  Defaults to \code{lower = 0}.
#' @param upper The upper end of the range to be used in conditionally sampling the GP covariance
#'   parameters (\code{kappa,tau_e}) in the slice sampler.  Defaults to \code{upper = 1e10}.
#' @param sub_size Integer vector whose length, \code{n}, equals the number of progressively coarser
#'   GP covariance matrices to use for tempered sampling steps in an alternative space to sample the
#'   GP covariance parameters.  Each entry denotes the number of sub-sample time points in \code{T}
#'   to draw under a latin hypercube design from the \emph{N x T} data matrix, \code{y} to employ
#'   for that distribution; for example, suppose \code{T = 300}.  \code{sub_size = c(100,50)}, 
#'   would randomly select first \code{100} and then \code{50} time points from \code{y} to use
#'   for each of the \code{n = 2} distributions.  Defaults to 
#'   \code{sub_size = c(floor(0.25*T),floor(0.1*T))}.
#' @param w_star Integer value denoting the number of cluster locations to sample ahead of 
#'   observations in the auxiliary Gibbs sampler used to sample the number of clusters
#'   and associated cluster assignments.  A higher value reduces samplin auto-correlation, but
#'   increases computational burden.  Defaults to \code{w_star = 2}.
#' @param w Numeric value denoting the step width used to construct the interval from
#'   which to draw a sample for each GP covariance parameter in the slice sampler.  This
#'   value is adaptively updated in the sampler tuning stage for each parameter to be equal
#'   to the difference in the 0.95 and 0.05 sample quantiles for each of 5 block updates.
#'   Defaults to \code{w = 1.5}.
#' @param n.iter Total number of MCMC iterations.  
#' @param n.burn Number of MCMC iterations to discard.  
#'   \code{gpdpgrow} will return \code{(n.iter - n.burn)} posterior samples.
#' @param n.thin Gap between successive sampling iterations to save.
#' @param n.tune Number of iterations (before ergodic chain instantiated) to adapt \code{w}, separately,
#'   for each covariance term, \code{p = 1,...,P}.  Sets each {w_p} to lie in the 90 percent credible 
#'   interval computed from the tuning sample (that is divided into 5 blocks so that \code{w_p} is 
#'   successively updated in each block of runs).
#' @param progress A boolean value denoting whether to display a progress bar during model execution.
#'        Defaults to \code{progress = TRUE}
#' @param b_move A boolean value denoting whether to sample the GP function, \code{bb}, in \emph{T x 1}
#'        Gibbs steps \code{b_move = TRUE} or through elliptical slice sampling.  
#'        Defaults to \code{b_move = TRUE}. Only used in the case there is any intermittent missingness;
#'        otherwise \code{bb} is marginalized out of the sampler and post-sampled from the its predictive
#'        distribution.
#' @param cluster A boolean value denoting whether to employ DP mix model over set of GP functions or
#'        to just use GP model with no clustering of covariance function parameters.  
#'        Defaults to \code{cluster = TRUE}
#' @param s An \emph{N x 1} integer vector that inputs a fixed clustering, rather than sampling it.  
#'        Defaults to \code{s = NULL}
#' @return S3 \code{gpdpgrow} object, for which many methods are available to return and view results.  Generic functions applied
#'	to an object, \code{res} of class \code{gpdpgrow}, includes:
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
#' ## subselect the columns of N x T, y, associated with 
#' ## the years 2011 - 2013
#' ## to examine the state level employment 
#' ## levels during the "great recession"
#' y_short   <- cps$y[,(cps$yr_label %in% c(2011:2013))]
#' 
#' ## uses default setting of a single "rational quadratic" covariance
#' ## run for 500 iterations, with half discarded as burn-in to 
#' ## obtain a more useful result.
#' res_gp               <- gpdpgrow(y = y_short, 
#'                                  n.iter = 4, 
#'                                  n.burn = 1, 
#'                                  n.thin = 1, 
#'                                  n.tune = 0)  
#'
#' ## Two plots of estimated functions, 
#' ## 1. faceted by cluster 
#' ## 2. fitted functions vs noisy observations
#' ## first plot will plot estimated denoised function, 
#' ## bb_i, for a single (randomly-selected) "state"
#' fit_plots_gp        <- cluster_plot( object = res_gp,  
#'                            units_name = "state", 
#'                            units_label = cps$st, 
#'                            single_unit = TRUE, 
#'                            credible = TRUE )
#' ## second plot will randomly select 6 states 
#' ## and plot their estimated denoised functions, bb_i.
#' ## with setting "single_unit = FALSE". 
#' ## (Option "num_plot" may be set to plot 
#' ## any integer number of 
#' ## randomly-selected units.)
#' fit_plots_gp        <- cluster_plot( object = res_gp,  
#'                                      units_name = "state", 
#'                                      units_label = cps$st, 
#'                                      single_unit = FALSE, 
#'                                      credible = TRUE )
#'
#' }
#' @aliases gpdpgrow
#' @aliases gpdpgrow.default
#' @author Terrance Savitsky \email{tds151@@gmail.com} Daniell Toth \email{danielltoth@@yahoo.com}
#' @references 
#'	T. D. Savitsky and D. Toth (2014) Bayesian Non-parametric Models for Collections of Time-
#'     indexed Functions. submitted to: JRSS Series A (Statistics in Society).
#' @references
#'     T. D. Savitsky (2014) Bayesian Non-parametric Functional Mixture
#'     Estimation for Time-indexed data. submitted to: Annals of Applied Statistics.
#' @references
#'     T. D. Savitsky (2014) Bayesian Non-Parametric Mixture Estimation for Time-Indexed Functional
#'     Data for \code{R}. Submitted to: Journal of Statistical Software. 
#' @export 
gpdpgrow			<- function(y, ipr, time_points, gp_cov, sn_order, jitter, gp_shape, gp_rate, 
                                   noise_shape, noise_rate, dp_shape, dp_rate, 
                                   M_init, lower, upper, sub_size, w_star, w, n.iter, n.burn, 
                                   n.thin, n.tune, progress, b_move, cluster,s)
     UseMethod("gpdpgrow")

################################################
## default dispatch method for mm-session models
################################################
#' @export
gpdpgrow.default		<- function(y, ipr = NULL, time_points = NULL, gp_cov = c("rq"), sn_order = NULL, 
                            jitter = 0.001, gp_shape = 0.5, 
                            gp_rate = 0.5, noise_shape = 2, noise_rate = 0.1,
                            dp_shape = 1, dp_rate = 1, M_init = 5, lower = 0, 
                            upper = 1e3, sub_size = c(floor(0.35*ncol(y)),floor(0.18*ncol(y))),
                            w_star = 2, w = 1.5, n.iter = 8000, n.burn = 4000, n.thin = 5, 
                            n.tune = 1500, progress = TRUE, b_move = TRUE, cluster = TRUE,
                            s = NULL)
{ ## start function gpdpgrow.default
     ########################################################################
     ## check user inputs for correctness
     ########################################################################
     if( is.null(y) )
     {
          stop("\nYo, you gotta input an N x T data matrix to get this
               party started.\n")
     }else{
          if( !(class(y) %in% c("data.frame","matrix")) )
          {
               stop("\nData, y, must be an N x T numeric data matrix.\n")
          } ## end check of class for inputted data.
     } ## end check of inputted data.
     
     
     if( !all(gp_cov %in% c("rq","se","sn")) )
     {
          stop("\nThe covariance kernel for each term must be specified as either,
               'se' (squared exponential), 'rq' (rational quadratic),
               or 'sn' (seasonal, quasi-periodic).\n")
     } ## end check of allowed entries for q_type
     
     ## ensure that 'sn' in gp_cov triggers entry in sn_order
     if( (is.null(sn_order)) )
     {
          if( "sn" %in% gp_cov ) 
          {
               stop("\nIf choose 'sn' (seasonal covariance) in 'gp_cov', must also specify
                    an integer for 'sn_order' for each 'sn' entry in 'gp_cov' to denote
                    the integer order associated to the term.\n")
          }
          
     }else{ ## !is.null(sn_order)
          if( length(sn_order) != length(grep("sn",gp_cov)) )
          {
               stop("\nThe length of 'sn_order' (to input integers for each seasonal term)
                    must be equal to the number of 'sn' entries in 'gp_cov'.\n")
          } ## end condition that length(sn_order) must equal number of sn entries in gp_cov
     
          if( any(sn_order != floor(sn_order)) ) ## ensure each input is an integer
          {
               stop("\nOrder of each seasonal term entered in 'sn_order' must be an integer.\n")
          } ## end to ensure entries in sn_order are integers
     } ## end check of allowable entries for sn_order
     
     
     if( length(sub_size) > 2 )
     {
          warning("\nWhile it is permitted to use more than 2 coarse terms, computational performance
                  in terms of time per effective sample wouldn't expect to improve.\n")
     }else{
          if( any(sub_size != floor(sub_size)) )
          {
               stop("\nInputs to 'sub_size' must be integers since they subselect time
                    points to build progressively coarser approximations to 
                    the full covariance matrix.\n")
          }
     } ## end conditions to ensure sub_size is entered properly
     sub_size = sort(sub_size, decreasing = TRUE) ## sort in decreasing order
     
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
     ## create boolean indicating whether clustering is user input
     if(is.null(s))
     {
          fix = FALSE
     }else{
          fix = TRUE
     }
     ## render gp_cov to lower case
     gp_cov         <- tolower(gp_cov)
     T              <- ncol(y) # number of time observations per experimental unit
     L_t            <- length(grep("rq",gp_cov)) + length(grep("se",gp_cov)) 
     L_s            <- length(grep("sn",gp_cov))
     L              <- L_t + L_s
     ## there will always be a trend component - even if only seasonality selected
     ## the multiplicative trend component of seasonality allows the seasonality to vary
     ## over the time-indexed observations
     if( is.null(time_points) )
     {
          time_points    <- 1:T
     }
     
     Omega_t        <- sapply(1:T,function(i){
                              sapply(1:T,function(j){
                              d = (time_points[i]-time_points[j])^2
                              })
                         })
     if(L_s > 0)
     {
          Omega_s        <- vector("list",L_s)
          for( i in 1:L_s )
          {
               Omega_s[[i]]           <- sapply(1:T,function(j){
                                             sapply(1:T,function(k){
                                             d = 2*( sin(pi*(time_points[j]-time_points[k])/sn_order[i]) )^2 
                                             })
                                        })
          } ## end loop i over L_s seasonality terms
     }else{ ## no seasonality - compose dummy Omega_s for input to c++.  won't be used.
          Omega_s             <- vector("list",1)
          Omega_s[[1]]        <- Omega_t
     } ## end conditional statement on seasonality to construct Euclidean distance matrices 
     
     #########################################################################
     ## convert input R data objects to equivalents used with C++
     #########################################################################
     ## indexing GP covariance function used for each of L terms
     gp_mod                   <- vector("numeric",L)
     gp_mod[gp_cov == "rq"]   <- 1
     gp_mod[gp_cov == "se"]   <- 2
     gp_mod[gp_cov == "sn"]   <- 3
     
     ## set the maximum number of iterations to widen slice interval for sampling
     ## covariance parameters
     n_slice_iter             <- 1e3 ## large
     
     ## conduct stratified sample of columns of y for each i in (1,...n) tempered steps
     ## ensure columns being sampled contain observations (versus test columns)
     sel_obs   <- (apply(y,2,function(x){all(is.na(x))}) == FALSE) ## T/F of length ncol(y)
     y_obs     <- y[,sel_obs]
     T_obs     <- length(sel_obs[sel_obs == TRUE]) ## number of columns from which to conduct sampling
     if( T_obs < max(sub_size) )
          stop(paste("The number of columns in y with one or more observations ",T_obs," is less than
                      the max(sub_size) ", max(sub_size), sep=""))
     cols_obs  <- (1:ncol(y))[sel_obs] ## the columns eligible for sampling
     n         <- length(sub_size) ## defaults to n = 2
     y_index   <- vector("list",n)  # the indices sampled for each sub_size
     ## define strata from which to draw the n stratified samples   
     length_stratum    <- min(T_obs,15)
     n_stratum         <- floor( T_obs / length_stratum )
     if(n_stratum == 0){n_stratum = 1}
     if(n_stratum == 1) ## just an SRS
     {
          for( i in 1:n ) 
          {
               ## the "-1" accounts for the C++ indexing that starts at 0 for
               ## selecting columns of y, Omega_t and Omega_s[[1:L_s]]
               ## y_index[[i]] is the index of time_points (columns), not the response values
               y_index[[i]] <- sort(sample(cols_obs,sub_size[i],replace = FALSE)) - 1
          }
     }else{ ## n_stratum >= 2
          for( i in 1:n )
          {
               ## re-set number of strata to a lower value if number sampled per stratum == 1
               size_ih             <- floor(sub_size[i]/n_stratum)
               if( size_ih == 1 ) ## reduce number of stratum if only 1 obs per stratum (except last)
               {
                    size_ih_i                     <- 3
               }else{ ## more than one value sampled per stratum
                    size_ih_i                     <- size_ih 
               } ## end conditional statement on whether there is only one sampled value per stratum
               n_stratum_i                   <- floor( sub_size[i]/size_ih_i )
               length_stratum_i              <- floor( T_obs / n_stratum_i )
               cols_sel_i                    <- vector("list", n_stratum_i)
               cut_stratum_i                 <- 1:(n_stratum_i-1)*length_stratum_i 
               
               set_ih              <- cols_obs[cols_obs < cut_stratum_i[1]]
               cols_sel_i[[1]]     <- sample( set_ih, size_ih_i, replace = FALSE)
               if( n_stratum > 2 ) ## more than 2 cut_points, so fill in middle strata
               {
                    for(h in 2:(n_stratum-1))
                    {
                         set_ih              <- cols_obs[(cols_obs >= cut_stratum_i[h-1]) & 
                                                                   (cols_obs < cut_stratum_i[h])]
                         cols_sel_i[[h]]     <- sample( set_ih, size_ih_i, replace = FALSE )
                    }
               }## end sampling of middle strata (for each i) for n_stratum > 2
               set_ih                   <- cols_obs[cols_obs > cut_stratum_i[(n_stratum_i-1)]]
               size_ih                  <- sub_size[i] - size_ih_i*(n_stratum_i-1) ## what's left
               cols_sel_i[[n_stratum_i]]  <- sample( set_ih, size_ih_i, replace = FALSE)
               ## the "-1" accounts for the C++ indexing that starts at 0 for
               ## selecting columns of y, Omega_t and Omega_s[[1:L_s]]
               ## y_index[[i]] is the index of time_points (columns), not the response values
               y_index[[i]]             <- sort(unlist(cols_sel_i)) - 1
          } ## end loop i over n sub_sizes (tempered transition steps)
     } ## end drawing of stratified samples for each i in (1,...,n), split by n_stratum
     w_cols    <- apply(y_obs,2,function(x){var(x,na.rm=TRUE)})
     
     
     ## indicator to display progress bar
     progress            <- as.integer(progress)
     b_move              <- as.integer(b_move)

     ################################################################
     ## conduct posterior sampling and capture results
     ################################################################
     
     if( cluster == TRUE ) ## DP(GP) mix
     {
          if( any(is.na(y)) ) ## co-sample {bb_i} in the case data are missing.
          {
               y[is.na(y)] <- -9 
               ## this value will trigger an update step in c++ for missing y values in each row.
               res       <- gpdpbPost(y, ipr, Omega_t, Omega_s, gp_mod, jitter, b_move, gp_shape, gp_rate, 
                                      noise_shape, noise_rate, lower, 
                                      upper, w_star, w, n_slice_iter, y_index, n.iter, n.burn, n.thin, 
                                      n.tune, M_init, dp_shape, dp_rate, progress)
               res$optpartition$y[y == -9]   <- NA ## convert missing values back from -9 to NA
               ## need this object in cluster_plot() for plotting data vs. fitted functions
          }else{ ## cluster == TRUE and no missing values: 
               ## marginalize over {bb_i} for sampling the GP parameters.  May be quicker.
               res       <- gpdpPost(y, ipr, Omega_t, Omega_s, gp_mod, jitter, gp_shape, gp_rate, 
                                     noise_shape, noise_rate, lower, upper,
                                     w_star, w, n_slice_iter, y_index, n.iter, n.burn, n.thin, n.tune, M_init,
                                     dp_shape, dp_rate, progress)    
          } ## end conditional statement on whether any missing response values, y      
     }else{ ## no clustering
          if( any(is.na(y)) ) ## co-sample {bb_i} in the case data are missing.
          {
               y[is.na(y)] <- -9 
               if( fix == TRUE ) ## user inputs clustering
               {
                    s         <- s - 1 ## c++ index
                    
               }else{ ## user does not input a fixed clustering 
                    ## create a single cluster - M = 1
                    N         <- nrow(y) ## number of observation units
                    s         <- rep(1,N) - 1 ## c++ index
               } ## end conditional statement on whether user inputs a fixed clustering structure
               res       <- gpBFixPost(y, ipr, Omega_t, Omega_s, gp_mod, jitter, gp_shape, gp_rate, 
                                       noise_shape, noise_rate, lower, upper,
                                       w, n_slice_iter, y_index, n.iter, n.burn, n.thin, n.tune, 
                                       progress, s)
               res$optpartition$y[y == -9]   <- NA ## convert missing values back from -9 to NA
          }else{ ## no clustering and no missing values
               if( fix == TRUE )
               {
                    s         <- s - 1 ## c++ index
               }else{ ## no clustering
                    ## create a single cluster - M = 1
                    N         <- nrow(y) ## number of observation units
                    s         <- rep(1,N) - 1 ## c++ index
               } ## end conditional statement on whether user inputs a fixed clustering structure 
               res       <- gpFixPost(y, ipr, Omega_t, Omega_s, gp_mod, jitter, gp_shape, gp_rate, 
                                      noise_shape, noise_rate, lower, upper,
                                      w, n_slice_iter, y_index, n.iter, n.burn, n.thin, n.tune, 
                                      progress, s)
          } ## end conditional statement on whether there are any missing values.
          
     } ## end conditional statement on whether user wants clustering (e.g. a DP prior)
     cat(paste("Your chosen GP covariance terms includes term = ", c(gp_cov), "\n",sep=""))
     ## add in sn_order to construct test set covariance matrix used for prediction
     res$sn_order                  <- sn_order
     res$residuals                 <- colMeans(res$Residuals)
     
     ## add some column names
     ## Theta
     T                             <- ncol(y)
     N                             <- nrow(y)
     P                             <- sum(res$gp_indices$P_vec)
     if( cluster == TRUE ) ## each i in 1,..N has its own unique P x 1, Theta_i
     {
          colnames(res$Theta)           <- paste("theta[",rep(1:P,times=N),",",rep(1:N,each=P),"]",sep="")
     }else{ ## just P global parameters
          if( fix == TRUE ) ## fixed clustering scheme
          {
               colnames(res$Theta)           <- paste("theta[",rep(1:P,times=N),",",rep(1:N,each=P),"]",
                                                      sep="") ## P is fast-moving; N is slow
               M                             <- length(unique(s))
               B                             <- nrow(res$Theta)
               thetastar_vec                 <- vector("list",B)
               for(b in 1:B)
               {
                    thetastar_vec[[b]]      <- as.vector(res$Theta_star[[b]]) ## P is fast, M is slow
               }
               res$Thetastar                <- do.call("rbind",thetastar_vec)
               colnames(res$Thetastar)      <- paste("theta_star[",rep(1:P,times=M),",",rep(1:M,each=P),"]",
                                                     sep="")
          }else{ ## no clustering
               res$Theta                     <- res$Theta[,1:P] ## columns are redundant due 
                                                                ## to single cluster
               colnames(res$Theta)           <- paste("theta[",1:P,"]",sep="")
          }
          
     }
     
     ## bb
     colnames(res$bb)              <- paste("bb[",rep(1:N,times=T),",",rep(1:T,each=N),"]",sep="")
     ## Tau_e
     colnames(res$Tau_e)           <- "tau_e"
     ## assign column labels to each function, f_1, ... , f_K, where bb = f_1 + ... + f_K
     ## only draw additive f's if don't sample bb during the MCMC.
     if( all(!is.na(res$optpartition$y)) )
     {
          for( l in 1:L )  
          {
               colnames(res$f[[l]])     <- paste("f_",l,"[",rep(1:N,times=T),",",rep(1:T,each=N),"]",
                                                 sep="")
          } ## end loop over K additive function terms
          
     } ## finish assign column labels to additive f's.
     
     class(res)     	          <- c("gpdpgrow")
     return(res)
     
} ## end function gpdpgrow()

#####################################################
## .Call statements to C++ functions
#####################################################

#####################################################
## function to perform MCMC sampling for the GP parameter
## that marginalizes over the GP function, {bb_i}.
#####################################################
#' Run a Bayesian functional data model under a GP prior whose parameters employ a DP prior
#'
#' An internal function to \code{\link{gpdpgrow}}
#'
#' @export
#' @aliases gpdpPost
#' @param y An \emph{N x T} matrix of N observations of \emph{T x 1} functions
#' @param ipr An optional input vector of inclusion probabilities for each observation unit in the case
#'        the observed data were acquired through an informative sampling design, so that unbiased
#'        inference about the population requires adjustments to the observed sample.  Defaults to
#'        \code{ipr = rep(1,nrow(y))} indicating an iid sample.
#' @param Omega_t A \emph{T x T} matrix of squared Eucidean distances for \code{T} time points
#' @param Omega_s A \code{list} object of length \code{L_s}, where each contains the \emph{T x T}
#'        matrix of Euclidean distances associated to each seasonal covariance term.
#' @param gp_mod An \emph{L x 1} numeric vector denoting the selected covariance function for each
#'        of \code{L} terms.  \code{gp_mod = 1} is \code{"rq"}.  \code{gp_mod = 2} is \code{"se"}.
#'        \code{gp_mod = 3} is \code{"sn"}.
#' @param jitter Numeric value added to diagonals of GP covariance matrix to stabilize inversion.
#' @param gp_shape The shape parameter of the Gamma base distribution for the \code{kappa_star}
#'        locations used to sample the DP prior on the \code{P} GP covariance parameters, 
#'        \code{kappa}, for each experimental unit.
#' @param gp_rate The rate parameter of the Gamma base distribution for the \code{kappa_star}
#'        locations used to sample the DP prior on the \code{P} GP covariance parameters, 
#'        \code{kappa}, for each experimental unit.  
#' @param noise_shape The shape parameter of the Gamma base distribution on \code{tau_e}, the
#'        model noise precision parameter. Defaults to \code{noise_shape = 3}.
#' @param noise_rate The rate parameter of the Gamma base distribution on \code{tau_e}, the model
#'        noise precision parameter.  Defaults to \code{noise_rate = 1}.
#' @param lower Minimum in range of support for GP covariance parameters, \code{kappa}.
#' @param upper Maximum in range of support for GP covariance parameters, \code{kappa}. 
#' @param w_star Tuning parameter for number of locations to sample not linked to observations
#'        in the auxiliary Gibbs sampler for cluster assignments. 
#' @param w Tuning parameter for slice sampling interval width used for GP 
#'        covariance parameters, \code{kappa}. 
#' @param n_slice_iter Maximum number of steps to widen slice samplind width for
#'        GP covariance parameters, \code{kappa}.
#' @param y_index List object where each contains index of time points to use in \code{n}
#'        progressively coarser distribution for sampling \code{kappa} in tempered update steps.
#' @param n.iter The number of MCMC sampling iterations
#' @param n.burn The number of warm-up iterations to discard
#' @param n.thin The interval or step size of post-burn-in samples to return
#' @param n.tune The number of tuning iterations to update the slice sampler width, \code{w}.
#' @param M_init Starting value of number of clusters for sampling cluster assignments.
#' @param dp_shape The shape parameter for the \eqn{\Gamma} prior on the DP concentration parameter.  
#'     The rate parameter is set of \code{1}.
#' @param dp_rate The rate parameter for the \eqn{\Gamma} prior on the DP concentration parameter. 
#'        Default value is \code{1}.
#' @param progress An indicator in \code{{0,1}} denoting whether to display a progress bar during model execution.
#'        \code{progress = 1} displays a progress bar. Defaults to \code{progress = 1}.
#' @return res A list object containing MCMC runs for all model parameters.
#' @seealso \code{\link{gpdpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{gpdpgrow}}
gpdpPost = function (y, ipr, Omega_t, Omega_s, gp_mod, jitter, gp_shape, gp_rate, noise_shape, noise_rate, 
                        lower, upper, w_star, w, n_slice_iter, y_index, n.iter, n.burn, n.thin, 
                        n.tune, M_init, dp_shape, dp_rate, progress) 
{
                              y  <- as.matrix(y)
                              stopifnot(length(Omega_s) >= length(gp_mod[gp_mod == 3]))
                              stopifnot(ncol(y) == nrow(Omega_t))
                              res <- .Call("GPDPMIX", y, Omega_t, Omega_s, gp_mod, jitter, 
                                           gp_shape, gp_rate, noise_shape, noise_rate,
                                           lower, upper, w_star, w, 
                                           n_slice_iter, y_index, n.iter, n.burn, 
                                           n.thin, n.tune, M_init,
                                           dp_shape, dp_rate, progress, ipr, PACKAGE = "growfunctions")
} ## end function gpdpPost

#####################################################
## function to perform MCMC sampling for the GP parameter
## that marginalizes over the GP function, {bb_i}, under
## no clustering - just a simple GP model
#####################################################
#' Run a Bayesian functional data model under a GP prior whose parameters employ a DP prior
#'
#' An internal function to \code{\link{gpdpgrow}}
#'
#' @export
#' @aliases gpPost
#' @param y An \emph{N x T} matrix of N observations of \emph{T x 1} functions
#' @param ipr An optional input vector of inclusion probabilities for each observation unit in the case
#'        the observed data were acquired through an informative sampling design, so that unbiased
#'        inference about the population requires adjustments to the observed sample.  Defaults to
#'        \code{ipr = rep(1,nrow(y))} indicating an iid sample.
#' @param Omega_t A \emph{T x T} matrix of squared Eucidean distances for \code{T} time points
#' @param Omega_s A \code{list} object of length \code{L_s}, where each contains the \emph{T x T}
#'        matrix of Euclidean distances associated to each seasonal covariance term.
#' @param gp_mod An \emph{L x 1} numeric vector denoting the selected covariance function for each
#'        of \code{L} terms.  \code{gp_mod = 1} is \code{"rq"}.  \code{gp_mod = 2} is \code{"se"}.
#'        \code{gp_mod = 3} is \code{"sn"}.
#' @param jitter Numeric value added to diagonals of GP covariance matrix to stabilize inversion.
#' @param gp_shape The shape parameter of the Gamma base distribution for the \code{kappa_star}
#'        locations used to sample the DP prior on the \code{P} GP covariance parameters, 
#'        \code{kappa}, for each experimental unit.
#' @param gp_rate The rate parameter of the Gamma base distribution for the \code{kappa_star}
#'        locations used to sample the DP prior on the \code{P} GP covariance parameters, 
#'        \code{kappa}, for each experimental unit.  
#' @param noise_shape The shape parameter of the Gamma base distribution on \code{tau_e}, the
#'        model noise precision parameter. Defaults to \code{noise_shape = 3}.
#' @param noise_rate The rate parameter of the Gamma base distribution on \code{tau_e}, the model
#'        noise precision parameter.  Defaults to \code{noise_rate = 1}.
#' @param lower Minimum in range of support for GP covariance parameters, \code{kappa}.
#' @param upper Maximum in range of support for GP covariance parameters, \code{kappa}.  
#' @param w Tuning parameter for slice sampling interval width used for GP 
#'        covariance parameters, \code{kappa}. 
#' @param n_slice_iter Maximum number of steps to widen slice samplind width for
#'        GP covariance parameters, \code{kappa}.
#' @param y_index List object where each contains index of time points to use in \code{n}
#'        progressively coarser distribution for sampling \code{kappa} in tempered update steps.
#' @param n.iter The number of MCMC sampling iterations
#' @param n.burn The number of warm-up iterations to discard
#' @param n.thin The interval or step size of post-burn-in samples to return
#' @param n.tune The number of tuning iterations to update the slice sampler width, \code{w}.
#' @param progress An indicator in \code{{0,1}} denoting whether to display a progress bar during model execution.
#'        \code{progress = 1} displays a progress bar. Defaults to \code{progress = 1}.
#' @return res A list object containing MCMC runs for all model parameters.
#' @seealso \code{\link{gpdpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{gpdpgrow}}
gpPost = function (y, ipr, Omega_t, Omega_s, gp_mod, jitter, gp_shape, gp_rate, noise_shape, noise_rate,
                     lower, upper, w, n_slice_iter, y_index, n.iter, n.burn, n.thin, n.tune, progress) 
{
     y  <- as.matrix(y)
     stopifnot(length(Omega_s) >= length(gp_mod[gp_mod == 3]))
     stopifnot(ncol(y) == nrow(Omega_t))
     res <- .Call("GP", y, Omega_t, Omega_s, gp_mod, jitter, 
                  gp_shape, gp_rate, noise_shape, noise_rate, lower, upper, w, 
                  n_slice_iter, y_index, n.iter, n.burn, 
                  n.thin, n.tune, progress, ipr, PACKAGE = "growfunctions")
} ## end function gpPost


#####################################################
## function to perform MCMC sampling for the GP parameter
## that marginalizes over the GP function, {bb_i}, under
## fixed / input clustering - user sets cluster structure - not estimated
#####################################################
#' Run a Bayesian functional data model under a GP prior whose parameters employ a DP prior
#'
#' An internal function to \code{\link{gpdpgrow}}
#'
#' @export
#' @aliases gpFixPost
#' @param y An \emph{N x T} matrix of N observations of \emph{T x 1} functions
#' @param ipr An optional input vector of inclusion probabilities for each observation unit in the case
#'        the observed data were acquired through an informative sampling design, so that unbiased
#'        inference about the population requires adjustments to the observed sample.  Defaults to
#'        \code{ipr = rep(1,nrow(y))} indicating an iid sample.
#' @param Omega_t A \emph{T x T} matrix of squared Eucidean distances for \code{T} time points
#' @param Omega_s A \code{list} object of length \code{L_s}, where each contains the \emph{T x T}
#'        matrix of Euclidean distances associated to each seasonal covariance term.
#' @param gp_mod An \emph{L x 1} numeric vector denoting the selected covariance function for each
#'        of \code{L} terms.  \code{gp_mod = 1} is \code{"rq"}.  \code{gp_mod = 2} is \code{"se"}.
#'        \code{gp_mod = 3} is \code{"sn"}.
#' @param jitter Numeric value added to diagonals of GP covariance matrix to stabilize inversion.
#' @param gp_shape The shape parameter of the Gamma base distribution for the \code{kappa_star}
#'        locations used to sample the DP prior on the \code{P} GP covariance parameters, 
#'        \code{kappa}, for each experimental unit.
#' @param gp_rate The rate parameter of the Gamma base distribution for the \code{kappa_star}
#'        locations used to sample the DP prior on the \code{P} GP covariance parameters, 
#'        \code{kappa}, for each experimental unit.  
#' @param noise_shape The shape parameter of the Gamma base distribution on \code{tau_e}, the
#'        model noise precision parameter. Defaults to \code{noise_shape = 3}.
#' @param noise_rate The rate parameter of the Gamma base distribution on \code{tau_e}, the model
#'        noise precision parameter.  Defaults to \code{noise_rate = 1}.
#' @param lower Minimum in range of support for GP covariance parameters, \code{kappa}.
#' @param upper Maximum in range of support for GP covariance parameters, \code{kappa}.  
#' @param w Tuning parameter for slice sampling interval width used for GP 
#'        covariance parameters, \code{kappa}. 
#' @param n_slice_iter Maximum number of steps to widen slice samplind width for
#'        GP covariance parameters, \code{kappa}.
#' @param y_index List object where each contains index of time points to use in \code{n}
#'        progressively coarser distribution for sampling \code{kappa} in tempered update steps.
#' @param n.iter The number of MCMC sampling iterations
#' @param n.burn The number of warm-up iterations to discard
#' @param n.thin The interval or step size of post-burn-in samples to return
#' @param n.tune The number of tuning iterations to update the slice sampler width, \code{w}.
#' @param progress An indicator in \code{{0,1}} denoting whether to display a progress bar during model execution.
#'        \code{progress = 1} displays a progress bar. Defaults to \code{progress = 1}.
#' @param s An integer vector inputting cluster membership structure if select \code{fix == TRUE}.
#' @return res A list object containing MCMC runs for all model parameters.
#' @seealso \code{\link{gpdpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{gpdpgrow}}
gpFixPost = function (y, ipr, Omega_t, Omega_s, gp_mod, jitter, gp_shape, gp_rate, noise_shape, noise_rate, 
                      lower, upper, w, n_slice_iter, y_index, n.iter, n.burn, n.thin, n.tune, progress, s) 
{
     y  <- as.matrix(y)
     stopifnot(length(Omega_s) >= length(gp_mod[gp_mod == 3]))
     stopifnot(ncol(y) == nrow(Omega_t))
     res <- .Call("GPFIX", y, Omega_t, Omega_s, gp_mod, jitter, 
                  gp_shape, gp_rate, noise_shape, noise_rate, lower, upper, w, 
                  n_slice_iter, y_index, n.iter, n.burn, 
                  n.thin, n.tune, progress, s, ipr, PACKAGE = "growfunctions")
} ## end function gpFixPost


#####################################################
## function to perform MCMC sampling for the GP parameter
## that co-samples the GP function, {bb_i}.
## Needed in the case of intermittant missingness
## in data matrix, y
#####################################################

#' Run a Bayesian functional data model under a GP prior whose parameters employ a DP prior
#'
#' An internal function to \code{\link{gpdpgrow}}
#'
#' @export
#' @aliases gpdpbPost
#' @param y An \emph{N x T} matrix of N observations of \emph{T x 1} functions. \code{y} may\
#'        have intermittant missing values. They should be input with \code{NA}.
#' @param ipr An optional input vector of inclusion probabilities for each observation unit in the case
#'        the observed data were acquired through an informative sampling design, so that unbiased
#'        inference about the population requires adjustments to the observed sample.  Defaults to
#'        \code{ipr = rep(1,nrow(y))} indicating an iid sample.
#' @param Omega_t A \emph{T x T} matrix of squared Eucidean distances for \code{T} time points
#' @param Omega_s A \code{list} object of length \code{L_s}, where each contains the \emph{T x T}
#'        matrix of Euclidean distances associated to each seasonal covariance term.
#' @param gp_mod An \emph{L x 1} numeric vector denoting the selected covariance function for each
#'        of \code{L} terms.  \code{gp_mod = 1} is \code{"rq"}.  \code{gp_mod = 2} is \code{"se"}.
#'        \code{gp_mod = 3} is \code{"sn"}.
#' @param jitter Numeric value added to diagonals of GP covariance matrix to stabilize inversion.
#' @param b_move An indicator in \code{{0,1}} denoting whether to sample GP functions, \code{(bb_i)}
#'        in a \emph{T x 1} Gibbs step or through elliptical slice sampling.
#'        \code{b_move = 1} samples in a Gibbs step. Defaults to \code{b_move = 1}.
#' @param gp_shape The shape parameter of the Gamma base distribution for the \code{kappa_star}
#'        locations used to sample the DP prior on the \code{P} GP covariance parameters, 
#'        \code{kappa}, for each experimental unit.
#' @param gp_rate The rate parameter of the Gamma base distribution for the \code{kappa_star}
#'        locations used to sample the DP prior on the \code{P} GP covariance parameters, 
#'        \code{kappa}, for each experimental unit.  
#' @param noise_shape The shape parameter of the Gamma base distribution on \code{tau_e}, the
#'        model noise precision parameter. Defaults to \code{noise_shape = 3}.
#' @param noise_rate The rate parameter of the Gamma base distribution on \code{tau_e}, the model
#'        noise precision parameter.  Defaults to \code{noise_rate = 1}.
#' @param lower Minimum in range of support for GP covariance parameters, \code{kappa}.
#' @param upper Maximum in range of support for GP covariance parameters, \code{kappa}. 
#' @param w_star Tuning parameter for number of locations to sample not linked to observations
#'        in the auxiliary Gibbs sampler for cluster assignments. 
#' @param w Tuning parameter for slice sampling interval width used for GP 
#'        covariance parameters, \code{kappa}. 
#' @param n_slice_iter Maximum number of steps to widen slice samplind width for
#'        GP covariance parameters, \code{kappa}.
#' @param y_index List object where each contains index of time points to use in \code{n}
#'        progressively coarser distribution for sampling \code{kappa} in tempered update steps.
#' @param n.iter The number of MCMC sampling iterations
#' @param n.burn The number of warm-up iterations to discard
#' @param n.thin The interval or step size of post-burn-in samples to return
#' @param n.tune The number of tuning iterations to update the slice sampler width, \code{w}.
#' @param M_init Starting value of number of clusters for sampling cluster assignments.
#' @param dp_shape The shape parameter for the \eqn{\Gamma} prior on the DP concentration parameter.  
#'     The rate parameter is set of \code{1}.
#' @param dp_rate The rate parameter for the \eqn{\Gamma} prior on the DP concentration parameter. 
#'        Default value is \code{1}.
#' @param progress An indicator in \code{{0,1}} denoting whether to display a progress bar during model execution.
#'        \code{progress = 1} displays a progress bar. Defaults to \code{progress = 1}.
#' @return res A list object containing MCMC runs for all model parameters.
#' @seealso \code{\link{gpdpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{gpdpgrow}}
gpdpbPost = function (y, ipr, Omega_t, Omega_s, gp_mod, jitter, b_move, gp_shape, gp_rate, noise_shape,
                     noise_rate, lower, upper, w_star, w, n_slice_iter, y_index, n.iter, n.burn, 
                     n.thin, n.tune, M_init, dp_shape, dp_rate, progress) 
{
     y  <- as.matrix(y)
     stopifnot(length(Omega_s) >= length(gp_mod[gp_mod == 3]))
     stopifnot(ncol(y) == nrow(Omega_t))
     res <- .Call("GPDPMIXMIS", y, Omega_t, Omega_s, gp_mod, jitter, 
                  b_move, gp_shape, gp_rate, noise_shape, noise_rate, lower, upper, w_star, w, 
                  n_slice_iter, y_index, n.iter, n.burn, 
                  n.thin, n.tune, M_init,
                  dp_shape, dp_rate, progress, ipr, PACKAGE = "growfunctions")
} ## end function gpdpbPost

#####################################################
## function to perform MCMC sampling for the GP parameter
## that co-samples the GP function, {bb_i}, under
## fixed / input clustering - user sets cluster structure - not estimated
#####################################################
#' Run a Bayesian functional data model under a GP prior with a fixed clustering structure
#' that co-samples latent functions, {bb_i}.
#'
#' An internal function to \code{\link{gpdpgrow}}
#'
#' @export
#' @aliases gpBFixPost
#' @param y An \emph{N x T} matrix of N observations of \emph{T x 1} functions. \code{y} may\
#'        have intermittant missing values. They should be input with \code{NA}.
#' @param ipr An optional input vector of inclusion probabilities for each observation unit in the case
#'        the observed data were acquired through an informative sampling design, so that unbiased
#'        inference about the population requires adjustments to the observed sample.  Defaults to
#'        \code{ipr = rep(1,nrow(y))} indicating an iid sample.
#' @param Omega_t A \emph{T x T} matrix of squared Eucidean distances for \code{T} time points
#' @param Omega_s A \code{list} object of length \code{L_s}, where each contains the \emph{T x T}
#'        matrix of Euclidean distances associated to each seasonal covariance term.
#' @param gp_mod An \emph{L x 1} numeric vector denoting the selected covariance function for each
#'        of \code{L} terms.  \code{gp_mod = 1} is \code{"rq"}.  \code{gp_mod = 2} is \code{"se"}.
#'        \code{gp_mod = 3} is \code{"sn"}.
#' @param jitter Numeric value added to diagonals of GP covariance matrix to stabilize inversion.
#' @param gp_shape The shape parameter of the Gamma base distribution for the \code{kappa_star}
#'        locations used to sample the DP prior on the \code{P} GP covariance parameters, 
#'        \code{kappa}, for each experimental unit.
#' @param gp_rate The rate parameter of the Gamma base distribution for the \code{kappa_star}
#'        locations used to sample the DP prior on the \code{P} GP covariance parameters, 
#'        \code{kappa}, for each experimental unit.  
#' @param noise_shape The shape parameter of the Gamma base distribution on \code{tau_e}, the
#'        model noise precision parameter. Defaults to \code{noise_shape = 3}.
#' @param noise_rate The rate parameter of the Gamma base distribution on \code{tau_e}, the model
#'        noise precision parameter.  Defaults to \code{noise_rate = 1}.
#' @param lower Minimum in range of support for GP covariance parameters, \code{kappa}.
#' @param upper Maximum in range of support for GP covariance parameters, \code{kappa}.  
#' @param w Tuning parameter for slice sampling interval width used for GP 
#'        covariance parameters, \code{kappa}. 
#' @param n_slice_iter Maximum number of steps to widen slice samplind width for
#'        GP covariance parameters, \code{kappa}.
#' @param y_index List object where each contains index of time points to use in \code{n}
#'        progressively coarser distribution for sampling \code{kappa} in tempered update steps.
#' @param n.iter The number of MCMC sampling iterations
#' @param n.burn The number of warm-up iterations to discard
#' @param n.thin The interval or step size of post-burn-in samples to return
#' @param n.tune The number of tuning iterations to update the slice sampler width, \code{w}.
#' @param progress An indicator in \code{{0,1}} denoting whether to display a progress bar during model execution.
#'        \code{progress = 1} displays a progress bar. Defaults to \code{progress = 1}.
#' @param s An integer vector inputting cluster membership structure if select \code{fix == TRUE}.
#' @return res A list object containing MCMC runs for all model parameters.
#' @seealso \code{\link{gpdpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{gpdpgrow}}
gpBFixPost = function (y, ipr, Omega_t, Omega_s, gp_mod, jitter, gp_shape, gp_rate, noise_shape, 
                       noise_rate, lower, upper, w, n_slice_iter, y_index, n.iter, n.burn, 
                       n.thin, n.tune, progress, s) 
{
     y  <- as.matrix(y)
     stopifnot(length(Omega_s) >= length(gp_mod[gp_mod == 3]))
     stopifnot(ncol(y) == nrow(Omega_t))
     res <- .Call("GPBFIX", y, Omega_t, Omega_s, gp_mod, jitter, 
                  gp_shape, gp_rate, noise_shape, noise_rate, lower, upper, w, 
                  n_slice_iter, y_index, n.iter, n.burn, 
                  n.thin, n.tune, progress, s, ipr, PACKAGE = "growfunctions")
} ## end function gpBFixPost


##########################
## S3 multiple dispatch for predict_functions()
##########################
#' Use the model-estimated covariance parameters from gpdpgrow() or gmrdpgrow to predict the function at
#' future time points.  
#'
#' This is the generic predict_functions method. See the following functions
#' for the details about different data structures:
#'
#' \itemize{
#'   \item \code{\link{predict_functions.gpdpgrow}} for objects of class \code{gpdpgrow}
#'   \item \code{\link{predict_functions.gmrfdpgrow}} for objects of class \code{gmrfdpgrow}
#' }
#'
#' @keywords manip
#' @param object Object of class \code{gpdpgrow} or\code{gmrfdpgrow()}.
#' @param J Scalar denoting number of draws to take from posterior predictive for each unit.
#'          Defaults to \code{J = 500}.
#' @param ... further arguments passed to or from other methods.
#' @export
predict_functions           <- 	function(object, J, ...){
          if(is.null(attr(object, "class"))){
               print("object must be of class gpdpgrow or gmrfdpgrow")
          }
     else  UseMethod("predict_functions", object)
}

#####################################################
## function to predict GP at future time values
#####################################################
#' Use the model-estimated GP covariance parameters from gpdpgrow() to predict the GP function at
#' future time points.  Inputs the \code{gpdpgrow} object of estimated parameters.
#'
#' A companion function to \code{\link{gpdpgrow}}
#'
#' @export 
#' @param object Object of class \code{gpdpgrow} returned from model run of \code{gpdpgrow()}
#' @param J Scalar denoting number of draws to take from posterior predictive for each unit.
#'          Defaults to \code{J = 500}.
#' @param test_times A numeric vector holding test times at which to predict GP function values
#'             Will use the estimated covariance parameters from the training data to predict
#'             functions at the test_times for the \code{N} observation units.
#' @param time_points Inputs a vector of common time points at which the collections of functions were
#'        observed (with the possibility of intermittent missingness).  The length of \code{time_points}
#'        should be equal to the number of columns in the data matrix, \code{y}.  Defaults to 
#'        \code{time_points = 1:ncol(y)}.
#' @param sn_order An integer vector of length, \code{L_s}, equal to the number of seasonal terms.
#'                  Conveys the order of the seasonality for each term on the scale of T; for example,
#'                  if T is dimensioned in months, and one wishes to model quarterly seasonality, then
#'                  the applicable seasonality term would be of order \code{3}.         
#' @param ... further arguments passed to or from other methods.
#' @return out  A list object containing containing two matrices; the first is a K x (N*T)
#'                  matrix of predicted function values for each of K sampled iterations.  N is 
#'                  slow index and denotes the number of experimental units.  The second matrix is
#'                  an N x T average over the K sampled draws, composed in Rao-Blackwellized fashion.
#' @aliases predict_functions.gpdpgrow 
#' @export predict_functions gpdpgrow 
#' @seealso \code{\link{gpdpgrow}}
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
#' res_gp              = gpdpgrow(y = y_short
#'                               n.iter = 50, 
#'                               n.burn = 25, 
#'                               n.thin = 1, 
#'                               n.tune = 0) 
#'
#' ## Prediction Model Runs
#' T_test             <- 4
#' T_yshort           <- ncol(y_short)
#' pred_gp            <- predict_functions( object = res_gp, 
#'                        test_times = (T_yshort+1):(T_yshort+T_test) )
#'
#' ## plot estimated and predicted functions
#' plot_gp         <- predict_plot(object = pred_gp, 
#'                                units_label = cps$st, 
#'                                single_unit = FALSE, 
#'                                credible = TRUE)  
#' }
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as a companion function for \code{\link{gpdpgrow}} for prediction
predict_functions.gpdpgrow = function (object, J = 500, test_times, time_points = NULL,  
                                       sn_order = NULL,...) 
{
     ## generate trend and seasonal asymmetric Euclidean distance matrices
     ## for input to the c++ predict_bb() function
     y              	     <- object$optpartition$y
     y[is.na(y)]              <- -9 ## change back missing values to -9 for c++ prediction function
     object$optpartition$y    <- y
     T              	     <- ncol(y)
     T_test         	     <- length(test_times)
     if( is.null(time_points) )
     {
          time_points   <- 1:T
     }
     
     ## trend covariance test-train
     Omegat_tetr        <- t(sapply(1:T_test,function(l){
                         	sapply(1:T,function(j){
                              	d = (test_times[l]-time_points[j])^2
                              	})
                           }))
     
     ## trend covariance test-test
     Omegat_tete        <- t(sapply(1:T_test,function(l){
                              sapply(1:T_test,function(j){
                                   d = (test_times[l]-test_times[j])^2
                                   })
                            }))
     
     
     L_s            <- length(sn_order) ## will be 0 if is.null(sn_order) == TRUE
     if(L_s > 0)
     {
          ## test-train seasonal covariance matrix
          Omegas_tetr        <- array(0,dim=c(T_test,T,L_s))
          ## test-test seasonal covariance matrix
          Omegas_tete        <- array(0,dim=c(T_test,T_test,L_s))
          for( i in 1:L_s )
          {         ## test-train
                    Omegas_tetr[,,i]   <- t(sapply(1:T_test,function(l){
                                             sapply(1:T,function(j){
                                                  d = 2*( sin(pi*(test_times[l]-time_points[j])
                                                              /sn_order[i]) )^2 
                                             })
                                          }))
                    
                    ## test-test
                    Omegas_tete[,,i]   <- t(sapply(1:T_test,function(l){
                                             sapply(1:T_test,function(j){
                                                  d = 2*( sin(pi*(test_times[l]-test_times[j])
                                                            /sn_order[i]) )^2 
                                             })
                                         }))
          } ## end loop i over L_s seasonality terms
     }else{ ## no seasonality - compose dummy Omega_s for input to c++.  won't be used.
          ## te-tr
          Omegas_tetr              <- array(0,dim=c(T_test,T,1))
          Omegas_tetr[,,1]         <- Omegat_tetr
          ## te-te
          Omegas_tete              <- array(0,dim=c(T_test,T_test,1))
          Omegas_tete[,,1]         <- Omegat_tete
     } ## end conditional statement on seasonality to construct Euclidean distance matrices 

     out <- .Call("predict_bb", object, Omegas_tetr, Omegas_tete, Omegat_tetr, 
                    Omegat_tete, J, PACKAGE = "growfunctions")
     
     class(out)                    <- c("predict_functions.gpdpgrow")
     return(out)
} ## end function predict_functions()

     
####################################
## accessor methods
####################################
#' Produce samples of MCMC output
#'
#' provides posterior sampled values for every model parameter of a
#' \code{gpdpgrow} object
#'
#' @param object A \code{gpdpgrow} object
#' @param ... Ignored
#' @export 
#' @aliases samples.gpdpgrow  
samples      	<- 	function(object,...){
     if(is.null(attr(object, "class"))){
          print("object must be of class gpdpgrow or gmrfdpgrow")
     }
     else  UseMethod("samples", object)
}

#' @export 
samples.gpdpgrow <- function(object,...)
{
     res            <- list(Theta = object$Theta, bb = object$bb, Tau_e = object$Tau_e,
                            M = object$M, Conc = object$Conc)
     if( !any(is.na(object$optpartition$y)) ) ## in which case, sampled individual f's
     {
          res$f     <- object$f
     }
     class(res) 	<- "samples.gpdpgrow"
     return(res)
}
     

     