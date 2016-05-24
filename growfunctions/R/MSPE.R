#' Compute normalized mean squared prediction error based on accuracy to impute missing data values
#'
#' Uses as input the output object from the gpdpgrow() and gmrfdpgrow() functions.
#'
#' @param object A \code{gpdpgrow} or \code{gmrfdpgrow} object. 
#' @param y_true An \code{N x T} numeric matrix of test set values.
#' @param pos An \code{N x T} matrix with all entries either \code{0} or {1}, where a \code{1}
#'        indexes a missing entry or test point in \code{y_true}.
#' @return A list object containing various MSPE fit statistics that measure the accuracy of 
#'         of predicting the values in \code{y_true} indexed by \code{pos}.
#'     \item{SSE}{Sum of squared errors based on \emph{full} \code{N x T}, \code{y_true - y_hat}.}
#'     \item{MSE}{Mean squared error computed from \code{SSE}.}
#'     \item{RMSE}{Square root of \code{MSE}.}
#'     \item{SSPE}{Sum of squared prediction error based on missing values.}
#'     \item{MSPE}{Mean squared prediction error based on missing values.}
#'     \item{nMSPE}{Mean squared prediction error based on missing values that is 
#'             normalized by the variance of the test observations to produce a value in 
#'             \code{[0,1]}.}
#'     \item{RMSPE}{Square root of \code{MSPE}.}
#' @seealso \code{\link{gpdpgrow}}, \code{\link{gmrfdpgrow}}
#' @examples 
#' {
#' library(growfunctions)
#' 
#' ## load the monthly employment count data for a collection of 
#' ## U.S. states from the Current 
#' ## Population Survey (cps)
#' data(cps)
#' ## subselect the columns of N x T, y, associated with 
#' ## the years 2009 - 2013
#' ## to examine the state level employment 
#' ## levels during the "great recession"
#' y_short   <- cps$y[,(cps$yr_label %in% c(2009:2013))]
#' 
#' ## dimensions
#' T         <- ncol(y_short)  ## time points per unit
#' N         <- nrow(y_short)  ## number of units
#' 
#' ## Demonstrate estimation of intermittent missing data 
#' ## from posterior predictive distribution by randomly
#' ## selecting 10 percent of entries in y_short and 
#' ## setting them to NA.
#' 
#' ## randomly assign missing positions in y.
#' ## assume every unit has equal number of 
#' ## missing positions
#' ## randomly select number of missing 
#' ## observations for each unit
#' m_factor  <- .1
#' M         <- floor(m_factor*N*T)
#' m_vec     <- rep(floor(M/N),N)
#' if( sum(m_vec) < M )
#' {
#'     m_left              <- M - sum(m_vec)
#'     pos_i               <- sample(1:N, m_left, 
#'                                 replace = FALSE)
#'      m_vec[pos_i]        <- m_vec[pos_i] + 1
#' } # end conditional statement on whether all 
#'   # missing cells allocated
#'   # randomly select missing 
#'   # positions for each unit
#' pos       <- matrix(0,N,T)
#' for( i in 1:N )
#' {
#'     sel_ij              <- sample(3:(T-3), m_vec[i], 
#'                            replace = FALSE) 
#'                            ## avoid edge effects
#'     pos[i,sel_ij]       <- 1
#' }
#'
#' ## configure N x T matrix, y_obs, with 10 percent missing 
#' ## observations (filled with NA)
#' ## to use for sampling.  Entries in y_short 
#' ## that are set to missing (NA) are
#' ## determined by entries of "1" in the 
#' ## N x T matrix, pos.
#' y_obs               <- y_short
#' y_obs[pos == 1]     <- NA       
#'
#' ## Conduct dependent GP model estimation under 
#' ## missing observations in y_obs.
#' ## We illustrate the ability to have multiple 
#' ## function or covariance terms
#' ## by seting gp_cov = c("se","sn"), which indicates 
#' ## the first term is a
#' ## squared exponential ("se") trend covariance term 
#' ## and the "sn" is a seasonality
#' ## term.  The setting, sn_order = 4, indicates the 
#' ## length scale of the seasonality
#' ## term is 4 month.  The season term is actually 
#' ## "quasi" seasonal, in that the
#' ## seasonal covariance kernel is multiplied by a 
#' ## squared exponential, which allows
#' ## the pattern of seasonality to evolve over time.
#' res_gp_2               <- gpdpgrow(y = y_obs, 
#'                                   gp_cov=c("se","sn"), 
#'                                   sn_order = 4, 
#'                                   n.iter = 10, 
#'                                   n.burn = 4, 
#'                                   n.thin = 1, 
#'                                   n.tune = 0) 
#'
#' ## 2 plots of estimated functions: 1. faceted by cluster and fit;
#' ## 2.  data for experimental units.
#' ## for a group of randomly-selected functions
#' fit_plots_gp_2        <- cluster_plot( object = res_gp_2,  
#'                                       units_name = "state", 
#'                                       units_label = cps$st, 
#'                                       single_unit = TRUE, 
#'                                       credible = TRUE )
#'                                       
#' ## Compute out-of-sample fit statistic, normalized mean-square 
#' ## prediction error (MSPE)
#' ## The normalized MSPE will take the predicted values 
#' ## for the entries in y_short purposefully
#' ## set to NA and will subtract them from the known true 
#' ## values in y_short (and square them).
#' ## This MSE on predicted (test) data is then 
#' ## divided by the variance of the test observations
#' ## to output something akin to a percent error.
#' ( nmspe_gp  <- MSPE(res_gp_2, y_short, pos)$nMSPE )    
#' }
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases MSPE
#' @export
MSPE      <-  function(object, y_true, pos)
{
     ## check inputs
     if( !(class(object) %in% c("gpdpgrow","gmrfdpgrow")) ) 
          stop("\n'object' must be a return object from either gpdpgrow() or gmrfdpgrow().\n")
     ## set dimensions 
     N            = nrow(y_true)
     T            = ncol(y_true)
     M            = sum(sum(pos)) ## number of missing observations 
     
     ## extract parameter and data objects 
     y_hat                   = object$optpartition$y_bar ## N x T
        
     ## compute MSE 
     ## extract true and predicted values for missing units 
     y_hat_test                 = y_hat[pos == 1];
     y_true_test                = y_true[pos == 1];
     y_true_diff                = y_true_test - mean(y_true_test);
     ## compute MSE and MSPE 
        SSE       = sum(sum((y_true-y_hat)^2));
        MSE       = 1/(N*T) * SSE;
        RMSE      = sqrt(MSE);
        SSPE      = sum((y_true_test - y_hat_test)^2);
        MSPE      = 1/(N*M - 1) * SSPE;
        totVar    = 1/(N*M - 1) * sum((y_true_diff)^2);
        nMSPE     = MSPE / totVar;
        RMSPE     = sqrt(MSPE);
        
        return(list(SSE = SSE, MSE = MSE, RMSE = RMSE, SSPE = SSPE, MSPE = MSPE,
                    nMSPE = nMSPE, RMSPE = RMSPE))
         
} ## end function MSPE 