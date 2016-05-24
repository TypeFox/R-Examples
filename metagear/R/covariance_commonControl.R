#' Generates a sampling variance-covariance matrix for modeling dependencies
#' among effect sizes due to sharing a common control.
#'
#' Generates K by K sampling variance-covariance (VCV) matrix that models the
#' dependencies that arise due to using the same control group study parameters 
#' when estimating multiple effect sizes.  This VCV matrix can then be used in
#' meta-analysis.  Currently only supports VCV calculation for log response
#' ratios (see Lajeunesse 2011).
#'
#' @param aDataFrame A data frame containing columns with all study parameters
#'    used to estimate effect sizes (e.g., means, SD, N's for treatment and 
#'    control groups).  Must also contain a column that codes which effect sizes
#'    share a common control.  See example below.
#' @param control_ID Label of the column that codes groups of effect sizes that 
#'  share the mean, SD, and N of a control group.  
#' @param X_t Column label for the means of (t)reatment group used to estimate the
#'    effect size.
#' @param SD_t Column label for the standard deviations (SD) of the treatment
#'    group used to estimate the effect size.
#' @param N_t Column label for the sample size (N) of the treatment group
#'    used to estimate the effect size.  
#' @param X_c Column label for the means of (c)ontrol group used to estimate the
#'    effect size.
#' @param SD_c Column label for the standard deviations (SD) of the control
#'    group used to estimate the effect size.
#' @param N_c Column label for the sample size (N) of the control group
#'    used to estimate the effect size. 
#' @param metric Option to designate which effect size metric for which the 
#'    common control VCV matrix is to be estimated.  Default is "RR" for log 
#'    response ratio.
#'
#' @return A K by K sampling variance-covariance matrix and a data frame aligned
#'    with the block diagonal design of the sampling matrix. 
#'
#' @references Lajeunesse, M.J. 2011. On the meta-analysis of response ratios for
#'     studies with correlated and multi-group designs. Ecology 92: 2049-2055. 
#'
#' @import Matrix
#' @export covariance_commonControl

covariance_commonControl <- function (aDataFrame, 
                                      control_ID, 
                                      X_t, 
                                      SD_t, 
                                      N_t,
                                      X_c, 
                                      SD_c, 
                                      N_c,
                                      metric = "RR") {
  
  ## generate list of control groups in dataframe
  controlList <- split(aDataFrame, as.factor(aDataFrame[, control_ID]))
  listV <- list(); dataAlignedWithV <- data.frame();
  
  for(i in 1:length(controlList)) { 
    ## stack dataframes in V order
    dataAlignedWithV <- rbind(dataAlignedWithV, controlList[[i]])
        
    if(metric == "RR") {    
      ## common control covariance and variance of response ratio based Lajeunesse 2011
      covar <- (controlList[[i]][, SD_c] ^ 2) / (controlList[[i]][, N_c] * (controlList[[i]][, X_c] ^ 2)) 
      var <- covar + (controlList[[i]][, SD_t] ^ 2) / (controlList[[i]][, N_t] * (controlList[[i]][, X_t] ^ 2)) 
    } 
        
    ## calculate V for ith element in the controlList
    V <- matrix(covar, nrow = length(var), ncol = length(var))
    diag(V) <- var
    
    ## collect V's with a list
    listV <- unlist(list(listV, list(V)), recursive = FALSE)
  } 
  
  ## convert list of V's into single matrix
  V <- as.matrix(bdiag(listV))
  
  ## return V matrix paired with aligned dataset
  return(list(V, dataAlignedWithV))
}
