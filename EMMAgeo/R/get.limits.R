#' Create lower and upper mode position limits to define robust end-members.
#' 
#' This function identifies the lower and upper limits within which robust 
#' end-members have clustered mode positions. It uses a kernel density estimate
#' of the mode positions of all input end-member loadings, clips it at a 
#' user-defined minimum density and returns the resulting rising and falling 
#' shoulders of the kde peaks as limits.
#' 
#' Note that the threshold above which a mode cluster is identified is an 
#' arbitrary, user-defined value and probably needs to be adjusted iteratively
#' to get reasonable results. The default value may or may not be adequate! 
#' 
#' 
#' @param loadings Numeric matrix with m loadings (rows) and n classes 
#' (columns).
#' @param bw Numeric scalar, bandwidth of the kernel, which is moved over the 
#' data set. If omitted, the default value of 1 % of the number of classes is
#' used.
#' @param threshold Numeric scalar, threshold quantile which is used to 
#' identify mode clusters. Only kde densities above this values are kept and 
#' used to derieve mode cluster limits.
#' @return Numeric matrix with lower and upper mode limits.
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{EMMA}}, \code{\link{model.em}}
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.
#' @keywords EMMA
#' @examples
#' 
#' ## load example data set
#' data(X, envir = environment())
#' 
#' ## define parameters
#' l <- c(0, 0.1)
#' q <- rbind(c(2, 3),
#'            c(3, 4))
#' 
#' ## model all possible end-members
#' em.pot <- model.em(X = X, 
#'                    q = q, 
#'                    l = l)
#'                    
#' ## infer mode cluster limits
#' limits <- get.limits(loadings = em.pot$loadings)
#' 
#' @export get.limits
get.limits <- function(
  loadings,
  bw,
  threshold = 0.7
) {
  
  ## create mode vector
  loadings_mode <- numeric(nrow(loadings))
  
  ## fill mode vector
  for(i in 1:length(loadings_mode)) {
    loadings_mode[i] <- seq(from = 1, 
                            to = ncol(loadings))[
                              loadings[i,] == max(loadings[i,], 
                                                  na.rm = TRUE)]
  }
  
  ## check/set bw
  if(missing(bw) == TRUE) {
    bw <- (max(loadings_mode) - min(loadings_mode)) / 100
  }
  
  
  ## create kde of modes
  kde <- density(x = loadings_mode,
                 bw = bw)
  
  ## keep kde parts above threshold value and convert to limits
  kde.ok <- kde$y >= quantile(x = kde$y, probs = threshold)
  
  kde.limits.1 <- diff(x = kde.ok) == 1
  kde.limits.2 <- diff(x = kde.ok) == -1
  
  ## create limits matrix
  limits <- cbind(kde$x[kde.limits.1], 
                  kde$x[kde.limits.2])
  
  ## sort limits row-wise
  limits <- t(apply(X = limits, 
                    MARGIN = 1, 
                    FUN = sort))
  
  ## print threshold value
  print(paste("Threshold is", quantile(x = kde$y, probs = threshold)))
  
  ## return result
  return(limits)
}