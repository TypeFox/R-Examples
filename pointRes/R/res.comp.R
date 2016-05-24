#' Calculate resilience components: resistance, recovery, resilience and relative resilience
#'
#' @description The function calculates resilience components on a \code{data.frame} of tree-ring series after Lloret et al. (2011), useful to analyze growth of individual trees prior, during and after extreme events / disturbances. The component 'resistance' is conceptually identical to 'abrupt growth changes' as described in Schweingruber et al. (1990). To identify negative event and pointer years, thresholds can be set as for the function \code{\link{pointer.rgc}}. 'Recovery' is the ability of tree growth to recover after disturbance, whereas 'resilience' reflects the ability of trees to reach pre-disturbance growth levels. Weighting of the resilience by the experienced growth reduction results in 'relative resilience'.
#'
#' @usage res.comp(data, nb.yrs = 4, res.thresh.neg = 40, series.thresh = 75)
#'
#' @param data a \code{data.frame} with tree-ring series as columns and years as rows (e.g., output of \code{read.rwl}, \code{bai.in} or \code{bai.out} of package dplR)
#' @param nb.yrs an \code{integer} specifying the number of years for pre- and post-disturbance periods to be considered in calculating resilience components. Defaults to 4.
#' @param res.thresh.neg a \code{numeric} specifying the threshold below which the resistance, expressed as a percentual change (i.e. relative growth reduction), is considered a negative event year for individual trees and years. Defaults to 40.
#' @param series.thresh a \code{numeric} specifying the minimum percentage of trees that should display a negative event year for that year to be considered as negative pointer year. Defaults to 75.
#' 
#' @details The function calculates the resilience components resistance, recovery, resilience and relative resilience as described in Lloret et al. (2011). A threshold on resistance can be set to identify negative event years for trees (cf. \code{\var{rgc.thresh.neg}} in function \code{\link{pointer.rgc}}), which are used to define negative pointer years for the site.
#'
#' If \code{\var{nb.yrs}}, \code{\var{res.thresh.neg}} and \code{\var{series.thresh}} are set to 4, 40 and 75 respectively, a negative pointer year will be defined when at least 75\% of the tree-ring series display an event year with resistance values indicating a growth decrease of at least 40\%, relative to the average growth in the 4 preceding years. The output provides the resilience components for all possible years, as well as for the selected pointer years separately.
#' 
#' Note that the resulting time series are truncated by \code{\var{nb.yrs}} at both ends inherent to the calculation methods. 
#'
#' @return 
#' #' The function returns a \code{list} containing the following components:
#' \item{resist}{a \code{matrix} with resistance values (i.e. relative growth changes) for individual tree-ring series}
#' \item{EYvalues}{a \code{matrix} indicating negative (-1) and non-event years (0) for individual tree-ring series}
#' \item{recov}{a \code{matrix} with recovery values for individual tree-ring series}
#' \item{resil}{a \code{matrix} with resilience values for individual tree-ring series}
#' \item{rel.resil}{a \code{matrix} with relative resilience values for individual tree-ring series}
#' \item{out}{a \code{data.frame} containing the following columns:}
#' \item{}{\code{year} - time stamp}
#' \item{}{\code{nb.series} - number of series considered}
#' \item{}{\code{perc.neg} - percentage of trees showing a negative event year}
#' \item{}{\code{nature} - number indicating whether the year is a negative (-1) or no pointer year (0)}
#' \item{}{\code{resist_mean} - mean resistance as percentual change over the available series}
#' \item{}{\code{resist_sd} - standard deviation of the resistance}
#' \item{}{\code{recov_mean} - mean recovery as percentual change over the available series}
#' \item{}{\code{recov_sd} - standard deviation of the recovery}
#' \item{}{\code{resil_mean} - mean resilience as percentual change over the available series}
#' \item{}{\code{resil_sd} - standard deviation of the resilience}
#' \item{}{\code{rel.resil_mean} - mean relative resilience calculated over the available series}
#' \item{}{\code{rel.resil_sd} - standard deviation of the relative resilience}
#' \item{out.select}{a \code{data.frame} containing a subset of rows from \code{out} that provide all statistics for years that were identified as negative pointer years}
#' \item{spec.param}{a \code{data.frame} specifying the arguments used in the calculation}
#' 
#' @author Marieke van der Maaten-Theunissen and Ernst van der Maaten.
#'
#' @references Lloret, F., Keeling, E.G. and Sala, A. (2011) Components of tree resilience: effects of successive low-growth episodes in old ponderosa pine forests. \emph{Oikos} 120: 1909-1920.
#' @references Schweingruber, F.H., Eckstein, D., Serre-Bachet, F. and \enc{Bräker}{Braker}, O.U. (1990) Identification, presentation and interpretation of event years and pointer years in dendrochronology. \emph{Dendrochronologia} 8: 9-38.
#'
#' @examples ## Calculate resilience components on tree-ring series
#' data(s033)
#' res <- res.comp(s033, nb.yrs = 4, res.thresh.neg = 40, series.thresh = 75)
#' res$out
#' res$out.select
#' 
#' @import stats
#' 
#' @export res.comp
#' 
res.comp <- function(data, nb.yrs = 4, res.thresh.neg = 40, series.thresh = 75)
{
  stopifnot(is.numeric(nb.yrs), length(nb.yrs) == 1, is.finite(nb.yrs))
  if(nb.yrs < 1) {
    stop("'nb.yrs' must be > 0")
  }
  stopifnot(is.numeric(res.thresh.neg), length(res.thresh.neg) == 1, 
            is.finite(res.thresh.neg))
  if(res.thresh.neg < 0) {
    stop("'res.thresh.neg' must be > 0")
  }
  if(res.thresh.neg > 100) {
    warning("res.thresh.neg' > 100 is unusual")
  }
  stopifnot(is.numeric(series.thresh), length(series.thresh) == 
              1, is.finite(series.thresh))
  if(series.thresh < 0 || series.thresh > 100) {
    stop("'series.thresh' must range from 0 to 100")
  }
  data2 <- as.matrix(data)
  if (!is.matrix(data2)) {
    stop("'data' must be coercible to a matrix")
  }
  if(ncol(data2) == 1) {
    stop("'data' must contain more than one series")
  }
  rnames <- rownames(data2)
  if (is.null(rnames)) {
    stop("'data' must have explicit row names")
  }
  yrs <- as.numeric(rnames)
  nyrs <- length(yrs)
  if (nyrs < 2 * nb.yrs + 1) {
    stop("'data' must be longer than the full calculation window (2 * nb.yrs + 1)")
  }
  
  start <- nb.yrs + 1
  
  avg.pre <- matrix(nrow = nrow(data2) - 2 * nb.yrs, ncol = ncol(data2))
  for(i in start:(nyrs-nb.yrs)) {
    avg.pre[i - nb.yrs,] <- colMeans(data2[(i - nb.yrs):(i - 1),])
  }
  rownames(avg.pre) <- yrs[start : (length(yrs) - nb.yrs)]
  resist <- data2[start : (nrow(data2) - nb.yrs), , drop = FALSE] / avg.pre[, , drop = FALSE]
  
  neg.thresh <- 1 - res.thresh.neg/100
  EYvalues <- ifelse(resist <= neg.thresh, -1, 0) 
  
  avg.post <- matrix(nrow = nrow(data2) - 2*nb.yrs, ncol = ncol(data2))
  for(i in start:(nyrs - nb.yrs)) {
    avg.post[i - nb.yrs,] <- colMeans(data2[(i + 1):(i + nb.yrs),])
  }
  rownames(avg.post) <- yrs[start : (length(yrs) - nb.yrs)]
  recov <- avg.post[, , drop = FALSE] / data2[start : (nrow(data2) - nb.yrs), , drop = FALSE]
  colnames(recov) <- colnames(data2)
  
  resil <- avg.post / avg.pre
  colnames(resil) <- colnames(data2)
  
  rel.resil <- (avg.post - data2[start:(nyrs - nb.yrs), ]) / avg.pre
  colnames(rel.resil) <- colnames(data2)
  
            year <- yrs[start : (length(yrs) - nb.yrs)]
       nb.series <- rowSums(!is.na(resist))
        perc.neg <- rowSums(resist <= neg.thresh, na.rm = TRUE)/nb.series * 100
             nat <- pmax(0, perc.neg - (series.thresh - 1e-07))
          nature <- ifelse(nat > 0, -1, 0)
     resist_mean <- (rowMeans(resist, na.rm = TRUE) - 1) * 100
       resist_sd <- apply(resist, 1, function(x) sd(x, na.rm = TRUE)) * 100
      recov_mean <- (rowMeans(recov, na.rm = TRUE) - 1) * 100
        recov_sd <- apply(recov, 1, function(x) sd(x, na.rm = TRUE)) * 100
      resil_mean <- (rowMeans(resil, na.rm = TRUE) - 1) * 100
        resil_sd <- apply(resil, 1, function(x) sd(x, na.rm = TRUE)) * 100
  rel.resil_mean <- rowMeans(rel.resil, na.rm = TRUE)
    rel.resil_sd <- apply(rel.resil, 1, function(x) sd(x, na.rm = TRUE))
  
  out <- data.frame(year, nb.series, perc.neg, nature, resist_mean, resist_sd,
                    recov_mean, recov_sd, resil_mean, resil_sd, rel.resil_mean,
                    rel.resil_sd, row.names = NULL)
  out.select <- subset(out, out[,"nature"] == -1)
  rownames(out.select) <- NULL
  
         out[,c(3, 5:12)] <- round(out[,c(3, 5:12)], 2)
  out.select[,c(3, 5:12)] <- round(out.select[,c(3, 5:12)], 2)

     resist <- round(resist, 2)
      recov <- round(recov, 2)
      resil <- round(resil, 2)
  rel.resil <- round(rel.resil, 2)
  
  spec.param <- data.frame(argument = c("nb.yrs","res.thresh.neg","series.thresh"), 
                           value = c(nb.yrs, res.thresh.neg, series.thresh))
  
  output <- list(resist = resist, EYvalues = EYvalues, recov = recov, resil = resil, 
                 rel.resil = rel.resil, out = out, out.select = out.select, spec.param = spec.param)
  class(output) <- "res.comp"
  return(output)
}