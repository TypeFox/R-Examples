l5 <- structure(list(expr = "Fluo ~ c + (d - c)/((1 + exp(b * (log(Cycles) - log(e))))^f)", 
                     fct = function (x, parm) 
                     {
                       b <- parm[1]
                       c <- parm[2]
                       d <- parm[3]
                       e <- parm[4]
                       f <- parm[5]
                       c + (d - c)/((1 + exp(b * (log(x) - log(e))))^f)
                     }, ssFct = function (x, y) 
                     {
                       d <- max(y) + 0.001
                       c <- min(y) - 0.001
                       x2 <- x[y > 0]
                       y2 <- y[y > 0]
                       logitTrans <- log((d - y2)/(y2 - c))
                       lmFit <- lm(logitTrans ~ log(x2))
                       coefVec <- coef(lmFit)
                       b <- coefVec[2]
                       e <- exp(-coefVec[1]/b)
                       f <- 1
                       ssVal <- as.numeric(c(b, c, d, e, f))
                       names(ssVal) <- l5$parnames
                       return(ssVal)
                     }, d1 = function (x, parm) 
                     {
                       b <- parm[1]
                       c <- parm[2]
                       d <- parm[3]
                       e <- parm[4]
                       f <- parm[5]
                       b * (c - d) * e^-b * f * x^(-1 + b) * (1 + e^-b * x^b)^(-1 - 
                                                                                 f)
                     }, d2 = function (x, parm) 
                     {
                       b <- parm[1]
                       c <- parm[2]
                       d <- parm[3]
                       e <- parm[4]
                       f <- parm[5]
                       -b * (c - d) * e^(-2 * b) * f * x^(-2 + b) * (1 + e^-b * 
                                                                       x^b)^(-2 - f) * (-(-1 + b) * e^b + (1 + b * f) * 
                                                                                          x^b)
                     }, inv = function (y, parm) 
                     {
                       b <- parm[1]
                       c <- parm[2]
                       d <- parm[3]
                       e <- parm[4]
                       f <- parm[5]
                       e * (1/(-1 + ((c - d)/(c - y))^(1/f)))^(-1/b)
                     }, expr.grad = expression(c + (d - c)/((1 + exp(b * (log(Cycles) - 
                                                                            log(e))))^f)), inv.grad = expression(e * (1/(-1 + ((c - 
                                                                                                                                  d)/(c - Fluo))^(1/f)))^(-1/b)), parnames = c("b", "c", 
                                                                                                                                                                               "d", "e", "f"), name = "l5", type = "five-parameter log-logistic"), .Names = c("expr", 
                                                                                                                                                                                                                                                              "fct", "ssFct", "d1", "d2", "inv", "expr.grad", "inv.grad", "parnames", 
                                                                                                                                                                                                                                                              "name", "type"))


# Example of an artificial chamber dPCR experiment using the test data set from
# qpcR. The function limit_cq is used to calculate the Cy0 value and converts 
# all values between a defined range to 1 and the remaining to 0.


#' Limit Cy0 values
#' 
#' Calculates the Cq values of a qPCR experiment
#' within a defined range of cycles. The function can be used to extract Cq
#' values of a chamber based qPCR for conversion into a dPCR experiment. All Cq
#' values are obtained by Second Derivative Maximum or by Cy0 method (Guescini
#' et al. (2008)).
#' 
#' @details 
#' The \code{Cq_range} for this function an be defined be the user. The default
#' is to take all amplification curves into consideration. However, under
#' certain circumstances it is recommended to define a range. For example if
#' amplifications are positive in early cycle numbers (less than 10).
#' 
#' Approximated second derivative is influenced both by how often interpolation
#' takes place in each data interval and by the smoothing method used. The user
#' is encouraged to seek optimal parameters for his data himself. See
#' \code{\link[chipPCR]{inder}} for details.
#' 
#' The calculation of the Cy0 value (equivalent of Cq) is based on a
#' five-parameter function. From experience this functions leads to good
#' fitting and avoids overfitting of critical data sets. Regardless, the user
#' is recommended to test for the optimal fitting function himself (see
#' \code{\link[qpcR]{mselect}} for details).
#' 
#' @param data a dataframe containing the qPCR data.
#' @param cyc the column containing the cycle data. Defaults to first column.
#' @param fluo the column(s) (runs) to be analyzed. If \code{NULL}, all runs will be
#' considered (equivalent of \code{(1L:ncol(data))[-cyc]}).
#' @param Cq_range is a user defined range of cycles to be used for the
#' determination of the Cq values.
#' @param model is the model to be used for the analysis for all runs. Defaults
#' to 'l5' (see \code{\link[qpcR]{pcrfit}}).
#' @param SDM if \code{TRUE}, Cq is approximated by the second derivative
#' method.  If \code{FALSE}, Cy0 method is used instead.
#' @param pb if \code{TRUE}, progress bar is shown.
#' @return A data frame with two columns and number of rows equal to the number
#' of runs analyzed. The column \code{Cy0} contains calculated Cy0 values. The
#' column \code{in.range} contains adequate logical constant if given Cy0 value
#' is in user-defined \code{Cq_range}.
#' @author Michal Burdukiewicz, Stefan Roediger.
#' @seealso SDM method: \code{\link[chipPCR]{inder}},
#' \code{\link[chipPCR]{summary.der}}.
#' 
#' Cy0 method: \code{\link[qpcR]{mselect}}, \code{\link[qpcR]{efficiency}}.
#' @references Guescini M, Sisti D, Rocchi MB, Stocchi L & Stocchi V (2008)
#' \emph{A new real-time PCR method to overcome significant quantitative
#' inaccuracy due to slight amplification inhibition}. BMC Bioinformatics, 9:
#' 326.
#' 
#' Ruijter JM, Pfaffl MW, Zhao S, et al. (2013) \emph{Evaluation of qPCR curve
#' analysis methods for reliable biomarker discovery: bias, resolution,
#' precision, and implications}. Methods, San Diego Calif 59:32--46.
#' @keywords Cy0 qPCR dPCR
#' @examples
#' 
#' library(qpcR)
#' test <- cbind(reps[1L:45, ], reps2[1L:45, 2L:ncol(reps2)], reps3[1L:45, 
#' 	      2L:ncol(reps3)])
#' 
#' # results.dPCR contains a column with the Cy0 values and a column with 
#' # converted values.
#' Cq.range <- c(20, 30)
#' ranged <- limit_cq(data = test, cyc = 1, fluo = NULL,
#'                            Cq_range = Cq.range, model = l5)
#' # Same as above, but without Cq.range
#' no_range <- limit_cq(data = test, cyc = 1, fluo = NULL, model = l5)
#' 
#' # Same as above, but only three columns
#' no_range234 <- limit_cq(data = test, cyc = 1, fluo = c(2:4), model = l5)
#' 
#' @export limit_cq
limit_cq <- function(data, cyc = 1, fluo = NULL,
                     Cq_range = c(1, max(data[cyc])), model = l5, SDM = TRUE, pb = FALSE) {
  if (Cq_range[1] > Cq_range[2]) {
    warning("First value of Cq_range is greater than second. Cq_range reversed.")
    Cq_range <- rev(Cq_range)
  }
  
  if (is.null(fluo))
    fluo <- (1L:ncol(data))[-cyc]
  
  if(pb) 
    pb <- txtProgressBar(min = 1, max = length(fluo), initial = 0, style = 3)
  
  if (SDM) {
    Cy0 <- vapply(fluo, function(fluo_col) {
      if(pb) setTxtProgressBar(pb, fluo_col)
      summary(inder(x = data[, cyc], y = data[, fluo_col]), print = FALSE)["SDM"]
    }, 0)
  } else {
    Cy0 <- vapply(fluo, function(fluo_col) {
      if(pb) setTxtProgressBar(pb, fluo_col)
      
      fit <- pcrfit(data = data, cyc = cyc, fluo = fluo_col, model = model)
      
      if(all(is.nan(fit$MODEL$d2(eff(fit)[["eff.x"]], coef(fit))))) {
        stop("The derivative cannot be computed using the chosen method. Consider runnning with SDM = TRUE.")
      } else {
        efficiency(fit, type = "Cy0", plot = FALSE)[["Cy0"]]
      }
    }, 0)
  }
  
  Cy0.res <- vapply(Cy0, function(Cy0_i)
    Cq_range[1] <= Cy0_i & Cy0_i <= Cq_range[2], TRUE)
  
  data.frame(Cy0 = Cy0, in.range = Cy0.res)
}


