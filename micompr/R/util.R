# Copyright (c) 2016 Nuno Fachada
# Distributed under the MIT License (http://opensource.org/licenses/MIT)

#' Associate colors to \emph{p}-values
#'
#' Associate colors to \emph{p}-values according to their value.
#'
#' @param pvals Vector of \emph{p}-values to which associate colors.
#' @param col Colors Vector of colors to associate with the \emph{p}-values
#' given in \code{pvals} according to the limits specified in \code{pvlims}.
#' @param pvlims Vector of \emph{p}-value upper limits, first value should be 1.
#'
#' @return A vector of colors associated with \emph{p}-values given in
#' \code{pvals}.
#'
#' @keywords internal
#'
#' @examples
#' micompr:::pvalcol(c(0.06, 0.9, 0.0001, 0.3, 0.2, 0.02),
#'                   c("darkgreen", "yellow", "red"))
#' # [1] "darkgreen" "darkgreen" "red"       "darkgreen" "darkgreen" "yellow"
#'
pvalcol <- function(pvals, col, pvlims = c(0.05, 0.01)) {

  idxs <- sapply(pvals, function(p) match(F, p < c(pvlims, -Inf)))
  colors <- col[idxs]
  colors
}

#' Default colors for plots in \code{micomp} package
#'
#' Default colors for plots in \code{micomp} package.
#'
#' @return Default colors for plots in \code{micomp} package.
#'
#' @keywords internal
#'
#' @examples
#' micompr:::plotcols()
#' # [1] "blue"   "red"    "green"  "gold"   "violet" "cyan"
#'
plotcols <- function() {
  c("blue", "red", "green", "gold", "violet", "cyan")
}

#' Center and scale vector
#'
#' Center and scale input vector using the specified method.
#'
#' @param v Vector to center and scale.
#' @param type Type of scaling: "center", "auto", "range", "iqrange", "vast",
#' "pareto", "level" or "none".
#'
#' @return Center and scaled vector using the specified method.
#'
#' @references Berg, R., Hoefsloot, H., Westerhuis, J., Smilde, A., and Werf, M.
#' (2006). Centering, scaling, and transformations: improving the biological
#' information content of metabolomics data. \emph{BMC Genomics} \bold{7}, 142.
#' DOI: 10.1186/1471-2164-7-142
#'
#' @export
#'
#' @examples
#'
#' v <- c(-100, 3, 4, 500, 10, 25, -8, -33, 321, 0, 2)
#'
#' centerscale(v, "center")
#' # [1] -165.81818  -62.81818  -61.81818  434.18182  -55.81818  -40.81818
#' # [7]  -73.81818  -98.81818  255.18182  -65.81818  -63.81818
#'
#' centerscale(v, "auto")
#' # [1] -0.9308937 -0.3526577 -0.3470437  2.4374717 -0.3133601 -0.2291509
#' # [7] -0.4144110 -0.5547596  1.4325760 -0.3694995 -0.3582716
#'
#' centerscale(v, "range")
#' # [1] -0.2763636 -0.1046970 -0.1030303  0.7236364 -0.0930303 -0.0680303
#' # [7] -0.1230303 -0.1646970  0.4253030 -0.1096970 -0.1063636
#'
#' centerscale(v, "iqrange")
#' # [1] -6.085071 -2.305254 -2.268557 15.933278 -2.048374 -1.497915 -2.708924
#' # [8] -3.626355  9.364470 -2.415346 -2.341952
#'
#' centerscale(v, "vast")
#' # [1] -0.34396474 -0.13030682 -0.12823247  0.90064453 -0.11578638 -0.08467115
#' # [7] -0.15312466 -0.20498338  0.52933609 -0.13652987 -0.13238117
#'
#' centerscale(v, "pareto")
#' # [1] -12.424134  -4.706731  -4.631804  32.531614  -4.182247  -3.058353
#' # [7]  -5.530919  -7.404075  19.119816  -4.931509  -4.781657
#'
#' centerscale(v, "level")
#' # [1] -2.5193370 -0.9544199 -0.9392265  6.5966851 -0.8480663 -0.6201657
#' # [7] -1.1215470 -1.5013812  3.8770718 -1.0000000 -0.9696133
#'
#' centerscale(v, "none")
#' # [1] -100    3    4  500   10   25   -8  -33  321    0    2
#'
centerscale <- function(v, type) {
  switch(type,
         center = v - mean(v),
         auto = (v - mean(v)) / stats::sd(v),
         range = (v - mean(v)) / (max(v) - min(v)),
         # Keep results equal to MATLAB by specifying type = 5
         iqrange = (v - mean(v)) / stats::IQR(v, type = 5),
         vast = (v - mean(v)) * mean(v) / stats::var(v),
         pareto = (v - mean(v)) / sqrt(stats::sd(v)),
         level = (v - mean(v)) / mean(v),
         none = v)
}

#' Concatenate multiple outputs with multiple observations
#'
#' Concatenate multiple outputs with multiple observations.
#'
#' @param outputlist List of outputs. Each output is a \emph{n} x \emph{m}
#' matrix, where \emph{n} is the number of observations and \emph{m} is the
#' number of variables (i.e. output length).
#' @param centscal Centering and scaling method: "center", "auto", "range",
#' "iqrange", "vast", "pareto", "level" or "none". This task is delegated to the
#' \code{\link{centerscale}} function.
#'
#' @return An \emph{n} x \emph{p} matrix, representing the \emph{n} observations
#' of the concatenated output, each observation of length \emph{p}, which is the
#' sum of individual output lengths.
#'
#' @export
#'
#' @examples
#'
#' # Collect 20 observations of 3 outputs with different scales and lengths
#'
#' # Output 1, length 100
#' out1 <- matrix(rnorm(2000, mean = 0, sd = 1), nrow = 20)
#'
#' # Output 2, length 200
#' out2 <- matrix(rnorm(4000, mean = 100, sd = 200), nrow = 20)
#'
#' # Output 1, length 50
#' out3 <- matrix(rnorm(1000, mean = -1000, sd = 10), nrow = 20)
#'
#' # Concatenate and range scale outputs, resulting matrix dimensions will be
#' # 20 x 350
#' outconcat <- concat_outputs(list(out1, out2, out3), "range")
concat_outputs <- function(outputlist, centscal = "none") {

  # Check if it's a list
  if (!is.list(outputlist)) {
    stop("'outputlist' argument is not a list")
  }

  # Check if it's not an empty list
  if (length(outputlist) == 0) {
    stop("'outputlist' is an empty list")
  }

  # Determine number of observations and output lengths
  odims <- sapply(outputlist, dim)
  nobs = unique(odims[1, ])
  outlen = odims[2, ]

  # Check that the number of observations is constant
  if (length(nobs) > 1) {
    stop("Number of observations (rows) is not the same for each output matrix")
  }

  # Concatenate outputs
  outconcat <- matrix(nrow = nobs, ncol = sum(outlen))
  for (i in 1:nobs) {
    outconcat[i, ] <- unlist(sapply(outputlist, function(x)
      centerscale(x[i, ], type = centscal)))
  }

  # Return concatenated output
  outconcat

}
