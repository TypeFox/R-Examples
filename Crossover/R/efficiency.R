
getCounts <- function(design, long=FALSE, omit.balanced=TRUE) {
  s <- summary(as.factor(design))
  if (length(unique(s))==1) return("")
  if (long) return(paste(paste(names(s), " (", s,")", sep=""), collapse=", "))
  return(paste(s, collapse=", "))
}

#' Create the design matrix, variance-covariance matrix, the variance of each
#' pairwise comparison and the efficicency of each pairwise comparison for a
#' cross-over design
#' 
#' Function to read in a cross-over design and create the design matrix X, 
#' the variance of each pairwise comparison and the efficicency of each pairwise comparison.
#' 
#' See the vignette of this package for further details.
#' 
#' @param design Cross-over design.
#' @param model Model - one of the following: 1) "Standard additive model",
#' 2) "Second-order carry-over effects", 3) "Full set of interactions",
#' 4) "Self-adjacency model", 5) "Placebo model", 6) "No carry-over into self
#' model", 7) "Treatment decay model", 8) "Proportionality model", 9) "No carry-over effects". 
#' @param model.param List of additional model specific parameters. In the
#' moment these are \code{ppp}, the proportionality parameter for the
#' proportionality model, and \code{placebos}, the number of placebo treatments
#' in the placebo model.
#' @param v Number of treatments
#' @return A list with the following elements:
#' \itemize{
#' \item xmat Design matrix for the given model (including subject and period effects)
#' \item var.trt.pair.adj Matrix of treament difference variances
#' \item eff.trt.pair.adj Matrix of treament difference efficiencies
#' }
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @references Jones, B., & Kenward, M. G. (2003). Design and analysis of
#' cross-over trials (Vol. 98). Chapman & Hall.
#' @keywords misc
#' @examples
#' 
#' design.efficiency(getDesign("fletcher1"))
#' design.efficiency(getDesign("fletcher1"), model=7)
#' design.efficiency(getDesign("switchback4t"), model=7)
#' 
#' @export design.efficiency
design.efficiency <- function(design, model=1, model.param=list(), v=length(levels(as.factor(design)))) {
  p <- dim(design)[1]
  s <- dim(design)[2]
  
  #if(missing(C)) 
  {
    Csub <- contrMat(n=rep(1, v), type="Tukey")
    class(Csub) <- "matrix"
    C <- appendZeroColumns(Csub, model, v)
  }
  # TODO DO we have to add model.param to estimable?
  if (!estimable(design, v, model, C)) {
    warning(paste("Not all treatment contrasts are estimable for model", model, "in this design."))
    m <- matrix(NA, v, v)
    return(list(xmat=cbind(rcdMatrix(rcd(design, v, model), v=v, model=model),  
                           getZ(s=dim(design)[2],p=dim(design)[1])), 
                var.trt.pair.adj=m, 
                eff.trt.pair.adj=m))
  }
  m <- matrix(0, v, v)
  # Actual variances
  variances <- do.call(getValues, c(list(design, model, C=C, v=v), model.param))
  m[lower.tri(m)] <- variances
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  # Ideal design
  im <- matrix(0, v, v)
  vn <- sapply(1:v, function(x) {sum(design==x)})
  for (i in 1:v) {
    for (j in 1:v) {
      if (i!=j) {
        im[i, j] <- 1/vn[i]+1/vn[j]
      }
    }
  }
  # Efficiency:
  em <- im/m
  diag(em) <- 0
  return(list(xmat=cbind(rcdMatrix(rcd(design, v, model), v=v, model=model),  
                         getZ(s=dim(design)[2],p=dim(design)[1])), 
              var.trt.pair.adj=m, 
              eff.trt.pair.adj=em))
}