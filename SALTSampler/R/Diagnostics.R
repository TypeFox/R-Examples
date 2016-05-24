Diagnostics <- function(mhOut) {
  #Diagnostic graphs from RunMh output. The acceptance rate is outputted for 
  #each dimension along with a trace plot. For type dirichlet, qqplots 
  #of the theoretical versus empirical marginal distributions are also
  #provided for each dimension.
  
  #Check input
  if (!(is.list(mhOut))){
    stop("mhOut is not a list.")
  }
  
  if (("Y" %in% names(mhOut)) & !(is.matrix(mhOut$Y))){
    stop("mhOut$Y is not a matrix")
  }
  
  if (!("S" %in% names(mhOut))) {
    stop("mhOut must contain S")
  }
  if (!is.matrix(mhOut$S)){
    stop("mhOut$S is not a matrix")
  }
  
  if (!("runTime" %in% names(mhOut))) {
    stop("mhOut must contain runTime")
  }
  if (!(is(mhOut$runTime) == "proc_time")) {
    stop("mhOut$runTime is not a 'proc_time'")
  }
  
  if (!("moveCount" %in% names(mhOut))) {
    stop("mhOut must contain moveCount")
  }
  if (!is.numeric(mhOut$moveCount)) {
    stop("mhOut$moveCount is not numeric")
  }
  if (!is.vector(mhOut$moveCount)) {
    stop("mhOut$moveCount is not a vector")
  }
  
  if (!("p" %in% names(mhOut))) {
    stop("mhOut must contain p")
  }
  if (!is.numeric(mhOut$p)) {
    stop("mhOut$p is not numeric")
  }
  if (length(mhOut$p) != 1) {
    stop("mhOut$p is not of length 1")
  }
  if (mhOut$p%%1 != 0){
    stop("mhOut$p is not an integer")
  }
  
  if (!("center" %in% names(mhOut))) {
    stop("mhOut must contain center")
  }
  if (!is.numeric(mhOut$center)) {
    stop("mhOut$center is not numeric")
  }
  if (!is.vector(mhOut$center)) {
    stop("mhOut$center is not a vector")
  }
  
  if (!("B" %in% names(mhOut))) {
    stop("mhOut must contain B")
  }
  if (!is.vector(mhOut$B)){
    stop("mhOut$B is not a vector")
  }
  if (mhOut$B%%1 != 0){
    stop("mhOut$B is not an integer")
  }
  if (length(mhOut$B) != 1){
    stop("mhOut$B is not of length 1")
  }
  
  if (!("concentration" %in% names(mhOut))) {
    stop("mhOut must contain concentration")
  }
  if (!is.vector(mhOut$concentration)){
    stop("mhOut$concentration is not a vector")
  }
  if (length(mhOut$concentration) != 1){
    stop("mhOut$concentration is not of length 1")
  }
  
  if (("h" %in% names(mhOut)) & (!is.numeric(mhOut$h))){
    stop("mhOut$h is not a numeric vector")
  }
  if (("h" %in% names(mhOut)) & (!is.vector(mhOut$h))){
    stop("mhOut$h is not a numeric vector")
  }
  
  if (!("type" %in% names(mhOut))) {
    stop("mhOut must contain type")
  }
  if (!(mhOut$type == 'dirichlet' || mhOut$type == 'user' || mhOut$type == 'multinom')) {
    stop("mhOut$type not recognized. Use 'dirichlet', 'multinom' or 'user'.")
  }
  
  if (!("dat" %in% names(mhOut))) {
    stop("mhOut must contain dat, although its value can be set to NULL")
  }
  if (!is.null(mhOut$dat) & !is.matrix(mhOut$dat)) {
    stop("mhOut$dat is not a matrix or NULL")
  }
  
  if (!("a" %in% names(mhOut))) {
    stop("mhOut must contain a")
  }
  if (!is.numeric(mhOut$a)) {
    stop("mhOut$a is not numeric")
  }
  if (!is.vector(mhOut$a)) {
    stop("mhOut$a is not a vector")
  }
  
  if (("h" %in% names(mhOut)) & (length(mhOut$center) != length(mhOut$h))){
    stop("Length of mhOut$center does not equal length of mhOut$h")
  }
  
  #Traceplots on true scale and logit scale
  par(oma = c(2,2,2,2))
  acceptRate <- mhOut$moveCount/mhOut$B 
  matplot(mhOut$S, type='l', lty=1, main = "Trace Plots on True Scale", ylab = "S", 
          xlab = "Iteration")
  matplot(mhOut$Y, type='l', lty=1, main = "Trace Plots on Logit Scale", ylab = "Y", 
          xlab = "Iteration")
  
  #For type = 'dirichlet' qq plots for each dimension
  if (mhOut$type == 'dirichlet') {
    for (i in 1:mhOut$p) {
      qq <- qqmath(mhOut$S[, i], distribution=function(x) qbeta(x, mhOut$a[i],
                                                                sum(mhOut$a[-i])), panel=function(...) {
                                                                  panel.qqmath(...); panel.abline(0, 1)
                                                                }
                   , main = sprintf("QQPlot of Marginal Distribution, Dimension: %i", i), 
                   xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles")
      print(qq)
      panel.text(525, 800, sprintf("Acceptance Rate: %s", round(acceptRate, 3)[i]), cex = 1.1) 
    }
  }
  
  #Return acceptance rates
  return(list(acceptRate = acceptRate))
}
