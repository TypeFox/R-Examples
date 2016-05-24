TriPlot <- function(mhOut, sumStat = FALSE) {
  #Plots samples from a 3-simplex projected into 
  #2-dimensions
  
  #Check inputs
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
    stop("mhOut$type not recognized. Use 'dirichlet', 'multinom' or user'.")
  }
  
  if (!("dat" %in% names(mhOut))) {
    stop("mhOut must contain dat, although it value can be set to NULL")
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
  #Plot samples on triangle
  par(pty='s', mar=c(0, 0, 0, 0))
  s <- sin(2*pi/6)
  plot(mhOut$S[, 2] + mhOut$S[, 3]/2, mhOut$S[, 3]*s, type = 'p', pch = 20,
       xlim = c(0, 1), ylim = c(0, s), pty = 's', axes = FALSE, xlab = '',
       ylab = '', main = "", cex = .2)
  segments(c(0, 1, 0.5), c(0, 0, 1)*s, c(1, 0.5, 0), c(0, 1, 0)*s)
  Sbar <- apply(mhOut$S, 2, mean)
  
  #Calculate and plot points corresponding to theortical mean & sample mean
  if (sumStat == TRUE){
    points(mhOut$center[2] + mhOut$center[3]/2, mhOut$center[3]*s, col = 'red', 
           pch = 16, cex = 3)
    points(Sbar[2] + Sbar[3]/2, Sbar[3]*s, col = 'blue', pch = 'x', cex = 3)
    
    #For dirichlet, calculate and plot theoretical mode
    if (mhOut$type == 'dirichlet') {
      mode <- (mhOut$a - 1)/(sum(mhOut$a) - mhOut$p)
      points(mode[2] + mode[3]/2, mode[3]*s, col = 'green', pch = 16, cex = 3)
      legend("topright", fill = c("red", "blue", "green"), 
             legend = c("Theoretical Mean", "Sample Mean", "Theoretical Mode"))
    } else {
      legend("topright", legend = c("Theoretical Mean", "Sample Mean"), 
             fill = c("red", "blue"))
    }
  }
}