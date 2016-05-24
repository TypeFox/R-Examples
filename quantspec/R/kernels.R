################################################################################
#' Kernel function.
#'
#' Implementations of kernel functions
#'
#' Daniell kernel function \code{W0}:
#' \deqn{\frac{1}{2\pi} I\{|x| \leq \pi\}.}{1/(2pi) I{|x|<=pi}.}
#'
#' @name kernels
#' @aliases W0
#' @export
#'
#' @param x real-valued argument to the function; can be a vector
#'
#' @examples
#' plot(x=seq(-8,8,0.05), y=W0(seq(-8,8,0.05)), type="l")
################################################################################
W0 <- function(x){
  W0.simple <- function(x) {if (abs(x) <= pi) {1/(2*pi)} else {0}}
  return(Vectorize(W0.simple)(x))
}

################################################################################
#' @details
#' Epanechnikov kernel \code{W1} (i. e., variance minimizing kernel function of order 2):
#' \deqn{\frac{3}{4\pi} (1-\frac{x}{\pi})^2 I\{|x| \leq \pi\}.}{3/(4pi) (1-x/pi)^2 I{|x|<=pi}.}
#'
#'
#' @name kernels
#' @aliases W1
#' @export
#'
#' @examples
#' plot(x=seq(-8,8,0.05), y=W1(seq(-8,8,0.05)), type="l")
################################################################################
W1 <- function(x){
  W1.simple <- function(x) {if (abs(x) <= pi) {.75*(1-(x/pi)^2)/pi } else {0}}
  return(Vectorize(W1.simple)(x))
}

################################################################################
#' @details
#' Variance minimizing kernel function \code{W2} of order 4:
#' \deqn{\frac{15}{32\pi} (7(x/\pi)^4 -10(x/\pi)^2+3) I\{|x| \leq \pi\}.}{(15/(32 pi) (7 (x/pi)^4 - 10 (x/pi)^2 + 3) I{|x|<=pi}.}
#'
#' @name kernels
#' @aliases W2
#' @export
#'
#' @examples
#' plot(x=seq(-8,8,0.05), y=W2(seq(-8,8,0.05)), type="l")
################################################################################
W2 <- function(x){
  W2.simple <- function(x) {if (abs(x) <= pi) {1/pi * (15/32) * (7*(x/pi)^4 -10*(x/pi)^2+3)} else {0}}
  return(Vectorize(W2.simple)(x))
}

################################################################################
#' @details
#' Variance minimizing kernel function \code{W3} of order 6:
#' \deqn{\frac{35}{256\pi} (-99(x/\pi)^6 + 189(x/\pi)^4 - 105(x/\pi)^2+15) I\{|x| \leq \pi\}.}{(35/(256 pi) (-99(x/pi)^6 + 189(x/pi)^4 - 105(x/pi)^2+15) I{|x|<=pi}.}
#'
#' @name kernels
#' @aliases W3
#' @export
#'
#' @examples
#' plot(x=seq(-8,8,0.05), y=W3(seq(-8,8,0.05)), type="l")
################################################################################
W3 <- function(x){if (abs(x) <= pi) {1/pi * (35/256) * (-99*(x/pi)^6 + 189*(x/pi)^4 - 105*(x/pi)^2+15)} else {0}}
W3 <- function(x){
  W3.simple <- function(x) {if (abs(x) <= pi) {1/pi * (35/256) * (-99*(x/pi)^6 + 189*(x/pi)^4 - 105*(x/pi)^2+15)} else {0}}
  return(Vectorize(W3.simple)(x))
}

################################################################################
#' @details
#' Kernel yield by convolution of two Daniell kernels:
#' \deqn{\frac{1}{\pi+a} \Big(1-\frac{|x|-a}{\pi-a} I\{a \leq |x| \leq \pi\}\Big).}
#'
#' @name kernels
#' @aliases WDaniell
#' @export
#'
#' @param a real number between 0 and \eqn{\pi}{pi}
#'
#' @examples
#' plot(x=seq(-pi,pi,0.05), y=WDaniell(seq(-pi,pi,0.05),a=(pi/2)), type="l")
################################################################################
WDaniell <- function(x,a=(pi/2)){
  WDaniell.simple <- function(x) {
    if (abs(x) <= pi) {(1/(pi+a))*(1-(a <= abs(x) & abs(x) <= pi)*(abs(x)-a)/(pi-a))} else {0}
  }
  return(Vectorize(WDaniell.simple)(x))
}
################################################################################
#' @details
#' Parzen Window for lagEstimators
#'
#' @name kernels
#' @aliases WParzen
#' @export
#'
#' @param u real number
#'
#' @examples
#' plot(x=seq(-2,2,0.05),y=WParzen(seq(-2,2,0.05)),type = "l")
################################################################################
WParzen <- function(u){
  WParzen.simple <- function(u) {
    if (abs(u) <= 1){
      if (abs(u) <= .5){
        (1-6*u^2+6*abs(u)^3)
      }
      else{
        (2*(1-abs(u))^3)
      }
    }else {(0)}
  }
  return(Vectorize(WParzen.simple)(u))
}
