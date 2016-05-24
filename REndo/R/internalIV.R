#'@title Constructs Internal Instrumental Variables From Data
#'@aliases internalIV
#'@description The function can be used to construct additional instruments to be supplied to \code{\link{hmlewbel}} as additional instruments
#' in the "EIV" argument.
#' @param    y   the vector or matrix containing the dependent variable.
#' @param    X   the data frame or matrix containing the exogenous regressors of the model.
#' @param    P   the endogenous variables of the model as columns of a matrix or dataframe.
#' @param    G   the functional form of G. It can take four values, \code{x2}, \code{x3},\code{lnx} or \code{1/x}. 
#' The last two forms are conditional on the values of the exogenous variables: greater than 0 or different from 0 respectively.
#' @param    IIV   the internal instrumental variable to be constructed. It can take six values, \code{"g","gp","gy","yp","p2"or"y2"}. See the "Details" section 
#' of \code{\link{hmlewbel}} for a description of the internal instruments.
#' @param    data  optional data frame or list containing the variables in the model.
#' @return    Returns a vector/matrix constructed from the data whcih can be used as instrumental variable either in \code{\link{hmlewbel}} or in any other function/algorithm making use
#' of instruments.
#' @references Lewbel, A. (1997). "Lewbel, A. (1997). 'Constructing Instruments for Regressions with Measurement Error when No Additional Data Are Available,
#' with An Application to Patents and R&D'. Econometrica, 65(5), 1201-1213."
#' @keywords endogeneity
#' @keywords instruments
#' @keywords lewbel
#' @seealso \code{\link{hmlewbel}}
#' @examples 
#' # load data
#' data(dataHMLewbel)
#'y <- dataHMLewbel$y
#'X <- cbind(dataHMLewbel$X1,dataHMLewbel$X2)
#'colnames(X) <- c("X1","X2")
#'P <- dataHMLewbel$P
#'# build an instrument gp = (G - mean(G))(P - mean(P))  using the internalIV() function 
#'# with G = "x3" meaning G(X) = X^3
#'eiv <- internalIV(y,X,P, G ="x3", IIV = "gp")
#'@export
# the endogenous variable should be on the last position in the formula
internalIV <- function(y,X,P, G = c("x2","x3","lnx","1/x"),IIV=c("g","gp","gy","yp","p2","y2"), data=NULL){

# function computing the de-meanned variable
s <- function(t){
    s1 <- t - mean(t)
    return(s1)
}   

if (G == "x3") {
  IIV1 <- s(X^3)    # q1i = s(X^3)
  IIV2 <- s(X^3) * s(P)  # q2i
  IIV3 <- s(X^3) * s(y)  # q3i
} else {

  if (G == "lnx"){
  if (any(X <= 0)) stop(" G cannot be computed: exogenous variables contain values less or equal to 0")
      else{
  IIV1 <- s(log10(X))    # q1i = s(X^3)
  IIV2 <- s(log10(X)) * s(P)  # q2i
  IIV3 <- s(log10(X)) * s(y)  # q3i
  }
} else {
  if (G == "1/x") {
  if (any(X==0)) stop("G cannot be computed: exogenous variables contain 0s")
  else
  {
  IIV1 <- s(1/X)    # q1i = s(X^3)
  IIV2 <- s(1/X) * s(P)  # q2i
  IIV3 <- s(1/X) * s(y)  # q3i
}
}
# here G = ()^2
else {
  if (G == "x2"){
IIV1 <- s(X^2)    # q1i = s(X^2)
IIV2 <- s(X^2) * s(P)  # q2i
IIV3 <- s(X^2) * s(y)  # q3i
}
}
}}
IIV4 <- s(y) * s(P)   # q41, 
IIV5 <- s(P)^2  # q5i
IIV6 <- s(y)^2  # q6i  
    
AIV <- NULL  
 
if (IIV=="g") {AIV <- IIV1} else {
    if (IIV == "gp") {AIV <- IIV2} else {
      if (IIV == "gy") {AIV <- IIV3} else {
        if (IIV == "yp") {AIV <- IIV4} else {
          if (IIV == "p2") {AIV <- IIV5} else {
            if (IIV == "y2") {AIV <- IIV6}
          }
        }
      }
    }
  }
  
return(AIV)

}

