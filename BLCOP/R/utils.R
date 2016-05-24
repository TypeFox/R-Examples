###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 11/08/2008 09:16:51
# newPMatrix
# Author: Francisco
###############################################################################
# DESCRIPTION: Creates a new pick matrix with appropriate column names
###############################################################################

newPMatrix <- function
(
    assetNames,       # set of assets referred to by the matrix
    numViews,         # number of views
    defaultValue = 0  # default value to use
)
{
    stopifnot(length(assetNames) > 0 && numViews > 0)
    matrix(defaultValue, ncol = length(assetNames), nrow = numViews, dimnames = list(NULL, assetNames))
}

.padMatrix <- function(x, targetRows, fillVal = 0)
{
  if(nrow(x) < targetRows)
  {

    fillRows <- matrix(fillVal, ncol = ncol(x), nrow = targetRows - nrow(x) )
    return(rbind(x, fillRows))
  }
  warning("x has enough rows")
  x
}

.padVector <- function(x, targetLength, fillVal = 0)
{
  if(length(x) < targetLength)
  {
    return(c(x, rep(fillVal, targetLength - length(x))))
  }
  return(x)  
}

.blockDiag <- function(A,B)
{
    stopifnot(class(A) == "matrix" && class(B) == "matrix")
    x <- ncol(A) + ncol(B)
    y <- nrow(A) + nrow(B)
    z <- matrix(0, ncol = x, nrow = y)
    z[1:nrow(A), 1:ncol(A)] <- A
    yOffset <- nrow(A)
    xOffset <- ncol(A)
    z[1:ncol(B)+xOffset, 1:nrow(B)+yOffset] <- B
    z
}

.assertClass <- function(object, classNames)
{
    if(! any(classNames %in% class(object)) )
    stop(paste("None of the classes:", classNames, ",were inherited by object"))
}

.removeZeroColumns <- function(mat) {
    
    isZeroColumn <- apply(mat == 0, 2, all)
    mat[,!isZeroColumn, drop = FALSE]
}

.correlationMatrix <- function(upperTriangle, dim)
{

    sigma <- matrix(0, nrow = dim, ncol = dim)
    sigma[upper.tri(sigma)] <- upperTriangle
    diag(sigma) <- 1   
    sigma <- t(sigma)
    sigma[upper.tri(sigma)] <- upperTriangle
    sigma
}

.varcovMatrix <- function(stdDeviations, correlations, dim)
{
    x <- .correlationMatrix(correlations, dim)
    x <- x * stdDeviations
    x <- t(t(x) * stdDeviations)
    x
}

.symmetricMatrix <- function(upperTriangle, dim)
{
    result <- matrix(NA, nrow = dim, ncol = dim)
    result[upper.tri(result, diag = TRUE)] <- upperTriangle
    result <- t(result)
    result[upper.tri(result, diag = TRUE)] <- upperTriangle
    result
}