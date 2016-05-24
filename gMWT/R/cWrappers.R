# Version: 28-06-2013, DF

# Changes:
# Changed general layout and added comments, 28-06-2013, DF
# Added the function: getP.Csub, 30-11-2012, DF

# Calculate the probabilistic indices for pairs and triples, using the C code
  getP.Cnaive <- function(x, y , z = NULL){
    if(is.null(z))
    {
      result <- .C("getPR", as.double(x), as.double(y), as.integer(length(x)), as.integer(length(y)),result=as.double(1))$result
    } else {
      result <- .C("getPTripR", as.double(x), as.double(y), as.double(z), as.integer(length(x)), as.integer(length(y)), as.integer(length(z)),result=as.double(1))$result
    }

    result
  } # end of function getP.Cnaive
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function permutes two vectors - basically the same as sample(c(x,y)), I could substitute it...
  permObs2.C <- function(x, y){
    nx <- length(x)
    ny <- length(y)
    res <- .C("permObs2R", as.double(x), as.double(y), as.integer(nx), as.integer(ny),result=double(nx+ny))$result
    xOut <- res[1:nx]
    yOut <- res[(nx+1):(nx+ny)]
    return(list(x=xOut,y=yOut))
  } # End of function permObs2.C
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function calculates the UIT test statistic for a given set of vectors x,y,z
  uit.C <- function(x, y, z){
    .C("uitR", as.double(x), as.double(y), as.double(z), as.integer(length(x)), as.integer(length(y)), as.integer(length(z)),result=double(5))$result[1]
  } # End of function uit.C
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function calcualtes the RHO fr a given set of vectors x,y,z
  getRho.C <- function(x, y, z){
    .C("getRhoR", as.integer(length(x)), as.integer(length(y)), as.integer(length(z)),result=double(6))$result[1]
  } # end of function getRho.C
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function calculates a vector of length 'nper' of permutation test statistics in the triple case using the submatric approach
  getP.Csub <- function(x,y,z,nper){
    .Call( "subMat", c(x, y, z), length(x), length(y),length(z), nper)$result
  } # end of function getP.Csub
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

