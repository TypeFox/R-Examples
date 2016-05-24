#' Convert Kernel Matrix to Distance Matrix
#'
#' \code{k2d} provides various methods for converting kernel matrix to distance matrix and vice versa. 
#'
#' There are various ways to convert kernel function values to distance between two points.
#' Normal-distance (when \code{method="norm"}) means a conversion
#' d_ND(x,y) = sqrt{ k(x,x) - 2k(x,y) + k(y,y)}.
#' 
#' Cauchy-Schwarz-type conversion (\code{method="CS"}) is more principled way:
#' d_CS(x,y) = arccos k^2(x,y)/k(x,x)k(y,y).
#' 
#' Other two simple ways are
#' d_exp(x,y) = exp(- k(x,y)/scale),
#' which is an exponential-type distance (\code{method="exp"}), and
#' d_n(x,y) = 1 - k(x,y)/sqrt{ k(x,x)k(y,y)},
#' which we call naive (\code{method="naive"}).
#'
#' For converting distance to kernel (\code{direction="d2k"}), it should be noted that we usually have distance between pairs of points only, and distances from "origin" are unknown.
#' Double-centering (\code{method="DC"}) is the most popular and simple way to convert distance to kernel. However, it does not make positive definite kernel in general,
#' and it sometimes require post-processing, e.g., cutting off negative eigenvalues (\code{pos=TRUE}).
#' Another simple way is exponential map (\code{method="exp"}), i.e., k(x,y) = exp( - d(x,y)/scale).
#'
#' @param Mat matrix, either kernel matrix or distance matrix to be converted.
#' @param direction a character string "k2d" or "d2k". The latter interpret the Mat as a distance matrix and convert it to a kernel matrix. Default "k2d".
#' @param method a character string to specify how the matrix is converted. "Default "norm".
#' @param scale a numeric parameter used to scale the matrix. Typically used in exp( - d(x,y)/scale).
#' @param pos logical. If \code{TRUE} when \code{direction="d2k"}, negative eigenvalues are round to zero to obtain positive semidefinite kernel matrix. Default \code{TRUE}.
#' @author Hideitsu Hino \email{hinohide@@cs.tsukuba.ac.jp}, Ken Takano, Yuki Yoshikawa, and Noboru Murata
#' @export
k2d <- function(Mat, direction="k2d",method="norm", scale=1, pos=TRUE){
  ## kernel to distance
  if(direction=="k2d"){
    num <- nrow(Mat)
    index <- 1:num
    comb.index <- combn(index,2)
    myD <- Mat-Mat

    if(method=="norm"){
      apply(comb.index,2,FUN=function(x){
            valxx <- Mat[x[1],x[1]]
            valxy <- Mat[x[1],x[2]]
            valyy <- Mat[x[2],x[2]]
            myD[x[1],x[2]]<<-myD[x[2],x[1]]<<- sqrt((valxx-2*valxy+valyy)^2)
          })
      diag(myD) <- 0

    }else if(method=="CS"){
      apply(comb.index,2,FUN=function(x){
            valxx <- Mat[x[1],x[1]]
            valxy <- Mat[x[1],x[2]]
            valyy <- Mat[x[2],x[2]]
            myD[x[1],x[2]]<<-myD[x[2],x[1]]<<- acos(valxy^2/(valxx*valyy))
})
      diag(myD) <- 0

    }else if(method=="exp"){
      apply(comb.index,2,FUN=function(x){
        valxx <- Mat[x[1],x[1]];valxy <- Mat[x[1],x[2]];valyy <- Mat[x[2],x[2]]
        ##myD[x[1],x[2]]<<-myD[x[2],x[1]]<<- exp(-valxy/(sqrt(valxx*valyy)*scale))}
        myD[x[1],x[2]]<<-myD[x[2],x[1]]<<- exp(-valxy/scale)}
      )
      diag(myD) <- 0

    }else if(method=="naive"){
      apply(comb.index,2,FUN=function(x){
        valxx <- Mat[x[1],x[1]];valxy <- Mat[x[1],x[2]];valyy <- Mat[x[2],x[2]]
            ## d_n(x,y) = 1 - k(x,y)/sqrt{ k(x,x)k(y,y)},
            myD[x[1],x[2]]<<-myD[x[2],x[1]]<<- 1- valxy/sqrt(valxx*valyy)}
      )
      diag(myD) <- 0

    }else{
      ## should return error with message
    }
    myD[which(is.na(myD))] <- max(myD,na.rm=TRUE)
    return(myD)
  }else{  ## distance to kernel
    if(method=="DC"){
      num <- nrow(Mat)
      index <- 1:num
      comb.index <- combn(index,2)
      myK <- Mat-Mat
      apply(comb.index,2,FUN=function(x){
            myK[x[1],x[2]]<<-myK[x[2],x[1]]<<- (-Mat[x[1],x[2]]+sum(Mat[,x[2]])/num+sum(Mat[,x[1]])/num - sum(Mat)/(num^2))/2
      })
      Mat <- myK; rm(myK)
    }else if(method=="exp"){
      apply(comb.index,2,FUN=function(x){
        valxx <- Mat[x[1],x[1]];valxy <- Mat[x[1],x[2]];valyy <- Mat[x[2],x[2]]
        myD[x[1],x[2]]<<-myD[x[2],x[1]]<<- exp(-valxy/scale)}
      )
      Mat <- myK; rm(myK)
    }else{
      ## should return error with message
    }
    
    ## keep the kernel matrix positive semidefinite
    if(pos){
      ##eig <- try(eigen(myK))
      eig <- eigen(Mat)
      ##if(class(eig)=="try-error"){browser()}
      eig$values[eig$values<0] <- 10^(-8)
      myK <- eig$vectors %*% diag(eig$values,length(eig$values)) %*% t(eig$vectors)
      myK <- .5*(myK+t(myK))
      return(myK)
    }
    return(myK)
  }
}
