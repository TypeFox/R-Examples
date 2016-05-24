
#Wolfgang Konen
#Cologne University of Applied Science
#RbfInter.R
#
# RBF interpolation model
#
 
# euclidean distance of x to all xp
# 
distLine <- function(x,xp) {
  z = outer(rep(1,nrow(xp)),x) - xp
  z = sqrt(rowSums(z*z))
}


#----------------------------------------------------------------------------------
#' Fit cubic RBF interpolation to training data for d>1.
#' 
#' The model for a point \eqn{z=(z_1,...,z_d)} is fitted using n sample points \eqn{x_1, ..., x_n} 
#' \cr
#'    \deqn{ s(z) = \lambda_1*\Phi(||z-x_1||)+... +\lambda_n*\Phi(||z-x_n||)
#'                  + c_0 + c_1*z_1 + ... + c_d*z_d  }
#' \cr
#' where \eqn{\Phi(r)=r^3} denotes the cubic radial basis function. The coefficients \eqn{\lambda_1, 
#' ..., \lambda_n, c_0, c_1, ..., c_d} are determined by this training procedure.
#' This is for the default case \code{squares==FALSE}. In case \code{squares==TRUE} 
#' there are d additional pure square terms and the model is
#' \cr
#'    \deqn{ s_{sq}(z) = s(z) + c_{d+1}*z_1^2 + ... + c_{d+d}*z_d^2 } 
#' \cr
#'   
#' The linear equation system is solved via SVD inversion. Near-zero elements 
#' in the diagonal matrix \eqn{D} are set to zero in \eqn{D^{-1}}. This is numerically stable 
#' for rank-deficient systems.
#' 
#' @param xp      n points \eqn{x_i} of dimension d are arranged in (n x d) matrix \code{xp}
#' @param U       vector of length n, containing samples \eqn{u(x_i)} of 
#'                the scalar function \eqn{u} to be fitted \cr
#'                - or - \cr
#'                (n x m) matrix, where each column 1,...,m contains one vector of samples
#'                \eqn{u_j(x_i)} for the m'th model, j=1,...,m
#' @param squares [FALSE] flag, see description
#' @param rho     [0.0] experimental: 0: interpolating, >0, approximating (spline-like) 
#'                Gaussian RBFs
#'                
#' @return \code{rbf.model},  an object of class \code{RBFinter}, which is basically a list 
#' with elements:
#'      \item{coef}{  (n+d+1 x m) matrix holding in column m the coefficients for the m'th 
#'                    model:      \eqn{\lambda_1, ..., \lambda_n, c_0, c_1, ..., c_d}.  
#'                    In case \code{squares==TRUE} it is an (n+2d+1 x m) matrix holding  
#'                    additionally the coefficients \eqn{c_{d+1}, ..., c_{d+d}}.}                    
#'      \item{xp}{  matrix xp   }
#'      \item{d}{  dimension d }
#'      \item{npts}{  number n of points \eqn{x_i} }
#'      \item{squares}{  TRUE or FALSE  }
#'      \item{type}{  "CUBIC"}
#'      
#' @seealso   \code{\link{trainGaussRBF}}, \code{\link{predict.RBFinter}}, \code{\link{interpRBF}}
#' @author Wolfgang Konen (\email{wolfgang.konen@@fh-koeln.de})
#----------------------------------------------------------------------------------
trainCubicRBF <- function(xp, U, squares=FALSE, rho=0.0) {
  ##DEBUG PK
  #save(xp, file="myXPvalue.RData")
  #save(U, file="myUvalue.RData")
  
  # calculate the inverse of M using SVD (numerically stable for rank-deficient M)
  svdInv <- function(M) {
    eps = 1e-14
    
    # alternate version PK
    #invM <- MASS::ginv(M)
    #return(invM)
    
    s = svd(M)
    invD <- 1/s$d
    invD[abs(s$d/s$d[1])<eps] <- 0
    invM <- s$v %*% diag(invD) %*% t(s$u)
    
    return(invM)
  }
  d=ncol(xp)
  if (squares) d=d+d
  npts=nrow(xp)
  edist=as.matrix(stats::dist(xp,upper=T))         # euclidean distance matrix
  phi=edist*edist*edist                     # cubic RBF (npts x npts matrix)
                                            # /WK/ experimental: rho>0 means spline approximation
  phi <- phi - diag(npts)*npts*rho          # /WK/ instead of exact interpolating RBF (rho=0)
  pMat=cbind(e=rep(1.0,npts),xp)            # linear tail LH(1,x1,x2,...)
  if (squares) pMat=cbind(pMat,xp*xp)       # ... plus direct squares x1^2, x2^2, ...
  nMat=matrix(0,d+1,d+1)                    # matrix of zeros
  
  M = rbind(cbind(phi,pMat),cbind(t(pMat),nMat))
  if (is.vector(U))
    rhs = c(U,rep(0,d+1))
  else if (is.matrix(U))
    rhs = rbind(U,matrix(0,d+1,ncol(U)))
  else 
    stop("U is neither vector nor matrix!")
  
  invM = svdInv(M)
  coef = invM %*% rhs
  ### 
  ### The linear-equation solution with 'solve' is more accurate when matrix M has full rank.
  ### But it crashes for rank-deficient or near-singular matrices M.
  ### This version only for vector-like input U:
  ###
  #coef2 = solve(M,rhs)
  #print(c(max(abs(M %*% coef - rhs)),max(abs(M %*% coef2 - rhs)),max(abs(rhs))))
  #browser()
  
  # this check does not yet work as expected, leave DEBUG==FALSE /WK/
  DEBUG=FALSE
  if (DEBUG) {
    inv2 = MASS::ginv(M)
    coef2 = inv2 %*% rhs
    testit::assert("rows and cols invM and inv2 do not match",nrow(invM)==nrow(inv2),ncol(invM)==ncol(inv2))
    browser()
    if (nrow(invM)<=ncol(invM)) 
      testit::assert("invM-check 1 failed",max(abs(invM %*% M - diag(rep(1,nrow(invM)))))<1e-5)
    else
      testit::assert("invM-check 2 failed",max(abs(M %*% invM - diag(rep(1,ncol(invM)))))<1e-5)    
  }
  
  rbf.model = list(coef=coef
                   ,xp=xp
                   ,d=d
                   ,npts=npts
                   ,squares=squares
                   ,type="CUBIC")  
  class(rbf.model) <- c("RBFinter","list")
  return(rbf.model)
}

#----------------------------------------------------------------------------------
#' Fit Gaussian RBF model to training data for d>1.
#' 
#' The model for a point \eqn{z=(z_1,...,z_d)} is fitted using n sample points \eqn{x_1, ..., x_n} 
#' \cr
#'    \deqn{ s(z) = \lambda_1*\Phi(||z-x_1||)+... +\lambda_n*\Phi(||z-x_n||)
#'                  + c_0 + c_1*z_1 + ... + c_d*z_d  }
#' \cr    
#' where \eqn{\Phi(r)=\exp(-r^2/(2*\sigma^2))} denotes the Gaussian radial basis function with width
#' \eqn{\sigma}. The coefficients \eqn{\lambda_1, ..., \lambda_n, c_0, c_1, ..., c_d} are determined 
#' by this training procedure.
#' This is for the default case \code{squares==FALSE}. In case \code{squares==TRUE} 
#' there are d additional pure square terms and the model is
#' \cr
#'    \deqn{ s_{sq}(z) = s(z) + c_{d+1}*z_1^2 + ... + c_{d+d}*z_d^2 } 
#'  
#' The linear equation system is solved via SVD inversion. Near-zero elements 
#' in the diagonal matrix \eqn{D} are set to zero in \eqn{D^{-1}}. This is numerically stable 
#' for rank-deficient systems.
#' 
#' @param xp      n points \eqn{x_i} of dimension d are arranged in (n x d) matrix \code{xp}
#' @param U       vector of length n, containing samples \eqn{u(x_i)} of 
#'                the scalar function \eqn{u} to be fitted \cr
#'                - or - \cr
#'                (n x m) matrix, where each column 1,...,m contains one vector of samples
#'                \eqn{u_j(x_i)} for the m'th model, j=1,...,m
#' @param squares [FALSE] flag, see description
#' @param RULE    ["One"] one out of ["One" | "Two" | "Three"], different rules for automatic width
#'                estimation.   
#' @param width   [-1] either a positive real value which is the constant width \eqn{\sigma} for all 
#'                Gaussians in all iterations, or -1. Then the appropriate width \eqn{\sigma} is 
#'                calculated anew in each iteration with one of the rules \code{RULE},
#'                based on the distribution of data points \code{xp}              
#' @param widthFactor [1.0] additional constant factor applied to each width \eqn{\sigma} 
#' @param rho     [0.0] experimental: 0: interpolating, >0, approximating (spline-like) 
#'                Gaussian RBFs
#'                
#' @return \code{rbf.model},  an object of class \code{RBFinter}, which is basically a list 
#' with elements:
#'      \item{coef}{  (n+d+1 x m) matrix holding in column m the coefficients for the m'th 
#'                    model:      \eqn{\lambda_1, ..., \lambda_n, c_0, c_1, ..., c_d}.  
#'                    In case \code{squares==TRUE} it is an (n+2d+1 x m) matrix holding  
#'                    additionally the coefficients \eqn{c_{d+1}, ..., c_{d+d}}.}
#'      \item{xp}{    matrix xp   }
#'      \item{d}{     dimension d }
#'      \item{npts}{  number n of points \eqn{x_i} }
#'      \item{squares}{TRUE or FALSE  }
#'      \item{width}{}
#'      \item{type}{  "GAUSS"}
#'      
#' @seealso   \code{\link{trainCubicRBF}}, \code{\link{predict.RBFinter}}, \code{\link{interpRBF}}
#' @author Wolfgang Konen, Samineh Bagheri
#'          
#----------------------------------------------------------------------------------
trainGaussRBF <- function(xp, U, width,squares=FALSE,RULE="One",widthFactor=1.0, rho=0.0) {
  if(width==-1){#When no width is given by user then width is selected automatically
    #RULE<-"One"
    #Automatic width adaptation is in testing phase therefore we would like to test several rules
    #from literature and then select teh most relevant ones
    
    #>> rule number One:
    #The rule is taken from Benoudjit,2002 (But originally from Haykin,1999)
    #It works based on  finding a compromise between locality and smoothness
    
    #>> rule number Two:
    #The rule is taken from Benoudjit,2002 (But originally from Moody and Darken,1989)
    #a vector is returned as width and width[i] is the width factor for the ith RBF
    
    #>> rule number Three:
    # This rule is taken from Sun and Jin,2013
    # The idea is basically adapting the width to the smallest interval of two points in each coordiante

    switch(RULE,
           "One" =     width<-(max(stats::dist(xp)))/sqrt(2*nrow(xp)),
           "Two" = {
              k<-2
              width<-sapply(1:nrow(xp),function(i){
                centroid<-matrix(xp[i,],nrow=1)
                ncentroid<-xp[-i,]
                #browser()
                #index<-FNN::knnx.index(ncentroid,centroid, k=k) #finding the index of k nearest neigbors
                mindists<-FNN::knnx.dist(ncentroid,centroid, k=k)
                #mindists<-mindists^2
                y<-(1/k)*(sqrt(sum(mindists*mindists)))
                return(y)
              })},
           "Three" = {width<-c()
                      for(i in 1:ncol(xp)){
                        interval<-max(xp[,i])-min(xp[,i])
                      width<-min(width,interval)  
             
           }}
           )
    #browser()
    verboseprint(verbose=0,important=FALSE,paste("Automatic adjustment of RBF width=",width[1]))
    
  }
  
  # calculate the inverse of M using SVD (numerically stable for rank-deficient M)
  svdInv <- function(M) {
    eps = 1e-14
    
    # alternate version PK
    #invM <- MASS::ginv(M)
    #return(invM)
    s = svd(M)
    invD <- 1/s$d
    invD[abs(s$d/s$d[1])<eps] <- 0
    invM <- s$v %*% diag(invD) %*% t(s$u)
    
    return(invM)
  }
  d=ncol(xp)
  if (squares) d=d+d
  npts=nrow(xp)
  edist=as.matrix(stats::dist(xp,upper=T))         # euclidean distance matrix
  #width=1
  width = width*widthFactor
  if (length(width)==1) {                   # /WK/ bug fix: former code was not correct for the case 
    wmat=width^2                            # that width may be a vector of length nrow(xp)
  } else {
    w = outer(rep(1,length(width)),width)
    wmat = w*w                              # this choice leads to non-symmetric phi, but it 
                                            # has smallest error on RbfFit2.R (multiply each width by 
                                            # factor 3 - 5)
    #wmat = w*t(w)                          # this is (wmat)_ij = width_i * width_j, it
                                            # guarantees that phi will be a symmetric matrix,
                                            # but it has bigger error  
  }
  phi <- exp(-0.5*(edist*edist)/wmat)       # Gaussian RBF with width=width
                                            # /WK/ experimental: rho>0 means spline approximation
  phi <- phi - diag(npts)*npts*rho          # /WK/ instead of exact interpolating RBF (rho=0)
  pMat=cbind(e=rep(1.0,npts),xp)            # linear tail LH(1,x1,x2,...)
  if (squares) pMat=cbind(pMat,xp*xp)       # ... plus direct squares x1^2, x2^2, ...
  nMat=matrix(0,d+1,d+1)                    # matrix of zeros

  
  M = rbind(cbind(phi,pMat),cbind(t(pMat),nMat))
  if (is.vector(U))
    rhs = c(U,rep(0,d+1))
  else if (is.matrix(U))
    rhs = rbind(U,matrix(0,d+1,ncol(U)))
  else 
    stop("U is neither vector nor matrix!")
  
  invM = svdInv(M)
  coef = invM %*% rhs
  ### 
  ### The linear-equation solution with 'solve' is more accurate when matrix M has full rank.
  ### But it crashes for rank-deficient or near-singular matrices M.
  ### This version only for vector-like input U:
  ###
  #coef2 = solve(M,rhs)
  #print(c(max(abs(M %*% coef - rhs)),max(abs(M %*% coef2 - rhs)),max(abs(rhs))))
  #browser()
  
  # this check does not yet work as expected, leave DEBUG==FALSE /WK/
  DEBUG=FALSE
  if (DEBUG) {
    inv2 = MASS::ginv(M)
    coef2 = inv2 %*% rhs
    testit::assert("rows and cols invM and inv2 do not match",nrow(invM)==nrow(inv2),ncol(invM)==ncol(inv2))
    browser()
    if (nrow(invM)<=ncol(invM)) 
      testit::assert("invM-check 1 failed",max(abs(invM %*% M - diag(rep(1,nrow(invM)))))<1e-5)
    else
      testit::assert("invM-check 2 failed",max(abs(M %*% invM - diag(rep(1,ncol(invM)))))<1e-5)    
  }
  
  rbf.model = list(coef=coef
                   ,xp=xp
                   ,d=d
                   ,npts=npts
                   ,squares=squares
                   ,width=width
                   ,type="GAUSS"
                   )  
  class(rbf.model) <- c("RBFinter","list")
  return(rbf.model)
}





#----------------------------------------------------------------------------------
#' Apply cubic or Gaussian RBF interpolation to new data for d>1.
#' 
#' @param x         vector holding a point of dimension d
#' @param rbf.model trained RBF model (or set of models), see \code{\link{trainCubicRBF}} 
#'                  or \code{\link{trainGaussRBF}}
#'                
#' @return          value \eqn{s(x)} of the trained model at x \cr
#'                  - or - \cr
#'                  vector \eqn{s_j(x)} with values for all trained models \eqn{j=1,...,m} at x
#' 
#' @seealso   \code{\link{trainCubicRBF}}, \code{\link{predict.RBFinter}}
#' @author Wolfgang Konen (\email{wolfgang.konen@@fh-koeln.de})
#----------------------------------------------------------------------------------
interpRBF <- function(x,rbf.model) {
  ##DEBUG PK
  #save(x, file="myXvaluePrediction.RData")
  #save(rbf.model, file="myRBFmodel.RData")
  
  #testit::assert("non-conform",length(x)==ncol(rbf.model$xp))
  if (length(x)!=ncol(rbf.model$xp)) {
    cat("problem in interpRBF\n")
    #browser()
    stop("Problem in interpRBF")
  }
  ed = distLine(x,rbf.model$xp)      # euclidean distance of x to all xp, 
                                     # this is up to 40x faster than dist() !!
                                     # ed is vector of length nrow(xp)  
  switch(rbf.model$type,
         "CUBIC" = {
           ph = ed*ed*ed           
         },
         "GAUSS" = {
           ph = exp(-0.5*(ed/rbf.model$width)^2)  # this works correctly even in case where width is a vector of length nrow(xp)
         })
    
  if (rbf.model$squares) 
    lhs = c(ph,1,x,x*x)
  else 
    lhs = c(ph,1,x)
  val = as.vector(lhs %*% rbf.model$coef)
  return (val)
}

#----------------------------------------------------------------------------------
#' Apply cubic or Gaussian RBF interpolation to new data for d>1.
#' 
#' Same as \code{\link{interpRBF}}, but it uses \code{dist} instead of 
#' \code{distLine} and it is 40x slower on data 
#' d=124, n_model=69, npts=249, see \code{testRbfInter2.R}
#' 
#' @param x         vector holding a point of dimension d
#' @param rbf.model trained RBF model (or set of models), see \code{\link{trainCubicRBF}} 
#'                  or \code{\link{trainGaussRBF}}
#'                  
#' @return val      value s(x) of the trained model at x 
#'                  - or - 
#'                  vector \eqn{s_j(x)} with values for all trained models \eqn{j=1,...,m} at x
#' 
#' @seealso   \code{\link{trainCubicRBF}}, \code{\link{predict.RBFinter}}
#' @keywords internal   
#' @author Wolfgang Konen (\email{wolfgang.konen@@fh-koeln.de})
#----------------------------------------------------------------------------------
interpRBF_old <- function(x,rbf.model) {
  ed = stats::dist(rbind(x,rbf.model$xp))    
  ed = ed[1:rbf.model$npts]    
  ph = ed*ed*ed
  if (rbf.model$squares) 
    lhs = c(ph,1,x,x*x)
  else 
    lhs = c(ph,1,x)
  val = as.vector(lhs %*% rbf.model$coef)
  return (val)
}

#----------------------------------------------------------------------------------
#' Apply cubic or Gaussian RBF interpolation
#' 
#' Apply cubic or Gaussian RBF interpolation to a set of new data points for d>1.
#' 
#' @param rbf.model trained RBF model (or set of models), see \code{\link{trainCubicRBF}} 
#'                  or \code{\link{trainGaussRBF}}
#' @param newdata   matrix or data frame with d columns. Each row contains a data point 
#'                  \eqn{x_i,\ i=1,\ldots,n}
#' @param ...       (not used)
#'                
#' @return          vector of model responses \eqn{s(x_i)}, one element for each data point \eqn{x_i} \cr
#'                  - or - \cr
#'                  if \code{rbf.model} is a set of \code{m} models, a \code{(n x m)}-matrix 
#'                  containing in each row the response \eqn{s_j(x_i)} of all models 
#'                  \eqn{j = 1,\ldots,m}  to \eqn{x_i}
#'  
#' @seealso   \code{\link{trainCubicRBF}}, \code{\link{trainGaussRBF}}, \code{\link{interpRBF}}
#' @author Wolfgang Konen (\email{wolfgang.konen@@fh-koeln.de})
#----------------------------------------------------------------------------------
predict.RBFinter <- function(rbf.model,newdata,...) {
  val = t(sapply(1:nrow(newdata),function(i)interpRBF(as.numeric(newdata[i,]),rbf.model)))
  # sapply binds the result vectors from interpRBF() by default **column-wise** together.
  # Therefore we use t() to return in the ith *row* of val 
  # the results for the ith row in newdata (row-wise).
 
#### 
#### This alternative version, which does the job of interpRBF right here (in-place),
#### is not really faster than the one sapply-line above (!)
####
#   x = as.matrix(newdata)
#   nx = nrow(x)
#   ed = t(sapply(1:nx,function(i)(dist(rbind(x[i,],rbf.model$xp)))[1:rbf.model$npts] ))
#   ph = ed*ed*ed
#   if (rbf.model$squares) 
#     lhs = cbind(ph,rep(1,nx),x,x*x)
#   else 
#     lhs = cbind(ph,rep(1,nx),x)
#   val = as.vector(lhs %*% rbf.model$coef)
  
  # if there is only one model in rbf.model:
  if (ncol(rbf.model$coef)==1) val = as.vector(val)
  return(val)
}
