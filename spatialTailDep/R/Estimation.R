#' @include gridUtilities.r
#' @include tailIntEmp.r
#' @include EstimationSmith.r
#' @include EstimationBR.r
NULL

#' Function \code{tailInt}
#' 
#' Integral of the bivariate parametric stable tail dependence function over the unit square, for the Smith model or the Brown-Resnick process.
#' 
#' @param loc A 2 x 2 matrix, where a row represents a location.
#' @param model Choose between "smith" and "BR". 
#' @param theta Parameter vector. For the Smith model, \code{theta} must be equal to the 2 x 2 covariance matrix. For the Brown-Resnick pocess, \code{theta} \eqn{ = (\alpha, \rho, \beta, c)}.
#' @return A scalar.
#' @export
#' @details This is an analytic implementation of the integral of the stable tail dependence function,
#' which is much faster than numerical integration. For the definitions of the parametric stable tail dependence 
#' functions, see Einmahl et al. (2014).
#' 
#' The parameter vector \code{theta} must be a positive semi-definite matrix if \code{model = "smith"} 
#' and a vector of length four if \code{model = "BR"}, where \eqn{0 < \alpha < 1}, \eqn{\rho > 0},
#' \eqn{0 < \beta \le \pi/2} and \eqn{c > 0}.
#' @references Einmahl, J.H.J., Kiriliouk, A., Krajina, A. and Segers, J. (2014), "An M-estimator of spatial tail dependence". See \url{http://arxiv.org/abs/1403.1975}. 
#' @seealso \code{\link{Mestimator}}, \code{\link{tailIntEmp}}
#' @examples
#' tailInt(loc = cbind(c(1,1),c(2,3)), model = "smith", theta = rbind(c(3,1),c(1,2)))
#' tailInt(loc = cbind(c(1,2),c(3,4)), model = "BR", theta = c(1.5,1,0.5,0.75))

tailInt <- function(loc, model, theta){
  if(model == "smith"){
    if(any(eigen(theta)$values<0)){
        warning("The covariance matix theta isn't positive semi-definite")
    } else{
        return(tailSmithInt(loc = c(t(loc)),siginv = solve(theta)))
    }  
  } else if(model == "BR"){
    if(theta[1] <= 0 || theta[1] >=2 || theta[2] <=0 || theta[3]<0 || theta[3]>=(pi/2) || theta[4] <=0){
      warning("Parameter values not valid")
    } else{
      alpha<-theta[1]
      taum<-parsBRtoTau(theta[2:4])
      return(tailBRInt(loc = c(t(loc)),tau=rbind(c(taum[1],taum[3]),c(taum[3],taum[2])),alpha=alpha))
    }
  }
}

tailFunc <- function(t, loc, model, theta){ 
  d<-nrow(loc)
  if(length(t) != d){
    warning("t must be a vector whose length is equal to the number of locations")
  } else{
      if(model == "smith"){
        if(any(eigen(theta)$values<0)){
          warning("The covariance matix theta isn't positive semi-definite")
        } else{
          return(tailSmith(t,t(loc),solve(theta),d))
        }
      } else if(model =="BR"){
        if(theta[1] <= 0 || theta[1] >=2 || theta[2] <=0 || theta[3]<0 || theta[3]>=(pi/2) || theta[4] <=0){
          warning("Parameter values not valid")
        } else{
          alpha<-theta[1]
          taum<-parsBRtoTau(theta[2:4])
          return(tailBR(t,t(loc),rbind(c(taum[1],taum[3]),c(taum[3],taum[2])),alpha,d))
      }
     }
  }    
}

#' Function \code{AsymVar}
#' 
#' Function to compute the asymptotic variance matrix of the pairwise M-estimator for the Smith model or the Brown-Resnick process.
#' 
#' @param pairs A \eqn{q} x 4 matrix giving the Cartesian coordinates of the \eqn{q} pairs of locations.
#' @param model Choose between "smith" and "BR". 
#' @param theta Parameter vector. For the Smith model, \code{theta} must be equal to the 2 x 2 covariance matrix. For the Brown-Resnick pocess, \code{theta} \eqn{ = (\alpha, \rho, \beta, c)}.
#' @param Tol The tolerance in the numerical integration procedure. Defaults to 1e-5.
#' @references Einmahl, J.H.J., Kiriliouk, A., Krajina, A. and Segers, J. (2014), "An M-estimator of spatial tail dependence". See \url{http://arxiv.org/abs/1403.1975}. 
#' @return A \eqn{q} x \eqn{q} matrix.
#' @seealso \code{\link{Mestimator}}, \code{\link{selectPairIndices}}, \code{\link{pairCoordinates}}
#' @details For a matrix of coordinates of pairs of locations, this function returns the asymptotic
#' variance matrix of the estimator. An optimal weight matrix can be defined as the inverse of the
#' asymptotic variance matrix. For a detailed description of this procedure, see Einmahl et al. (2014).
#' 
#' The parameter vector \code{theta} must be a positive semi-definite matrix if \code{model = "smith"} 
#' and a vector of length four if \code{model = "BR"}, where \eqn{0 < \alpha < 1}, \eqn{\rho > 0},
#' \eqn{0 < \beta \le \pi/2} and \eqn{c > 0}.
#' @export
#' @examples
#' ## Define the locations of three stations
#' (locations <- rbind(c(1.0,1.0),c(2.0,1.0),c(1.2,2.5)))
#' ## select pairs
#' (pairIndices <- selectPairIndices(locations, maxDistance = 3))
#' (pairs <- pairCoordinates(locations, pairIndices))

#' ## Smith model parameter matrix
#' (theta <- rbind(c(1.5, .5), c(.5, 1)))
#' ## The matrix. Takes a couple of seconds to compute.
#' ## AsymVar(pairs, model = "smith", theta = theta, Tol = 1e-04)
#' 
#' ## Parameters of the Brown-Resnick process
#' (theta <- c(1.5,1,0.5,0.25))
#' ## The matrix. Takes a couple of seconds to compute.
#' ## AsymVar(pairs, model = "BR", theta = theta, Tol = 1e-04)
AsymVar<-function(pairs, model, theta, Tol = 1e-5){
  if(model == "smith"){
    if(any(eigen(theta)$values<0)){
      warning("The covariance matix theta isn't positive semi-definite")
    } else{
        return(AsymVarSm(pairs,solve(theta),Tol))
      }
    } else if(model =="BR"){
      if(theta[1] <= 0 || theta[1] >=2 || theta[2] <=0 || theta[3]<0 || theta[3]>=(pi/2) || theta[4] <=0){
        warning("Parameter values not valid")
      } else{
          alpha<-theta[1]
          taum<-parsBRtoTau(theta[2:4])
          return(AsymVarBR(pairs,rbind(c(taum[1],taum[3]),c(taum[3],taum[2])),alpha,Tol))
      }  
  }
}

#' Function \code{Mestimator}
#' 
#' Function to compute the pairwise M-estimator for the parameters of the Smith model or the Brown-Resnick process.
#' 
#' @param x An \eqn{n} x \eqn{d} data matrix.
#' @param locations A \eqn{d} x 2 matrix containing the Cartesian coordinates of \eqn{d} points in the plane.
#' @param pairIndices A \eqn{q} x 2 matrix containing the indices of \eqn{q} pairs of points from the matrix \code{locations}.
#' @param k The threshold parameter in the definition of the empirical stable tail dependence function. 
#' @param model Choose between "smith" and "BR".
#' @param Tol The tolerance parameter in the numerical integration procedure; defaults to 1e-05.
#' @param startingValue Initial value of the parameters in the minimization routine. Defaults to diag(2) for the Smith model and (1, 1.5, 0.75, 0.75) for the BR process.
#' @param Omega A \eqn{q} x \eqn{q} matrix specifying the metric with which the distance between the parametric and nonparametric estimates will be computed. The default is the identity matrix, i.e., the Euclidean metric.
#' @param iterate A Boolean variable. If \code{TRUE} (the default), then the estimator is calculated twice, first with \code{Omega} specified by the user, and then a second time with the optimal \code{Omega} calculated at the initial estimate.
#' @param covMat A Boolean variable. If \code{TRUE} (the default), the covariance matrix is calculated.
#' @return A list with the following components:
#' \tabular{ll}{
#' \code{theta} \tab The estimator with estimated optimal weight matrix. \cr
#' \code{theta_pilot} \tab The estimator without the optimal weight matrix.\cr
#' \code{covMatrix} \tab The estimated covariance matrix for the estimator. \cr
#' \code{Omega} \tab The weight matrix with which the estimator was calculated.
#' }
#' @details  For a detailed description of the estimation procedure, see Einmahl et al. (2014).
#' Some tips for using this function: 
#' \itemize{
#' \item{\code{n} versus \code{d}: }{ if the number of locations \eqn{d} is small (\eqn{d < 8} say), a sufficiently
#' large sample size (eg \eqn{n > 2000}) is needed to obtain a satisfying result, especially for the Brown-Resnick process. 
#' However, if \eqn{d} is large, a sample size of \eqn{n = 500} should suffice.}
#' \item{\code{pairIndices}: }{ if the number of pairs \eqn{q} is large, \code{Mestimator} will be rather slow. This is due
#' to the calculation of \code{Omega} and \code{covMat}. Setting \code{iterate = FALSE} and 
#' \code{covMat = FALSE} will make this procedure fast even for several hundreds of pairs of locations.}
#' \item{\code{Tol}: }{ the tolerance parameter is used when calculating the three- and four-dimensional integrals
#' in the asymptotic covariance matrix (see Appendix B in Einmahl et al. (2014)). A tolerance of 1e-04 often suffices, although
#' the default tolerance is a safer choice.}
#' \item{\code{StartingValue}: }{ for the Smith model, the estimator usually doesn't depend on the starting value
#' at all. For the Brown-Resnick process, it is advised to try a couple of starting values if \eqn{d} 
#' is very small, preferably a starting value with \eqn{c < 1} and one with \eqn{c > 1}.}
#' \item{\code{iterate}: }{ if \code{iterate = TRUE}, the matrix \code{Omega} is calculated. This weight matrix tends to have a larger
#' effect when \eqn{d} is large and/or when the Smith model is used.}
#' \item{\code{covMat}: }{ if the resulting covariance matrix is incorrect (eg negative diagonal values), then \code{Tol} is set too high.
#' For the Smith model, the order of the parameters is \eqn{(\sigma_{11},\sigma_{22},\sigma_{12})}. }
#' }
#' @seealso \code{\link{selectPairIndices}}, \code{\link{pairCoordinates}}
#' @references Einmahl, J.H.J., Kiriliouk, A., Krajina, A. and Segers, J. (2014), "An M-estimator of spatial tail dependence". See \url{http://arxiv.org/abs/1403.1975}. 
#' @export
#' @examples
#' ## define the locations of 4 stations
#' (locations <- rbind(c(1,1),c(2,1),c(1,2),c(2,2)))
#' ## select the pairs of locations; here, we select all locations
#' (pairIndices <- selectPairIndices(locations, maxDistance = 2))
#' 
#' ## We use the rmaxstab function from the package SpatialExtremes to 
#' ## simulate from the Smith and the Brown-Resnick process.
#' 
#' ## The Smith model
#' set.seed(2)
#' x<-rmaxstab(n = 5000, coord = locations,cov.mod="gauss",cov11=1,cov22=2,cov12=0.5)
#' ## calculate the pairwise M-estimator. This may take up to one minute or longer.
#' ## Mestimator(x, locations, pairIndices, 100, model="smith",Tol = 5e-04)
#'
#' ## The Brown-Resnick process
#' set.seed(2)
#' x <- rmaxstab(n = 5000, coord = locations, cov.mod = "brown", range = 3, smooth = 1)
#' ## We can only simulate isotropic processes with rmaxstab, so we multiply the coordinates
#' ## of the locations with V^(-1) (beta,c). Here we choose beta = 0.25 and c = 1.5
#' (Vmat<-matrix(c(cos(0.25),1.5*sin(0.25),-sin(0.25),1.5*cos(0.25)),nrow=2))
#' (locationsAniso <- locations %*% t(solve(Vmat)))
#' ## calculate the pairwise M-estimator. This may take up to one minute or longer.
#' ## Mestimator(x, locationsAniso, pairIndices, 300, model="BR",Tol = 5e-04)

Mestimator <- function(x, locations, pairIndices, k, model, Tol = 1e-05, startingValue = NULL, 
                            Omega = diag(nrow(pairIndices)), iterate = TRUE, covMat=TRUE) {
  if (length(k) > 1) {
    k <- k[1]
    warning("Currently, only scalar k has been implemented. Keeping k[1] and ignoring the rest.")
  }
  if(model == "smith"){
    if(is.null(startingValue)){startingValue<-diag(2)}
    return(MestimatorSmith(x, locations, pairIndices,k, Tol=Tol,startingValue = startingValue,Omega=Omega,iterate=iterate,covMat=covMat))
  } else if(model == "BR"){
    if(is.null(startingValue)){startingValue<-c(1,1.5,0.75,0.75)}
    return(MestimatorBR(x, locations, pairIndices,k, Tol=Tol,startingValue = startingValue,Omega=Omega,iterate=iterate,covMat=covMat))
  }
}
