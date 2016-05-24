outlier <- function(x,k=3,location=mean,scale=sd)
    {
    LOCATION <- location(x)
    SCALE <- scale(x)
    OUTLIER <- ifelse(abs(x-LOCATION) >= k*SCALE, 1, 0)
    OUTLIER
    }




#' Function to Find Outliers for an epplab Object
#' 
#' Function to decide wether observations are considered outliers or not in
#' specific projection directions of an \code{epplab} object.
#' 
#' Denote \eqn{location_j} as the location of the jth projection direction and
#' analogously \eqn{scale_j} as its scale. Then an observation \eqn{x} is an
#' outlier in the jth projection direction, if \eqn{|x-location_j| \geq k \
#' scale_j}{|x-location_j| >= k scale_j}.
#' 
#' Naturally it is best to use for this purpose robust location and scale
#' measures like \code{\link{median}} and \code{\link{mad}} for example.
#' 
#' @param x An object of class \code{epplab}.
#' @param which The directions in which outliers should be searched. The
#' default is to look at all.
#' @param k Numeric value to decide when an observation is considered an
#' outlier or not. Default is 3. See details.
#' @param location A function which gives the univariate location as an output.
#' The default is \code{\link{mean}}.
#' @param scale A function which gives the univariate scale as an output. The
#' default is \code{\link{sd}}.
#' @return A list with class 'epplabOutlier' containing the following
#' components: \item{outlier}{A matrix with only zeros and ones. A value of 1
#' classifies the observation as an outlier in this projection direction.}
#' \item{k}{The factor \code{k} used.} \item{location}{The name of the
#' \code{location} estimator used.} \item{scale}{The name of the \code{scale}
#' estimator used.} \item{PPindex}{The name of the \code{PPindex} used.}
#' \item{PPalg}{The name of the \code{PPalg} used.}
#' @author Klaus Nordhausen
#' @seealso \code{\link{EPPlab}}
#' @references \cite{Ruiz-Gazen, A., Larabi Marie-Sainte, S. and Berro, A.
#' (2010), Detecting multivariate outliers using projection pursuit with
#' particle swarm optimization, \emph{COMPSTAT2010}, pp. 89-98.}
#' @keywords multivariate
#' @examples
#' 
#' # creating data with 3 outliers
#' n <-300 
#' p <- 10
#' X <- matrix(rnorm(n*p),ncol=p)
#' X[1,1] <- 9
#' X[2,4] <- 7 
#' X[3,6] <- 8
#' # giving the data rownames, obs.1, obs.2 and obs.3 are the outliers.
#' rownames(X) <- paste("obs",1:n,sep=".")
#' 
#' PP<-EPPlab(X,PPalg="PSO",PPindex="KurtosisMax",n.simu=20, maxiter=20)
#' OUT<-EPPlabOutlier(PP, k = 3, location = median, scale = mad)
#' OUT
#' 
#' @export EPPlabOutlier
EPPlabOutlier <- function(x, which=1:ncol(x$PPdir), k=3, location=mean, scale=sd)
    {
    LOCname <- deparse(substitute(location))
    SCAname <- deparse(substitute(scale))
    
    FITS <- fitted(x, which=which)
    OUTLIER <- apply(FITS,2,outlier,k=k,location=location, scale=scale)
    RES <- list(outlier=OUTLIER, k=k, location=LOCname, scale=SCAname, PPindex=x$PPindex, PPalg=x$PPalg)
    class(RES) <- "epplabOutlier"
    RES
    }
