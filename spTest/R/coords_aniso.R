#' Geometric Anisotropy Correction from geoR
#'
#' This is the function \code{\link[geoR]{coords.aniso}} from the R package \strong{geoR}. It performs the linear transformation or reverse transformation of a set of spatial coordinates in \eqn{R^2} according to geometric anisotropy parameters. See \href{http://cran.r-project.org/web/packages/geoR/geoR.pdf}{geoR documentation} for more details.
#'
#'@export
#'@keywords external
#'
#' @param coords An \eqn{n \times 2} matrix of spatial coordinates.
#' @param aniso.pars A vector of length two containing the anisotropy angle and anisotropy ratio, respectively.
#' @param reverse Logical. If \code{reverse = TRUE} the reverse transformation is performed.
#'
#' @details Details on the function can be found in the \href{http://cran.r-project.org/web/packages/geoR/geoR.pdf}{geoR documentation}. The function is included in this package to avoid the necesity of loading the \strong{geoR} package when simulating anisotropic spatial processes.
#' 
#'
#'@return An \eqn{n \times 2} matrix of transformed coordinates.
#'
#'@author Paulo Justiniano Ribeiro Jr.
#'@author Peter J. Diggle
#'
#'@references   Paulo J. Ribeiro Jr & Peter J. Diggle geoR: a package for geostatistical analysis R-NEWS, 1(2):15-18. June, 2001.
#'
#' @seealso \code{\link[geoR]{coords.aniso}}
#' @examples
#' #Set parameter values for exponential covariance function
#' sigma.sq <- 1
#' tau.sq <- 0.0
#' phi <- 1/4
#' #Set anisotropy parameters
#' aniso.angle <- pi/4
#' aniso.ratio <- 2
#' #function to compute distance corresponding to a correlation contour 
#' #for a spatial process in R^2 with an exponential covariance function
#' cor.dist = function(cv, ss, ts, phi)
#' {
#'	h <- -1*( log(cv*(ss+ts)/ss) ) /(phi)
#'	return(h)
#' }
#' #compute the distance of 0.5 correlation for isotropic process
#' r <-  cor.dist(0.5, sigma.sq, tau.sq, phi)
#' #create a sequence of points corresponding to equicorrelation 
#' #with a data point observed at the location (0,0)
#' angle <- seq(0, 2*pi, by = 0.01)
#' x <- r*sin(angle)
#' y <- r*cos(angle)
#' cor.coords <- cbind(x,y)
#' #plot the contour of equicorrelation for isotropic process
#' plot(cor.coords, type = "l", xlim = c(-5,5), ylim = c(-5,5), pty = "s", lwd = 2)
#' points(0,0, pch = 20, cex = 1.3)
#' #transform and plot the coordinates according to the anisotropy parameters
#' aniso.cor.coords <- coords.aniso(cor.coords, c(aniso.angle, aniso.ratio), rev = TRUE)
#' lines(aniso.cor.coords, lwd = 2, lty = 2)

coords.aniso = function (coords, aniso.pars, reverse = FALSE) 
{
    coords <- as.matrix(coords)
    n <- nrow(coords)
    if (length(aniso.pars) != 2) 
        stop("argument aniso.pars must be a vector with 2 elements: the anisotropy angle and anisotropy ratio, respectively")
    psiA <- aniso.pars[1]
    psiR <- aniso.pars[2]
    if (psiR < 1) {
        psiR <- round(psiR, digits = 8)
        if (psiR < 1) 
            stop("anisotropy ratio must be greater than 1")
    }
    rm <- matrix(c(cos(psiA), -sin(psiA), sin(psiA), cos(psiA)), 
        ncol = 2)
    tm <- diag(c(1, 1/psiR))
    if (reverse) 
        coords.mod <- coords %*% solve(rm %*% tm)
    else coords.mod <- coords %*% rm %*% tm
    return(coords.mod)
}