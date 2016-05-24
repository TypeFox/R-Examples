#' Maximum Entropy [De]Regularized Density Estimation
#' 
#' Density estimation based on maximum entropy methods
#' 
#' See the references for further details. And also Mosek "Manuals". The
#' acronym has a nice connection to
#' 
#' http://www.urbandictionary.com/define.php?term=medder
#' 
#' Term used in Bahamian dialect, mostly on the Family Islands like Eleuthera
#' and Cat Island meaning "mess with" "get involved" "get entangled" "fool
#' around" "bother"
#' 
#' "I don't like to medder up with all kinda people"
#' 
#' "Don't medder with people (chirren)"
#' 
#' "Why you think she medderin up in their business."
#' 
#' This version implements a class of penalized density estimators solving:
#
#' \deqn{\min_x \phi(x_1) | A_1 x_1 - A_2 x_2 = b,  0 \leq x_1, -\lambda \leq x_2 \leq \lambda }
#
#' where \eqn{x}x is a vector with two component subvectors: \eqn{x_1} is a 
#' vector of function values of the density \eqn{x_2} is a vector of dual values,
#' \eqn{\lambda} is typically positive, if negative then the \eqn{x_2} constraint
#' is replaced by \eqn{x_2 \geq 0}, which for \eqn{\alpha = 1}, 
#' constrains the fitted density to be log-concave, for
#' negative lambda fitting for other alphas and Dorders 
#' proceed at your own risk. See also below...
#' \eqn{\phi} is an additive convex function in the coordinates of \eqn{x_1},
#' interpretable as a negative Renyi entropy.
#' This formulation purports to solve a class of dual penalized maximum (Renyi) entropy
#' problems in which regularization is achieved by constraining the sup norm of the
#' dual vector.  In the primal representation of the problems this corresponds to
#' a roughness penalty on the total variation of the Dorder derivative of some
#' transformation of the fitted density.  
#'
#' @param x Data: either univariate or bivariate (not yet implemented in
#' Rmosek)
#' @param v Undata: either univariate or bivariate, by default there is an
#' equally spaced grid of 300 values
#' @param lambda total variation penalty parameter, if lambda is in [-1,0], a
#' concavity constraint is imposed. If lambda is in (-oo, -1) a convexity
#' constraint on .5 x^2 + log f is imposed. See Koenker and Mizera (2013) for
#' further details on this last option, and Koenker and Mizera (2010) for
#' further details on the concavity constrained options.
#' @param alpha Renyi entropy parameter characterizing fidelity criterion
#' @param Dorder Order of the derivative operator for the penalty
#' @param rtol Convergence tolerance for Mosek algorithm
#' @param verb Parameter controlling verbosity of solution, 0 for silent, 5
#' gives rather detailed iteration log.
#' @param control Mosek control list see KWDual documentation
#' @return An object of class "Medde" with components \item{x}{points of
#' evaluation on the domain of the density} \item{y}{estimated function values
#' at the evaluation points x} \item{phi}{n by 2 matrix of (sorted) original
#' data and estimated log density values these data points} \item{logLik}{log
#' likelihood value} \item{status}{exit status from Mosek}
#' @author Roger Koenker and Ivan Mizera
#' @seealso This function is based on an earlier function of the same name in
#' the deprecated package MeddeR that was based on an R-Matlab interface.
#' @references  Chen, Y. and R.J. Samworth, (2013) "Smoothed log-concave
#' maximum likelihood estimation with applications", \emph{Statistica Sinica},
#' 23, 1373--1398.
#' 
#' Koenker, R and I. Mizera, (2007) ``Density Estimation by Total Variation
#' Regularization,'' \emph{Advances in Statistical Modeling and Inference:
#' Essays in Honor of Kjell Doksum}, V.N. Nair (ed.), 613-634.
#' 
#' Koenker, R and I. Mizera, (2006) ``The alter egos of the regularized maximum
#' likelihood density estimators: deregularized maximum-entropy, Shannon,
#' Renyi, Simpson, Gini, and stretched strings,'' \emph{ Proceedings of the 7th
#' Prague Symposium on Asymptotic Statistics}.
#' 
#' Koenker, R and I. Mizera, (2010) ``Quasi-Concave Density Estimation''
#' \emph{Annals of Statistics}, 38, 2998-3027.
#' 
#' Koenker, R and I. Mizera, (2013) ``Convex Optimization, Shape Constraints,
#' Compound Decisions, and Empirical Bayes Rules,'' JASA, 109, 674--685.
#' 
#' Koenker, R and I. Mizera, (2014) ``Convex Optimization in R.'',
#' \emph{Journal of Statistical Software}, 60, 1-23.
#' @keywords nonparametric
#' @export
#' @import Matrix
#' @examples
#' 
#' #Maximum Likelihood Estimation of a Log-Concave Density
#' set.seed(1968)
#' x <- rgamma(50,10)
#' m <- medde(x, v = 50, lambda = -.5, verb = 5)
#' plot(m, type = "l")
#' lines(m$x,dgamma(m$x,10),col = 2)
#' points(x,m$g,cex = 0.5)
#' rug(x)
#' title(paste("log likelihood = ", round(m$logLik,2)))
#' legend(14,.12,c("ghat","true"),lty = 1, col = 1:2)
#' 
#' #Maximum Likelihood Estimation of a Gamma Density with TV constraint
#' set.seed(1968)
#' x <- rgamma(50,5)
#' f <- medde(x, v = 50, lambda = 0.005, verb = 5)
#' plot(f, type = "l")
#' lines(f$x,dgamma(f$x,5),col = 2)
#' legend(10,.15,c("ghat","true"),lty = 1, col = 1:2)
#' 
#' 
medde <- function(x, v = 300, lambda = 0.5, alpha = 1, Dorder = 1, 
	rtol = 1e-06, verb = 0, control = NULL){
############################################################################
#
#
# First version: 22 Dec 2006  (only does Dorder <- 1, and a few alphas)
#               26 Dec 2006  (does Dorder <- 1,2,3, and a few more alphas)
#               28 Dec 2006  (added merge flag) 
#               31 Dec 2006  (log-concave case) 
#                5 Jan 2007  (bivariate case) 
#                9 Apr 2008  (well-tempered case) 
#                1 Jul 2015  (Revised version without SparseM)
############################################################################


if(length(v) == 1){
	eps <- ifelse(lambda < 0, 0.0001,1)
	v <- seq(min(x) - eps, max(x) + eps, length = v)
	}
dimx <- NCOL(x)
dimv <- NCOL(v)

mesh1 <- function(x, v, Dorder = 1) {
############################################################################
# mesh <- mesh1(x, v, Dorder, merge)
#
# Meshing for univariate medde, old option to merge x and v removed in R version
#
# Input:
#    x      - vector of n observations (aka data)
#    v      - vector of m virtual observations (aka undata)
#    Dorder - order of the TV penalty, e.g. 0 for TV(phi(f)), 1 for TV(phi(f)'), ... 
#
# Output:  
#
#      XL  - the evaluation operator 
#      XC  - the Riemann weights used to integrate 
#      XJ  - the penalty block of the constraint matrix
#      XU  - the unique values at which the density is to be estimated
#
# Roger Koenker, last MeddeR revision 4 Jan 2007
#                    RMosek version 12 Nov 2012
############################################################################

n <- length(x)
p <- length(v)
k <- findInterval(x,v)
h <- diff(v)
ia <- c(1:n,1:n)
ja <- c(k+1,k)
ra <- c((x - v[k])/h[k], (v[k+1] - x)/h[k])
R <- t(sparseMatrix(ia, ja, x = ra, dims = c(n,p)))
s <- 1/h
d <- 0.5*(c(h, 0) + c(0,h))
H <- Diagonal(x = d)

if(Dorder %in% 0:2){
  switch(Dorder + 1,
  {#case 0
     q <-  p-1
     IA <- c(1:(p-1), 1:(p-1))
     JA <- c(1:(p-1), 2:p)
     XA <- c(-s[1:(p-1)], s[1:(p-1)])
     },
  {#case 1
     q <-  p-2
     IA <- c(1:(p-2), 1:(p-2), 1:(p-2))
     JA <- c(1:(p-2), 2:(p-1), 3:p)
     XA <- c(-s[1:(p-2)], s[1:(p-2)]+s[2:(p-1)], -s[2:(p-1)])
     D  <- .5 * (h[1:(p-2)]+h[2:(p-1)])
     XA <- XA / c(D, D, D)
     },
  {#case 2
     q <-  p-3
     IA <- c(1:(p-3), 1:(p-3), 1:(p-3), 1:(p-3))
     JA <- c(1:(p-3), 2:(p-2), 3:(p-1), 4:p)
     D1 <-  .5 * (h[1:(p-3)] + h[2:(p-2)])
     D2 <-  .5 * (h[2:(p-2)] + h[3:(p-1)])
     DA <- .5 * (D1 + D2)
     XA <- c(s[1:(p-3)]/D1, - (s[1:(p-3)] +  s[2:(p-2)])/D1 - s[2:(p-2)]/D2,
             s[2:(p-2)]/D1 + s[2:(p-2)]/D2 + s[3:(p-1)]/D2, -s[3:(p-1)]/D2)
     XA <- XA/c(DA, DA, DA, DA)
     }
  )
}
else
   stop("Dorder must be in {0,1,2}")

   XJ <- sparseMatrix(IA, JA, x = XA, dims = c(q,p))
   list(XJ = XJ, XC = d, XL = R, XU = v)
}
if(dimx != dimv) 
	stop('x and v of different dimensions')
if(dimx == 1)
	mesh <- mesh1(x,v,Dorder)
else if(dimx == 2){
   if(Dorder != 1)
	stop('Dorder must be 1 for bivariate data')
   else {
	stop("Bivariate smoothing not (yet) implemented for REBayes")
	#mesh <- mesh2(x,v,merge)
	#T <- mesh$T
        }
   }
else
   stop('x and v must be either 1d or 2d')


XL <- mesh$XL
XC <- mesh$XC
XJ <- mesh$XJ
XU <- mesh$XU

n <- NROW(x)
p <- NCOL(XJ)
q <- NROW(XJ)
H <- Diagonal(x = XC)
C <- rep(0,p+q)
if(lambda > -1)
   A <- cbind(H , t(XJ))
else {
   A <- cbind(H , -t(XJ))
   C[(p+1):(p+q)] <- rep(1,q)
   }
L <- as(XL %*% rep(1,n)/n, "vector")

if(lambda > 0) { # TV Constraint
   LX <- c(rep(0,p),  -rep(lambda,q))
   UX <- c(rep(Inf,p), rep(lambda,q))
   }
else {  # Concavity/Convexity Constraint
   LX <- rep(0,p+q)
   UX <- rep(Inf,p+q)
   }

P <- list(sense = "min")
P$c <- C
P$A <- A
P$bx <- rbind(LX,UX)
P$bc <- rbind(L,L)

opro <- matrix(list(),nrow = 5, ncol = p)
rownames(opro) <- c(" type ", "j", "f", "g", "h")

#switch solver case{'mosek'} # MOSEK solution  (For another day, only mosek for now)
if(alpha == 1){ # Shannon
    opro[1, ] <- as.list(rep("ent", p))
    opro[2, ] <- as.list(1:p)
    opro[3, ] <- as.list(XC)
    opro[4, ] <- as.list(rep(0, p))
    opro[5, ] <- as.list(rep(0, p))
   }
else if(alpha == 0.5){ # Hellinger
    opro[1, ] <- as.list(rep("pow", p))
    opro[2, ] <- as.list(rep(1,p))
    opro[3, ] <- as.list(-XC)
    opro[4, ] <- as.list(rep(0.5, p))
    opro[5, ] <- as.list(rep(0, p))
   }
else if(alpha == 0){ # Berg
    opro[1, ] <- as.list(rep("log", p))
    opro[2, ] <- as.list(rep(1,p))
    opro[3, ] <- as.list(-XC)
    opro[4, ] <- as.list(rep(0, p))
    opro[5, ] <- as.list(rep(0, p))
   }
else if(alpha == 2){ # Pearson
    opro[1, ] <- as.list(rep("pow", p))
    opro[2, ] <- as.list(rep(1,p))
    opro[3, ] <- as.list(XC)
    opro[4, ] <- as.list(rep(2, p))
    opro[5, ] <- as.list(rep(0, p))
   }
else if(alpha == 3){ # Silverman for Good
    opro[1, ] <- as.list(rep("pow", p))
    opro[2, ] <- as.list(rep(1,p))
    opro[3, ] <- as.list(XC)
    opro[4, ] <- as.list(rep(3, p))
    opro[5, ] <- as.list(rep(0, p))
   }
else
   stop('specified alpha not (yet) implemented')


P$scopt <- list(opro = opro)
P$dparam$intpnt_nl_tol_rel_gap <- rtol
if(length(control)){
    P$iparam <- control$iparam
    P$dparam <- control$dparam
    P$sparam <- control$sparam
}
z <- Rmosek::mosek(P, opts = list(verbose = verb))
status = z$sol$itr$solsta
if(status != "OPTIMAL") warning(paste("Solution status = ", status))
f <- z$sol$itr$xx[1:p]
g <- t(XL) %*% f
o <- order(x)
phi <- cbind(x[o],log(g[o]))
logLik <- sum(phi[,2])
z <- list(x = v, y = f, phi = phi, logLik = logLik, status = status)
class(z) <- "medde"
z
}
#' Plotting method for medde objects
#'
#' @param x object obtained from medde fitting
#' @param xlab label for horizontal axis
#' @param ylab label for vertical axis
#' @param ... other parameters to be passed to plot method
plot.medde <- function(x, xlab = "x", ylab = "f(x)", ...){
	plot.default(x, type = "l", xlab = xlab, ylab = ylab, ...)
	}


#' Quantile function for medde estimate
#'
#' Slightly modified version  borrowed from the package logcondens 
#' @param p vector of probabilities at which to evaluate the quantiles
#' @param medde fitted object from medde
#' @keywords nonparametric
#' @export
qmedde <- function (p, medde) { 
    if (any(p < 0 | p > 1)) stop("All p's must be in [0, 1]!\n")
    x <- medde$phi[,1]
    phi <- medde$phi[,2]

    Fhat <- function (x, phi) {
       Jexp <- function (x, y, v = 1) {
          m <- length(x)
          z <- exp(x)
          d <- y - x
          k <- (1:m)[abs(d) > 0.005]
          z[k] <- z[k] * (exp(v * d[k]) - 1)/d[k]
          k <- (1:m)[abs(d) <= 0.005]
          z[k] <- z[k] * 
              (v + d[k] * (v/2 + d[k] * (v/6 + d[k] * (v/24 + d[k] * v/120))))
          return(z)
          }
    n <- length(x)
    Fhat <- 1:n * 0
    dx <- diff(x)
    Fhat[2:n] <- cumsum(dx * Jexp(phi[1:(n - 1)], phi[2:n]))
    Fhat[n] <- max(Fhat[n], 1)
    Fhat
    }
    qloglin <- function (u, v) 
        ifelse(abs(v) > 1e-06, log(1 + ((exp(v) - 1) * u))/v, u + v * u * (1 - u)/2)

    Fhat <- Fhat(x,phi)
    n <- length(x)
    m <- length(p)
    qs <- rep(0,m)
    for (i in 1:m) {
        p0 <- p[i]
        if (p0 == 0) {
            q <- -Inf
        }
        if (p0 == 1) {
            q <- x[n]
        }
        if ((p0 > 0) && (p0 < 1)) {
            xj <- max(x[Fhat <= p0])
            j <- length(x[x <= xj])
            u <- (p0 - Fhat[j])/(Fhat[j + 1] - Fhat[j]) 
            v <- (x[j + 1] - x[j]) * (phi[j + 1] - phi[j])
            q <- xj + (x[j + 1] - x[j]) * qloglin(u,v)
        }
        qs[i] <- as.numeric(q)
    }
    qs
}
#' Random number generation from a medde estimate
#'
#' @param n number of observations desired in calls to rmedde
#' @param medde fitted medde object for calls in qmedde and rmedde
#' @param smooth option to draw random meddes from the smoothed density
#' @keywords nonparametric
#' @export
rmedde <- function(n, medde, smooth = TRUE) {
	z <- qmedde(sort(runif(n)), medde) 
	if(smooth){
	   x <- medde$x
	   fx <- medde$y
	   X <- medde$phi[,1]
	   dx <- diff(medde$x)[1]
	   v <- var(X) - sum(((x - mean(X))^2)*fx)*dx
	   z <- z + rnorm(n,sd = sqrt(v))
	   }
	z
}
