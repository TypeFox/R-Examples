#--------------------------------------------------------------------
#   svarmod.sb.R (npsp package)
#--------------------------------------------------------------------
#   kappasb(x, dk)
#   disc.sb(nx, dk, rmax)
#   fitsvar.sb.iso(esv, dk, nx, rmax, min.contrib, method, iter, tol)
#
#   (c) R. Fernandez-Casal         Last revision: Apr 2013
#--------------------------------------------------------------------
# PENDENTE:
#   - documentación
#   - @examples
#--------------------------------------------------------------------



#--------------------------------------------------------------------
#   kappasb(x, dk = 0)
#--------------------------------------------------------------------
#' Coefficients of an extended Shapiro-Botha variogram model 
#' 
#' Computes the coefficients of an extended Shapiro-Botha variogram model. 
#' @details If \code{dk >= 1}, the coefficients are computed as: 
#' \deqn{\kappa_d(x) = (2/x)^{(d-2)/2} \Gamma(d/2) J_{(d-2)/2}(x)} 
#' where \eqn{J_p} is the Bessel function of order \eqn{p}. 
#' \cr If \code{dk == 0}, the coefficients are computed as: 
#' \deqn{\kappa _\infty(x) = e^{-x^2}}
#' (corresponding to a model valid in any spatial dimension). 
#' \cr NOTE: some authors denote these functions as \eqn{\Omega_d}.
#' 
#' @param  x   numeric vector (on which the kappa function will be evaluated).
#' @param  dk  dimension of the kappa function.
#' @return
#' A vector with the coefficients of an extended Shapiro-Botha variogram model. 
#' @references
#' Shapiro, A. and Botha, J.D. (1991) Variogram fitting with a general class of 
#'   conditionally non-negative definite functions. \emph{Computational Statistics 
#'   and Data Analysis}, \bold{11}, 87-96. 
#' @seealso
#' \code{\link{svarmod.sb.iso}}, \code{\link{besselJ}}.
#' @examples
#' kappasb(seq(0, 6*pi, len = 10), 2)
#'   
#' curve(kappasb(x/5, 0), xlim = c(0, 6*pi), ylim = c(-1, 1), lty = 2)
#' for (i in 1:10) curve(kappasb(x, i), col = gray((i-1)/10), add = TRUE)
#' abline(h = 0, lty = 3)
#' @export
kappasb <- function(x, dk = 0) {
#--------------------------------------------------------------------
  if ( zeros <- any(index <- x < sqrt(.Machine$double.eps)) ) {
      dx <- dim(x)
      x <- x[!index]
      # Alternativamente se podría hacer pmax(x, .Machine$double.eps^0.5)
  }    
  res <- switch( min(dk + 1, 5),
      exp(-x*x),            # dk = 0
      cos(x),               # dk = 1
      besselJ(x, nu = 0),   # dk = 2
      sin(x)/x,             # dk = 3
      gamma(dk/2) * (2/x)^((dk-2)/2) * besselJ(x, nu = (dk-2)/2) # dk >=4    
    ) # switch
  if (zeros) { 
      res[!index] <- res
      res[index] <- 1
      dim(res) <- dx
  }     
  return(res)
} 



#--------------------------------------------------------------------
#   disc.sb(nx, dk = 0, rmax = 1)
#--------------------------------------------------------------------
#' Discretization nodes of a Shapiro-Botha variogram model
#' 
#' Computes the discretization nodes of a `nonparametric' extended Shapiro-Botha 
#' variogram model, following Gorsich and Genton (2004), as the scaled roots of 
#' Bessel functions.
#' @param  nx  number of discretization nodes.
#' @param  dk  dimension of the kappa function.
#' @param  rmax  maximum lag considered.
#' @details
#' If \code{dk >= 1}, the nodes are computed as: 
#' \deqn{x_i = q_i/rmax; i = 1,\ldots, nx,} where 
#' \eqn{q_i} are the first \eqn{n} roots of \eqn{J_{(d-2)/2}}, \eqn{J_p} 
#' is the Bessel function of order \eqn{p} and \eqn{rmax} 
#' is the maximum lag considered. The computation of the zeros of the Bessel  
#' function is done using the efficient algorithm developed by Ball (2000).
#' 
#' If \code{dk == 0} (corresponding to a model valid in any spatial dimension), 
#' the nodes are computed so the gaussian variogram models involved have
#' practical ranges: 
#    \deqn{r_i = ( 1 + 1.2(i-1))rmax/nx; i = 1,\ldots, nx.}
#'    \deqn{r_i = ( 1 + (i-1))rmax/nx; i = 1,\ldots, nx.}
#' @references
#' Ball, J.S. (2000) Automatic computation of zeros of Bessel functions and other
#'   special functions. \emph{SIAM Journal on Scientific Computing}, \bold{21}, 
#'   1458-1464.  
#'
#' Gorsich, D.J. and Genton, M.G. (2004) On the discretization of nonparametric 
#'   covariogram estimators. \emph{Statistics and Computing}, \bold{14}, 99-108.
#' @seealso
#' \code{\link{kappasb}}, \code{\link{fitsvar.sb.iso}}.
#' @examples
#' disc.sb( 12, 1, 1.0)
#' 
#' nx <- 1
#' dk <- 0
#' x <- disc.sb(nx, dk, 1.0)
#' h <- seq(0, 1, length = 100)
#' plot(h, kappasb(x * h, 0), type="l", ylim = c(0, 1))
#' abline(h = 0.05, lty = 2)
#' @export
disc.sb <- function(nx, dk = 0, rmax = 1) {   
#--------------------------------------------------------------------
    if (dk == 0) 
        # OJO: equiv. modelos gausianos, pueden aparecer inestabilidades
        # Nodos de discretización "geométricos"
        # return( 1.732/(seq(1, by = -1/nx, length = nx) * rmax) )
        # return( sqrt(3)/(seq(1/nx, by = 1.2/nx, length = nx) * rmax) )
        return( sqrt(3)/(seq(1/nx, 1, length = nx) * rmax) )
    # Let's go FORTRAN!
    #   subroutine disc_sbv(nx, x, dim, range)
    ret <-.Fortran( "disc_sbv", nx = as.integer(nx), x = double(nx), 
                  dk = as.integer(dk), as.double(rmax))
    return(ret$x)
}                  



#--------------------------------------------------------------------
#   fitsvar.sb.iso(esv, dk = ncol(esv$data$x), nx = NULL, rmax = esv$grid$max, 
#        min.contrib = 10, method = c("cressie", "equal", "npairs", "linear"), 
#        iter = 10, tol = sqrt(.Machine$double.eps)) {
#--------------------------------------------------------------------
#' Fit an isotropic Shapiro-Botha variogram model
#' 
#' Fits a `nonparametric' isotropic Shapiro-Botha variogram model by WLS through
#' quadratic programming.
#' Following Gorsich and Genton (2004), the nodes are selected as the scaled 
#' roots of Bessel functions (see \code{\link{disc.sb}}).
#' @param  esv pilot semivariogram estimate, a \code{\link{np.svar}}-\code{\link{class}} 
#'   (or \code{\link{svar.bin}}) object. Typically an output of the function
#'   \code{\link{np.svariso}}. 
#' @param  dk  dimension of the kappa function (\code{dk == 0} corresponds to a model 
#'   valid in any dimension; if \code{dk > 0}, it should be greater than 
#'   or equal to the dimension of the spatial process \code{ncol(esv$data$x)}).
#' @param  nx  number of discretization nodes. Defaults to \code{min(nesv - 1, 50)},
#' where \code{nesv} is the number of semivariogram estimates.
#' @param  rmax  maximum lag considered in the discretization
#'   (range of the fitted variogram on output).
#' @param  min.contrib  minimum number of (equivalent) contributing pairs 
#' (pilot estimates with a lower number are ignored, with a warning).
#' @param  method  string indicating the WLS fitting method to be used
#'   (e.g. \code{method = "cressie"}). See "Details" below.
#' @param  iter  maximum number of interations of the WLS algorithm (used only  
#'   if \code{method == "cressie"}).
#' @param  tol  absolute convergence tolerance (used only  
#'   if \code{method == "cressie"}).
#' @details
#' The fit is done using a (possibly iterated) weighted least squares criterion, minimizing: 
#' \deqn{WLS(\theta) = \sum_i w_i[(\hat{\gamma}(h_i)) -	\gamma(\theta; h_i)]^2.}
#' The different options for the argument \code{method} define the WLS algorithm used:
#' \describe{
#'  \item{\code{"cressie"}}{The default method. The procedure is
#'  iterative, with \eqn{w_i = 1} (OLS) used for the first step
#'  and with the weights recalculated at each iteration,
#'  following Cressie (1985), until convergence: \deqn{w_i =
#'  N(h_i)/\gamma(\hat{\theta}; h_i)^2,} where \eqn{N(h_i)}
#'  is the (equivalent) number of contributing pairs in the
#'  estimation at lag \eqn{h_i}.}
#'  \item{\code{"equal"}}{Ordinary least squares: \eqn{w_i = 1}.}
#'  \item{\code{"npairs"}}{\eqn{w_i = N(h_i).}} 
#'  \item{\code{"linear"}}{\eqn{w_i = N(h_i)/h_i^2} 
#'  (default fitting method in \pkg{gstat} package).} 
#' } 
#' Function \code{\link[quadprog]{solve.QP}} of \pkg{quadprog} package is used 
#' to solve the quadratic programming problem. If \code{nx} and/or \code{dim(esv)}
#' are large, this function may fail with error message "matrix D in quadratic 
#' function is not positive definite!". 
#' @return
#' Returns the fitted variogram model, an object of \code{\link{class}} \code{fitsvar} 
# (extending \code{\link{sb.iso}}: a \code{\link{svarmod}}) object
#' with an additional component \code{fit} containing:
#' \item{u}{vector of lags/distances.}
#' \item{sv}{vector of pilot semivariogram estimates.}
#' \item{fitted.sv}{vector of fitted semivariances.}
#' \item{wls}{value of the WLS objective function.}
#' \item{method}{string indicating the WLS fitting method used.}
#' \item{iter}{number of WLS iterations (if \code{method == "cressie"}).}   
#' 
#' @references
#' Ball, J.S. (2000) Automatic computation of zeros of Bessel functions and other
#'   special functions. \emph{SIAM Journal on Scientific Computing}, \bold{21}, 
#'   1458-1464.
#' 
#' Cressie, N. (1985) Fitting variogram models by weighted least squares.
#'   \emph{Mathematical Geology}, \bold{17}, 563-586. 
#' 
#' Cressie, N. (1993) \emph{Statistics for Spatial Data}. New York. Wiley.
#' 
#' Fernandez Casal R., Gonzalez Manteiga W. and  Febrero Bande M. (2003) 
#' Flexible Spatio-Temporal Stationary Variogram Models, 
#' \emph{Statistics and Computing}, \bold{13}, 127-136.
#'
#' Gorsich, D.J. and Genton, M.G. (2004) On the discretization of nonparametric 
#'   covariogram estimators. \emph{Statistics and Computing}, \bold{14}, 99-108.
#' 
#' Shapiro, A. and Botha, J.D. (1991) Variogram fitting with a general class of 
#'   conditionally non-negative definite functions. \emph{Computational Statistics 
#'   and Data Analysis}, \bold{11}, 87-96. 
#' @seealso
#' \code{\link{svarmod.sb.iso}}, \code{\link{disc.sb}}, \code{\link{plot.fitsvar}}.
#' @export
#--------------------------------------------------------------------
fitsvar.sb.iso <- function(esv, dk = ncol(esv$data$x), nx = NULL, rmax = esv$grid$max, 
        min.contrib = 10, method = c("cressie", "equal", "npairs", "linear"), 
        iter = 10, tol = sqrt(.Machine$double.eps)) {
#   PENDENTE:
#     - Rounding errors?  w <- w/sum(w)
#     - rematar documentación: details, examples, ...
#     - Versión preliminar, final 'fit.svar.sb' válida para modelos anisotrópicos
#     - verificar missing values
#     - nodes un vector de ptos de discretización
#--------------------------------------------------------------------
    if (!inherits(esv, "svar.bin"))
      stop("function only works for objects of class (or extending) 'svar.bin'.")
    # if (esv$svar$type != "isotropic")
    if (esv$grid$nd != 1)
      stop("pilot variogram estimates 'esv' must be isotropic.")
    method <- match.arg(method)    
    if (method != "cressie") iter <- 1         
    # if (!require(quadprog)) stop("'quadprog' package is required.")
    # Let's go...
    u <- as.numeric(coords(esv))
    if (inherits(esv, "np.svar")) { 
        # "np.svar" class 
        v <- esv$est
        if (!is.null(esv$locpol$hat)) {
            # Aproximación varianza estilo Cressie (suponiendo independencia) para ajuste wls
            # PENDIENTE: ESCRIBIR/REVISAR ESTAS CUENTAS            
            n <- 1 / with(esv, rowSums(locpol$hat^2 / 
                pmax(matrix(binw, nrow=grid$n, ncol=grid$n, byrow=TRUE), 1))) # nº equivalente de aportaciones
        } else {
            n <- esv$binw # nº de aportaciones            
        }    
    } else {
        # "svar.bin" class 
        v <- esv$biny
        n <- esv$binw # nº de aportaciones    
    }
    if (!all(index <- n >= min.contrib)) {
        warning("some pilot semivariogram estimates will be ignored (contrib < min.contrib)")
        u <- u[index]
        v <- v[index]
        n <- n[index]
    }
    n.esv <- length(u)
    # Discretization points
    if (is.null(nx)) {
        nx <- min(n.esv - 1, 50)
    } else {
        if (nx >= length(u)) {
          warning("'nx' must be less than the number of variogram estimates (corrected)")
          nx <- n.esv - 1
        }  
    }  
    x <- disc.sb( nx, dk, rmax)
    # M(1:nesv,1:npar)
    # M <- cbind(-outer(u, x, function(u, x) kappasb(u*x, dk)), 1)      
    M <- cbind( -kappasb(outer(u, x), dk), 1)
    # Reescalar pesos
    w <- switch(method,
        npairs =  n,
        linear =   n/pmax(u^2, .Machine$double.eps^0.5),
                  1 # default: "cressie", "equal"
        )
    w <- w/sum(w)
    # (iterative) quadratic programming
    n.par <- nx + 1 
    bvec <- rep(0, n.par)
    Amat <- diag(n.par)
    Amat[1:nx, n.par] <- -1   
    # wls loop
    i <- 0
    conver <- FALSE
    while( (i < iter) & !conver) {
        i <- i + 1
        # d = t(M) %*% diag(w) %*% v
        dvec <- drop((w * v) %*% M)
        # D = t(M) %*% diag(w) %*% M
        Dmat <- crossprod(sqrt(w) * M)
        # solve min(-d^T b + 1/2 b^T D b) with the constraints A^T·b >= b_0
        res <- solve.QP(Dmat, dvec, Amat, bvec)
        fit <- drop(M %*% res$solution)
        if (i > 1) {
            # Absolute parameter difference convergence criteria
            conver <-  max(abs(sol - res$solution)) < tol
            w <- n/pmax(fit^2, .Machine$double.eps^0.5)
            w <- w/sum(w)
        }
        wls <- sum(w*(v-fit)^2)
        sol <- res$solution
    }    
    if (method != "cressie") conver <- TRUE else iter <- i
    if (!conver) warning("the wls algorithm did not converge. \n",
          "   Maximum number of iterations exceeded; \n",
          "   the current values may be an approximate solution, \n",
          "   or 'tol' is too big or 'iter' is too small.")
    result <- svarmod.sb.iso( dk = dk, x = x, z = sol[-n.par], nu = sol[n.par], 
          range = rmax)
    result$fit <- list(u = u, sv = v, fitted.sv = fit, wls = wls, 
        method = method, iter = iter)  
    oldClass(result) <- c("fitsvar", oldClass(result))     
    return(result)
#--------------------------------------------------------------------
} # fitsvar.sb.iso




