##' Fit a linear model with elastic-net regularization
##'
##' Adjust a linear model with elastic-net regularization, mixing a
##' (possibly weighted) \eqn{\ell_1}{l1}-norm (LASSO) and a
##' (possibly structured) \eqn{\ell_2}{l2}-norm (ridge-like). The
##' solution path is computed at a grid of values for the
##' \eqn{\ell_1}{l1}-penalty, fixing the amount of \eqn{\ell_2}{l2}
##' regularization. See details for the criterion optimized.
##'
##' @param x matrix of features, possibly sparsely encoded
##' (experimental). Do NOT include intercept. When normalized os
##' \code{TRUE}, coefficients will then be rescaled to the original
##' scale.
##'
##' @param y response vector.
##'
##' @param lambda1 sequence of decreasing \eqn{\ell_1}{l1}-penalty
##' levels. If \code{NULL} (the default), a vector is generated with
##' \code{nlambda1} entries, starting from a guessed level
##' \code{lambda1.max} where only the intercept is included, then
##' shrunken to \code{min.ratio*lambda1.max}.
##'
##' @param lambda2 real scalar; tunes the \eqn{\ell_2}{l2} penalty in
##' the Elastic-net. Default is 0.01. Set to 0 to recover the Lasso.
##'
##' @param penscale vector with real positive values that weight the
##' \eqn{\ell_1}{l1}-penalty of each feature. Default set all weights
##' to 1.
##'
##' @param struct matrix structuring the coefficients (preferably
##' sparse). Must be at least positive semidefinite (this is checked
##' internally if the \code{checkarg} argument is \code{TRUE}). The
##' default uses the identity matrix. See details below.
##'
##' @param intercept logical; indicates if an intercept should be
##' included in the model. Default is \code{TRUE}.
##'
##' @param normalize logical; indicates if variables should be
##' normalized to have unit L2 norm before fitting.  Default is
##' \code{TRUE}.
##'
##' @param naive logical; Compute either 'naive' of classic
##' elastic-net as defined in Zou and Hastie (2006): the vector of
##' parameters is rescaled by a coefficient \code{(1+lambda2)} when
##' \code{naive} equals \code{FALSE}.  No rescaling otherwise.
##' Default is \code{FALSE}.
##'
##' @param nlambda1 integer that indicates the number of values to put
##' in the \code{lambda1} vector.  Ignored if \code{lambda1} is
##' provided.
##'
##' @param min.ratio minimal value of \eqn{\ell_1}{l1}-part of the
##' penalty that will be tried, as a fraction of the maximal
##' \code{lambda1} value. A too small value might lead to unstability
##' at the end of the solution path corresponding to small
##' \code{lambda1} combined with \eqn{\lambda_2=0}{lambda2=0}.  The
##' default value tries to avoid this, adapting to the
##' '\eqn{n<p}{n<p}' context. Ignored if \code{lambda1} is provided.
##'
##' @param max.feat integer; limits the number of features ever to
##' enter the model; i.e., non-zero coefficients for the Elastic-net:
##' the algorithm stops if this number is exceeded and \code{lambda1}
##' is cut at the corresponding level. Default is
##' \code{min(nrow(x),ncol(x))} for small \code{lambda2} (<0.01) and
##' \code{min(4*nrow(x),ncol(x))} otherwise. Use with care, as it
##' considerably changes the computation time.
##'
##' @param beta0 a starting point for the vector of parameter. When
##' \code{NULL} (the default), will be initialized at zero. May save
##' time in some situation.
##'
##' @param control list of argument controlling low level options of
##' the algorithm --use with care and at your own risk-- :
##' \itemize{%
##'
##' \item{\code{verbose}: }{integer; activate verbose mode --this one
##' is not too much risky!-- set to \code{0} for no output; \code{1}
##' for warnings only, and \code{2} for tracing the whole
##' progression. Default is \code{1}. Automatically set to \code{0}
##' when the method is embedded within cross-validation or stability
##' selection.}
##'
##' \item{\code{timer}: }{logical; use to record the timing of the
##' algorithm. Default is \code{FALSE}.}
##'
##' \item{\code{max.iter}: }{the maximal number of iteration used to
##' solve the problem for a given value of lambda1. Default is 500.}
##'
##' \item{\code{method}: }{a string for the underlying solver
##' used. Either \code{"quadra"}, \code{"pathwise"} or
##' \code{"fista"}. Default is \code{"quadra"}.}
##'
##' \item{\code{threshold}: }{a threshold for convergence. The
##' algorithm stops when the optimality conditions are fulfill up to
##' this threshold. Default is \code{1e-7} for \code{"quadra"} and
##' \code{1e-2} for the first order methods.}
##'
##' \item{\code{monitor}: }{indicates if a monitoring of the
##' convergence should be recorded, by computing a lower bound between
##' the current solution and the optimum: when \code{'0'} (the
##' default), no monitoring is provided; when \code{'1'}, the bound
##' derived in Grandvalet et al. is computed; when \code{'>1'}, the
##' Fenchel duality gap is computed along the algorithm.}
##' }
##'
##' @param checkargs logical; should arguments be checked to
##' (hopefully) avoid internal crashes? Default is
##' \code{TRUE}. Automatically set to \code{FALSE} when calls are made
##' from cross-validation or stability selection procedures.
##'
##' @return an object with class \code{quadrupen}, see the
##' documentation page \code{\linkS4class{quadrupen}} for details.
##'
##' @note The optimized criterion is the following: \if{latex}{\deqn{%
##' \hat{\beta}_{\lambda_1,\lambda_1} = \arg \min_{\beta} \frac{1}{2}
##' (y - X \beta)^T (y - X \beta) + \lambda_1 \|D \beta \|_{1} +
##' \frac{\lambda_2}{2} \beta^T S \beta, }} \if{html}{\out{ <center>
##' &beta;<sup>hat</sup>
##' <sub>&lambda;<sub>1</sub>,&lambda;<sub>2</sub></sub> =
##' argmin<sub>&beta;</sub> 1/2 RSS(&beta) + &lambda;<sub>1</sub>
##' &#124; D &beta; &#124;<sub>1</sub> + &lambda;/2 <sub>2</sub>
##' &beta;<sup>T</sup> S &beta;, </center> }}
##' \if{text}{\deqn{beta.hat(lambda1, lambda2) = argmin_beta 1/2
##' RSS(beta) + lambda1 |D beta|1 + lambda2 beta' S beta,}} where
##' \eqn{D}{D} is a diagonal matrix, whose diagonal terms are provided
##' as a vector by the \code{penscale} argument. The \eqn{\ell_2}{l2}
##' structuring matrix \eqn{S}{S} is provided via the \code{struct}
##' argument, a positive semidefinite matrix (possibly of class
##' \code{Matrix}).
##'
##' @seealso See also \code{\linkS4class{quadrupen}},
##' \code{\link{plot,quadrupen-method}} and \code{\link{crossval}}.
##' @name elastic.net
##' @rdname elastic.net
##' @keywords models, regression
##'
##' @examples
##' ## Simulating multivariate Gaussian with blockwise correlation
##' ## and piecewise constant vector of parameters
##' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
##' cor <- 0.75
##' Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
##' Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
##' Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo)
##' diag(Sigma) <- 1
##' n <- 50
##' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
##' y <- 10 + x %*% beta + rnorm(n,0,10)
##'
##' labels <- rep("irrelevant", length(beta))
##' labels[beta != 0] <- "relevant"
##' ## Comparing the solution path of the LASSO and the Elastic-net
##' plot(elastic.net(x,y,lambda2=0), label=labels) ## a mess
##' plot(elastic.net(x,y,lambda2=10), label=labels) ## a lot better
##'
##' @export
elastic.net <- function(x,
                        y,
                        lambda1   = NULL,
                        lambda2   = 0.01,
                        penscale  = rep(1,p),
                        struct    = NULL,
                        intercept = TRUE,
                        normalize = TRUE,
                        naive     = FALSE,
                        nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
                        min.ratio = ifelse(n<=p,1e-2,1e-4),
                        max.feat  = ifelse(lambda2<1e-2,min(n,p),min(4*n,p)),
                        beta0     = NULL,
                        control   = list(),
                        checkargs = TRUE) {

  ## ===================================================
  ## CHECKS TO (PARTIALLY) AVOID CRASHES OF THE C++ CODE
  p <- ncol(x) # problem size
  n <- nrow(x) # sample size
  if (checkargs) {
    if (is.data.frame(x))
      x <- as.matrix(x)
    if(!inherits(x, c("matrix", "dgCMatrix")))
      stop("x has to be of class 'matrix' or 'dgCMatrix'.")
    if(any(is.na(x)))
      stop("NA value in x not allowed.")
    if(!is.numeric(y))
      stop("y has to be of type 'numeric'")
    if(n != length(y))
      stop("x and y have not correct dimensions")
    if(length(penscale) != p)
      stop("penscale must have ncol(x) entries")
    if (any(penscale <= 0))
      stop("weights in penscale must be positive")
    if(!inherits(lambda2, "numeric") | length(lambda2) > 1)
      stop("lambda2 must be a scalar.")
    if(lambda2 < 0)
      stop("lambda2 must be a non negative scalar.")
    if (!is.null(lambda1)) {
      if(any(lambda1 <= 0))
        stop("entries inlambda1 must all be postive.")
      if(is.unsorted(rev(lambda1)))
        stop("lambda1 values must be sorted in decreasing order.")
      if (length(lambda1)>1 & !is.null(beta0))
        warning("providing beta0 for a serie of l1 penalties mught be inefficient.")
    }
    if(min.ratio < 0)
        stop("min.ratio must be non negative.")
    if (!is.null(struct)) {
      if (ncol(struct) != p | ncol(struct) != p)
          stop("struct must be a (square) positive semidefinite matrix.")
      if (any(eigen(struct,only.values=TRUE)$values<0))
          stop("struct must be a (square) positive semidefinite matrix.")
      if(!inherits(struct, "dgCMatrix"))
          struct <- as(struct, "dgCMatrix")
    }      
    if (!is.null(beta0)) {
      beta0 <- as.numeric(beta0)
      if (length(beta0) != p)
        stop("beta0 must be a vector with p entries.")
    }
    if (length(max.feat)>1)
        stop("max.feat must be an integer.")
    if(is.numeric(max.feat) & !is.integer(max.feat))
        max.feat <- as.integer(max.feat)
  }

  return(quadrupen(beta0 = beta0,
                   x=x,
                   y=y,
                   penalty   = "elastic.net",
                   lambda1   = lambda1,
                   lambda2   = lambda2,
                   penscale  = penscale,
                   struct    = struct,
                   intercept = intercept,
                   normalize = normalize,
                   naive     = naive,
                   nlambda1  = nlambda1,
                   min.ratio = min.ratio,
                   max.feat  = max.feat,
                   control   = control)
         )
}

##' Fit a linear model with infinity-norm plus ridge-like regularization
##'
##' Adjust a linear model penalized by a mixture of a (possibly
##' weighted) \eqn{\ell_\infty}{l-infinity}-norm (bounding the
##' magnitude of the parameters) and a (possibly structured)
##' \eqn{\ell_2}{l2}-norm (ridge-like). The solution path is computed
##' at a grid of values for the infinity-penalty, fixing the amount of
##' \eqn{\ell_2}{l2} regularization. See details for the criterion
##' optimized.
##'
##' @param x matrix of features, possibly sparsely encoded
##' (experimental). Do NOT include intercept. When normalized os
##' \code{TRUE}, coefficients will then be rescaled to the original
##' scale.
##'
##' @param y response vector.
##'
##' @param lambda1 sequence of decreasing \eqn{\ell_\infty}{l-infinity}
##' penalty levels. If \code{NULL} (the default), a vector is
##' generated with \code{nlambda1} entries, starting from a guessed
##' level \code{lambda1.max} where only the intercept is included,
##' then shrunken to \code{min.ratio*lambda1.max}.
##'
##' @param lambda2 real scalar; tunes the \eqn{\ell_2}{l2}-penalty in
##' the bounded regression. Default is 0.01. Set to 0 to regularize
##' only by the infinity norm (be careful regarding numerical
##' stability in that case, particularly in the high dimensional
##' setting).
##'
##' @param penscale vector with real positive values that weight the
##' infinity norm of each feature. Default set all weights to 1. See
##' details below.
##'
##' @param struct matrix structuring the coefficients.  Must be at
##' least positive semidefinite (this is checked internally if the
##' \code{checkarg} argument is \code{TRUE}). The
##' default uses the identity matrix. See details below.
##'
##' @param intercept logical; indicates if an intercept should be
##' included in the model. Default is \code{TRUE}.
##'
##' @param normalize logical; indicates if variables should be
##' normalized to have unit L2 norm before fitting.  Default is
##' \code{TRUE}.
##'
##' @param naive logical; Compute either 'naive' of 'classic' bounded
##' regression: mimicing the Elastic-net, the vector of parameters is
##' rescaled by a coefficient \code{(1+lambda2)} when \code{naive}
##' equals \code{FALSE}.  No rescaling otherwise. Default is
##' \code{FALSE}.
##'
##' @param nlambda1 integer that indicates the number of values to put
##' in the \code{lambda1} vector.  Ignored if \code{lambda1} is
##' provided.
##'
##' @param min.ratio minimal value of infinity-part of the penalty
##' that will be tried, as a fraction of the maximal \code{lambda1}
##' value. A too small value might lead to unstability at the end of
##' the solution path corresponding to small \code{lambda1}.  The
##' default value tries to avoid this, adapting to the
##' '\eqn{n<p}{n<p}' context. Ignored if \code{lambda1} is provided.
##'
##' @param max.feat integer; limits the number of features ever to
##' enter the model: in our implementation of the bounded regression,
##' it corresponds to the variables which have left the boundary along
##' the path.  The algorithm stops if this number is exceeded and
##' \code{lambda1} is cut at the corresponding level. Default is
##' \code{min(nrow(x),ncol(x))} for small \code{lambda2} (<0.01) and
##' \code{min(4*nrow(x),ncol(x))} otherwise. Use with care, as it
##' considerably changes the computation time.
##'
##' @param control list of argument controlling low level options of
##' the algorithm --use with care and at your own risk-- :
##' \itemize{%
##'
##' \item{\code{verbose}: }{integer; activate verbose mode --this one
##' is not too much risky!-- set to \code{0} for no output; \code{1}
##' for warnings only, and \code{2} for tracing the whole
##' progression. Default is \code{1}. Automatically set to \code{0}
##' when the method is embedded within cross-validation or stability
##' selection.}
##'
##' \item{\code{timer}: }{logical; use to record the timing of the
##' algorithm. Default is \code{FALSE}.}
##'
##' \item{\code{max.iter}: }{the maximal number of iteration used to
##' solve the problem for a given value of \code{lambda1}. Default is
##' 500.}
##'
##' \item{\code{method}: }{a string for the underlying solver
##' used. Either \code{"quadra"} or \code{"fista"} are available for
##' bounded regression. Default is \code{"quadra"}.}
##'
##' \item{\code{threshold}: }{a threshold for convergence. The
##' algorithm stops when the optimality conditions are fulfill up to
##' this threshold. Default is \code{1e-7} for \code{"quadra"} and
##' \code{1e-2} for \code{"fista"}.}
##'
##' \item{\code{bulletproof}: }{logical; indicates if the bulletproof
##' mode should be used while running the \code{"quadra"}
##' method. Default is \code{TRUE}. See details below.}}
##'
##' @param checkargs logical; should arguments be checked to
##' (hopefully) avoid internal crashes? Default is
##' \code{TRUE}. Automatically set to \code{FALSE} when calls are made
##' from cross-validation or stability selection procedures.
##'
##' @return an object with class \code{quadrupen}, see the
##' documentation page \code{\linkS4class{quadrupen}} for details.
##'
##' @note The optimized criterion is  \if{latex}{\deqn{%
##' \hat{\beta}_{\lambda_1,\lambda_1} = \arg \min_{\beta} \frac{1}{2}
##' (y - X \beta)^T (y - X \beta) + \lambda_1 \|D \beta \|_{\infty} +
##' \frac{\lambda_2}{2} \beta^T S \beta, }} \if{html}{\out{ <center>
##' &beta;<sup>hat</sup>
##' <sub>&lambda;<sub>1</sub>,&lambda;<sub>2</sub></sub> =
##' argmin<sub>&beta;</sub> 1/2 RSS(&beta) + &lambda;<sub>1</sub>
##' &#124; D &beta; &#124;<sub>&infin;</sub> + &lambda;/2 <sub>2</sub>
##' &beta;<sup>T</sup> S &beta;, </center> }}
##' \if{text}{\deqn{beta.hat(lambda1, lambda2) = argmin_beta 1/2
##' RSS(beta) + lambda1 max|D beta| + lambda2 beta' S beta,}} where
##' \eqn{D}{D} is a diagonal matrix, whose diagonal terms are provided
##' as a vector by the \code{penscale} argument. The \eqn{\ell_2}{l2}
##' structuring matrix \eqn{S}{S} is provided via the \code{struct}
##' argument, a positive semidefinite matrix (possibly of class
##' \code{Matrix}).
##'
##' Note that the quadratic algorithm for the bounded regression may
##' become unstable along the path because of singularity of the
##' underlying problem, e.g. when there are too much correlation or
##' when the size of the problem is close to or smaller than the
##' sample size. In such cases, it might be a good idea to switch to
##' the proximal solver, slower yet more robust. This is the strategy
##' adopted by the \code{'bulletproof'} mode, that will send a warning
##' while switching the method to \code{'fista'} and keep on
##' optimizing on the remainder of the path. When \code{bulletproof}
##' is set to \code{FALSE}, the algorithm stops at an early stage of
##' the path of solutions. Hence, users should be careful when
##' manipulating the resulting \code{'quadrupen'} object, as it will
##' not have the size expected regarding the dimension of the
##' \code{lambda1} argument.
##'
##' Singularity of the system can also be avoided with a larger
##' \eqn{\ell_2}{l2}-regularization, via \code{lambda2}, or a
##' "not-too-small" \eqn{\ell_\infty}{l-infinity} regularization, via
##' a larger \code{'min.ratio'} argument.
##'
##' @seealso See also \code{\linkS4class{quadrupen}},
##' \code{\link{plot,quadrupen-method}} and \code{\link{crossval}}.
##' @name bounded.reg
##' @rdname bounded.reg
##' @keywords models, regression
##'
##' @examples
##' ## Simulating multivariate Gaussian with blockwise correlation
##' ## and piecewise constant vector of parameters
##' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
##' cor <- 0.75
##' Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
##' Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
##' Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo)
##' diag(Sigma) <- 1
##' n <- 50
##' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
##' y <- 10 + x %*% beta + rnorm(n,0,10)
##'
##' ## Infinity norm without/with an additional l2 regularization term
##' ## and with structuring prior
##' labels <- rep("irrelevant", length(beta))
##' labels[beta != 0] <- "relevant"
##' plot(bounded.reg(x,y,lambda2=0) , label=labels) ## a mess
##' plot(bounded.reg(x,y,lambda2=10), label=labels) ## good guys are at the boundaries
##'
##' @export
bounded.reg <- function(x,
                        y,
                        lambda1   = NULL,
                        lambda2   = 0.01,
                        penscale  = rep(1,p),
                        struct    = NULL,
                        intercept = TRUE,
                        normalize = TRUE,
                        naive     = FALSE,
                        nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
                        min.ratio = ifelse(n<=p,1e-2,1e-3),
                        max.feat  = ifelse(lambda2<1e-2,min(n,p),min(4*n,p)),
                        control   = list(),
                        checkargs = TRUE) {

  ## ===================================================
  ## CHECKS TO (PARTIALLY) AVOID CRASHES OF THE C++ CODE
  p <- ncol(x) # problem size
  n <- nrow(x) # sample size
  if (checkargs) {
    if (is.data.frame(x))
      x <- as.matrix(x)
    if(!inherits(x, c("matrix", "dgCMatrix")))
      stop("x has to be of class 'matrix' or 'dgCMatrix'.")
    if(any(is.na(x)))
      stop("NA value in x not allowed.")
    if(!is.numeric(y))
      stop("y has to be of type 'numeric'")
    if(n != length(y))
      stop("x and y have not correct dimensions")
    if(length(penscale) != p)
      stop("penscale must have ncol(x) entries")
    if (any(penscale <= 0))
      stop("weights in penscale must be positive")
    if(!inherits(lambda2, "numeric") | length(lambda2) > 1)
      stop("lambda2 must be a scalar.")
    if(lambda2 < 0)
      stop("lambda2 must be a non negative scalar.")
    if (!is.null(lambda1)) {
      if(any(lambda1 <= 0))
        stop("entries inlambda1 must all be postive.")
      if(is.unsorted(rev(lambda1)))
        stop("lambda1 values must be sorted in decreasing order.")
    }
    if(min.ratio < 0)
      stop("min.ratio must be non negative.")
    if (!is.null(struct)) {
      if (ncol(struct) != p | ncol(struct) != p)
        stop("struct must be a (square) positive definite matrix.")
      if (any(eigen(struct,only.values=TRUE)$values<=.Machine$double.eps))
        stop("struct must be a (square) positive definite matrix.")
      if(!inherits(struct, "dgCMatrix"))
        struct <- as(struct, "dgCMatrix")
    }
    
    if (length(max.feat)>1)
      stop("max.feat must be an integer.")
    if(is.numeric(max.feat) & !is.integer(max.feat))
      max.feat <- as.integer(max.feat)
  }

  return(quadrupen(beta0 = NULL,
                   x=x,
                   y=y,
                   penalty   = "bounded.reg",
                   lambda1   = lambda1,
                   lambda2   = lambda2,
                   penscale  = penscale,
                   struct    = struct,
                   intercept = intercept,
                   normalize = normalize,
                   naive     = naive,
                   nlambda1  = nlambda1,
                   min.ratio = min.ratio,
                   max.feat  = max.feat,
                   control   = control)
         )
}

quadrupen <- function(beta0    ,
                      x        ,
                      y        ,
                      penalty  ,
                      lambda1  ,
                      lambda2  ,
                      penscale ,
                      struct   ,
                      intercept,
                      normalize,
                      naive    ,
                      nlambda1 ,
                      min.ratio,
                      max.feat ,
                      control) {

  n <- nrow(x)
  p <- ncol(x)
  ## ============================================
  ## RECOVERING LOW LEVEL OPTIONS
  quadra <- TRUE
  if (!is.null(control$method)) {
    if (control$method != "quadra") {
      quadra <- FALSE
    }
  }
  ctrl <- list(verbose      = 1, # default control options
               timer        =  FALSE,
               max.iter     = max(500,p),
               method       = "quadra",
               threshold    = ifelse(quadra, 1e-7, 1e-2),
               monitor      = 0,
               bulletproof  = TRUE,
               usechol      = TRUE)
  ctrl[names(control)] <- control # overwritten by user specifications
  if (ctrl$timer) {r.start <- proc.time()}

  ## ======================================================
  ## STARTING C++ CALL TO ENET_LS
  if (ctrl$timer) {cpp.start <- proc.time()}
  if (penalty == "elastic.net") {
    out <- .Call("elastic_net",
                 beta0        ,
                 x            ,
                 y            ,
                 struct       ,
                 lambda1      ,
                 nlambda1     ,
                 min.ratio    ,
                 penscale     ,
                 lambda2      ,
                 intercept    ,
                 normalize    ,
                 rep(1,n)     ,
                 naive        ,
                 ctrl$thresh  ,
                 ctrl$max.iter,
                 max.feat     ,
                 switch(ctrl$method,
                        quadra   = 0,
                        pathwise = 1,
                        fista    = 2, 0),
                 ctrl$verbose,
                 inherits(x, "sparseMatrix"),
                 ctrl$usechol,
                 ctrl$monitor,
                 package = "quadrupen")
    coefficients <- sparseMatrix(i = out$iA+1,
                                 j = out$jA+1,
                                 x = c(out$nzeros),
                                 dims=c(length(out$lambda1),p))
    active.set   <- sparseMatrix(i = out$iA+1,
                                 j = out$jA+1,
                                 dims=c(length(out$lambda1),p))
  }
  if (penalty == "bounded.reg") {
    out <- .Call("bounded_reg",
                 x,
                 y,
                 struct,
                 lambda1      ,
                 nlambda1     ,
                 min.ratio    ,
                 penscale     ,
                 lambda2      ,
                 intercept    ,
                 normalize    ,
                 rep(1,n)     ,
                 naive        ,
                 ctrl$thresh  ,
                 ctrl$max.iter,
                 max.feat     ,
                 ifelse(ctrl$method=="fista",1,0),
                 ctrl$verbose,
                 inherits(x, "sparseMatrix"),
                 ctrl$bulletproof,
                 package = "quadrupen")
    coefficients <- Matrix(out$coefficients)
    active.set   <- sparseMatrix(i = out$iB+1,
                                 j = out$jB+1,
                                 dims=c(length(out$lambda1),p))
    out$delta.hat  <- NULL ## not implemented
    out$delta.star <- NULL ##
  }
  ## END OF CALL
  if (ctrl$timer) {
    internal.timer <- (proc.time() - cpp.start)[3]
    external.timer <- (proc.time() - r.start)[3]
  } else {
    internal.timer <- NULL
    external.timer <- NULL
  }

  ## ======================================================
  ## BUILDING THE QUADRUPEN OBJECT
  out$converge[out$converge == 0] <- "converged"
  out$converge[out$converge == 1] <- "max # of iterate reached"
  out$converge[out$converge == 2] <- "max # of feature reached"
  out$converge[out$converge == 3] <- "system has become singular"
  monitoring  <- list()
  monitoring$it.active <- c(out$it.active)
  monitoring$it.optim  <- c(out$it.optim)
  monitoring$max.grad  <- c(out$max.grd)
  monitoring$status    <- c(out$converge)
  monitoring$pensteps.timer <- c(out$timing)
  monitoring$external.timer <- external.timer
  monitoring$internal.timer <- internal.timer
  monitoring$dist.to.opt    <- c(out$delta.hat)
  monitoring$dist.to.str    <- c(out$delta.star)
  dim.names <- list()
  dimnames(coefficients)[[1]] <- round(c(out$lambda1),3)
  dimnames(coefficients)[[2]] <- 1:p
  mu <- drop(out$mu)

  ## FITTED VALUES AND RESIDUALS...
  if (intercept) {
    fitted <- sweep(tcrossprod(x,coefficients),2L,-mu,check.margin=FALSE)
  } else {
    mu <- 0
    fitted <- tcrossprod(x,coefficients)
  }
  residuals <- apply(fitted,2,function(y.hat) y-y.hat)
  normy     <- ifelse(intercept,sum((y-mean(y))^2),sum(y^2))
  r.squared <- 1-colSums(residuals^2)/normy

  return(new("quadrupen",
             coefficients = coefficients   ,
             active.set   = active.set     ,
             intercept    = intercept      ,
             mu           = mu             ,
             meanx        = drop(out$meanx),
             normx        = drop(out$normx),
             fitted       = fitted         ,
             residuals    = residuals      ,
             r.squared    = r.squared      ,
             penscale     = penscale       ,
             penalty      = penalty        ,
             naive        = naive          ,
             lambda1      = c(out$lambda1) ,
             lambda2      = lambda2        ,
             struct       = struct         ,
             monitoring   = monitoring     ,
             control      = ctrl))

}

standardize <- function(x,y,intercept,normalize,penscale,zero=.Machine$double.eps) {

  n <- length(y)
  p <- ncol(x)
  ## ============================================
  ## INTERCEPT AND NORMALIZATION TREATMENT
  if (intercept) {
    xbar <- colMeans(x)
    ybar <- mean(y)
  } else {
    xbar <- rep(0,p)
    ybar <- 0
  }

  ## ============================================
  ## NORMALIZATION
  if (normalize) {
    normx <- sqrt(drop(colSums(x^2)- n*xbar^2))
    if (any(normx < zero)) {
      warning("A predictor has no signal: you should remove it.")
      normx[abs(normx) < zero] <- 1 ## dirty way to handle 0/0
    }
    ## normalizing the predictors...
    x <- sweep(x, 2L, normx, "/", check.margin = FALSE)
    ## xbar is scaled to handle internaly the centering of X for
    ## sparsity purpose
    xbar <- xbar/normx
  } else {
    normx <- rep(1,p)
  }
  normy <- sqrt(sum(y^2))

  ## and now normalize predictors according to penscale value
  if (any(penscale != 1)) {
    x <- sweep(x, 2L, penscale, "/", check.margin=FALSE)
    xbar <- xbar/penscale
  }
  ## Computing marginal correlation
  if (intercept) {
    xty   <- drop(crossprod(y-ybar,sweep(x,2L,xbar)))
  } else {
    xty   <- drop(crossprod(y,x))
  }

  return(list(xbar=xbar, ybar=ybar, normx=normx, normy=normy, xty=xty, x=x))
}

## ======================================================
## GENERATE A GRID OF PENALTY IF NONE HAS BEEN PROVIDED
get.lambda1.l1 <- function(xty,nlambda1,min.ratio) {
  lmax <- max(abs(xty))
  return(10^seq(log10(lmax), log10(min.ratio*lmax), len=nlambda1))
}

get.lambda1.li <- function(xty,nlambda1,min.ratio) {
  lmax <- sum(abs(xty))
  return(10^seq(log10(lmax), log10(min.ratio*lmax), len=nlambda1))
}
