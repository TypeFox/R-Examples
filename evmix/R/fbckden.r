#' @export
#' 
#' @title Cross-validation MLE Fitting of Boundary Corrected Kernel Density Estimation
#'   Using a Variety of Approaches
#'
#' @description Maximum likelihood estimation for fitting boundary corrected 
#' kernel density estimator using a variety of approaches (and many possible kernels),
#' by treating it as a mixture model.
#'
#' @inheritParams fkden
#' @inheritParams fkdengpd
#' @inheritParams bckden
#' @inheritParams fgpd
#' 
#' @details The boundary corrected kernel density estimator using a variety of approaches
#' (and many possible kernels) is fitted to the entire dataset using
#' cross-validation maximum likelihood estimation. The estimated bandwidth,
#' variance and standard error are automatically output. 
#' 
#' The log-likelihood and negative log-likelihood are also provided for wider
#' usage, e.g. constructing your own extreme value
#' mixture models or profile likelihood functions. The parameter
#' \code{lambda} must be specified in the negative log-likelihood
#' \code{\link[evmix:fbckden]{nlbckden}}.
#' 
#' Log-likelihood calculations are carried out in
#' \code{\link[evmix:fbckden]{lbckden}}, which takes bandwidths as inputs in
#' the same form as distribution functions. The negative log-likelihood is a
#' wrapper for \code{\link[evmix:fbckden]{lbckden}}, designed towards making
#' it useable for optimisation (e.g. \code{lambda} given as first input).
#' 
#' The alternate bandwidth definitions are discussed in the
#' \code{\link[evmix:kernels]{kernels}}, with the \code{lambda} used here but 
#' \code{bw} also output. The \code{bw} specification is the same as used in the
#' \code{\link[stats:density]{density}} function.
#' 
#' The possible kernels are also defined in \code{\link[evmix:kernels]{kernels}} help
#' documentation with the \code{"gaussian"} as the default choice.
#' 
#' Unlike the standard KDE, there is no general rule-of-thumb bandwidth for all these
#' estimators, with only certain methods having a guideline in the literature, so none
#' have been implemented. Hence, a bandwidth must always be specified.
#' 
#' The \code{simple}, \code{renorm}, \code{beta1}, \code{beta2} \code{gamma1} and \code{gamma2}
#' density estimates require renormalisation, achieved
#' by numerical integration, so is very time consuming.
#' 
#' Missing values (\code{NA} and \code{NaN}) are assumed to be invalid data so are ignored.
#' 
#' Cross-validation likelihood is used for kernel density component, obtained by
#' leaving each point out in turn and evaluating the KDE at the point left out:
#'    \deqn{L(\lambda)\prod_{i=1}^{n} \hat{f}_{-i}(x_i)}
#' where 
#'    \deqn{\hat{f}_{-i}(x_i) = \frac{1}{(n-1)\lambda} \sum_{j=1: j\ne i}^{n} K(\frac{x_i - x_j}{\lambda})}
#' is the KDE obtained when the \eqn{i}th datapoint is dropped out and then 
#' evaluated at that dropped datapoint at \eqn{x_i}.
#' 
#' Normally for likelihood estimation of the bandwidth the kernel centres and
#' the data where the likelihood is evaluated are the same. However, when using
#' KDE for extreme value mixture modelling the likelihood only those data in the
#' bulk of the distribution should contribute to the likelihood, but all the
#' data (including those beyond the threshold) should contribute to the density
#' estimate. The \code{extracentres} option allows the use to specify extra
#' kernel centres used in estimating the density, but not evaluated in the
#' likelihood. The default is to just use the existing data, so
#' \code{extracentres=NULL}.
#' 
#' The default optimisation algorithm is "BFGS", which requires a finite negative 
#' log-likelihood function evaluation \code{finitelik=TRUE}. For invalid 
#' parameters, a zero likelihood is replaced with \code{exp(-1e6)}. The "BFGS" 
#' optimisation algorithms require finite values for likelihood, so any user 
#' input for \code{finitelik} will be overridden and set to \code{finitelik=TRUE} 
#' if either of these optimisation methods is chosen.
#'
#' It will display a warning for non-zero convergence result comes from 
#' \code{\link[stats:optim]{optim}} function call.
#' 
#' If the hessian is of reduced rank then the variance (from inverse hessian)
#' and standard error of bandwidth parameter cannot be calculated, then by default 
#' \code{std.err=TRUE} and the function will stop. If you want the bandwidth estimate
#' even if the hessian is of reduced rank (e.g. in a simulation study) then
#' set \code{std.err=FALSE}. 
#' 
#' @section Warning:
#' Two important practical issues arise with MLE for the kernel bandwidth:
#' 1) Cross-validation likelihood is needed for the KDE bandwidth parameter
#' as the usual likelihood degenerates, so that the MLE \eqn{\hat{\lambda} \rightarrow 0} as
#' \eqn{n \rightarrow \infty}, thus giving a negative bias towards a small bandwidth.
#' Leave one out cross-validation essentially ensures that some smoothing between the kernel centres
#' is required (i.e. a non-zero bandwidth), otherwise the resultant density estimates would always
#' be zero if the bandwidth was zero.
#' 
#' This problem occassionally rears its ugly head for data which has been heavily rounded,
#' as even when using cross-validation the density can be non-zero even if the bandwidth is zero.
#' To overcome this issue an option to add a small jitter should be added to the data
#' (\code{x} only) has been included in the fitting inputs, using the 
#' \code{\link[base:jitter]{jitter}} function, to remove the ties. The default options red in the 
#' \code{\link[base:jitter]{jitter}} are specified above, but the user can override these.
#' Notice the default scaling \code{factor=0.1}, which is a tenth of the default value in the
#' \code{\link[base:jitter]{jitter}}
#' function itself.
#' 
#' A warning message is given if the data appear to be rounded
#' (i.e. more than 5% of data are tied). If the estimated bandwidth is too small, then
#' data rounding is the likely culprit. Only use the jittering when the MLE of
#' the bandwidth is far too small. 
#' 
#' 2) For heavy tailed populations the bandwidth is positively biased, giving oversmoothing
#' (see example). The bias is due to the distance between the upper (or lower) order statistics not
#' necessarily decaying to zero as the sample size tends to infinity. Essentially, as the distance
#' between the two largest (or smallest) sample datapoints does not decay to zero, some smoothing between
#' them is required (i.e. bandwidth cannot be zero). One solution to this problem is to splice
#' the GPD at a suitable threshold to remove the problematic tail from the inference for the bandwidth, 
#' using the \code{\link[evmix:fbckdengpd]{fbckdengpd}} function for a heavy upper tail. See MacDonald et al (2013).
#' 
#' @return \code{\link[evmix:fbckden]{fbckden}} gives leave one out cross-validation
#' (log-)likelihood and 
#' \code{\link[evmix:fbckden]{lbckden}} gives the negative log-likelihood. 
#' \code{\link[evmix:fbckden]{nlbckden}} returns a simple list with the following elements
#'
#' \tabular{ll}{
#' \code{call}: \tab \code{optim} call\cr
#' \code{x}: \tab (jittered) data vector \code{x}\cr
#' \code{kerncentres}: actual kernel centres used \code{x}\cr
#' \code{init}: \tab \code{linit} for lambda\cr
#' \code{optim}: \tab complete \code{optim} output\cr
#' \code{mle}: \tab vector of MLE of bandwidth\cr
#' \code{cov}: \tab variance of MLE of bandwidth\cr
#' \code{se}: \tab standard error of MLE of bandwidth\cr
#' \code{nllh}: \tab minimum negative cross-validation log-likelihood\cr
#' \code{n}: \tab total sample size\cr
#' \code{lambda}: \tab MLE of lambda (kernel half-width)\cr
#' \code{bw}: \tab MLE of bw (kernel standard deviations)\cr
#' \code{kernel}: \tab kernel name\cr
#' \code{bcmethod}: \tab boundary correction method\cr
#' \code{proper}: \tab logical, whether renormalisation is requested\cr
#' \code{nn}: \tab non-negative correction method\cr
#' \code{offset}: \tab offset for log transformation method\cr
#' \code{xmax}: \tab maximum value of scale beta or copula
#' }
#' 
#' The output list has some duplicate entries and repeats some of the inputs to both 
#' provide similar items to those from \code{\link[evd:fpot]{fpot}} and to make it 
#' as useable as possible.
#'  
#' @note An initial bandwidth must be provided, so \code{linit} and \code{bwinit} 
#' cannot both be \code{NULL}
#' 
#' The extra kernel centres \code{extracentres} can either be a vector of data or \code{NULL}.
#' 
#' Invalid parameter ranges will give \code{0} for likelihood, \code{log(0)=-Inf} for
#' log-likelihood and \code{-log(0)=Inf} for negative log-likelihood. 
#' 
#' Infinite and missing sample values are dropped.
#' 
#' Error checking of the inputs is carried out and will either stop or give warning message
#' as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Kernel_density_estimation}
#' 
#' \url{http://en.wikipedia.org/wiki/Cross-validation_(statistics)}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Bowman, A.W. (1984). An alternative method of cross-validation for the smoothing of
#' density estimates. Biometrika 71(2), 353-360.
#' 
#' Duin, R.P.W. (1976). On the choice of smoothing parameters for Parzen estimators of
#' probability density functions. IEEE Transactions on Computers C25(11), 1175-1179.
#' 
#' MacDonald, A., Scarrott, C.J., Lee, D., Darlow, B., Reale, M. and Russell, G. (2011).
#' A flexible extreme value mixture model. Computational Statistics and Data Analysis
#' 55(6), 2137-2157.
#' 
#' MacDonald, A., C. J. Scarrott, and D. S. Lee (2011). Boundary correction, consistency
#' and robustness of kernel densities using extreme value theory. Submitted.
#' Available from: \url{http://www.math.canterbury.ac.nz/~c.scarrott}.
#' 
#' Wand, M. and Jones, M.C. (1995). Kernel Smoothing. Chapman && Hall.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}.
#'
#' @section Acknowledgments: Based on code
#' by Anna MacDonald produced for MATLAB.
#' 
#' @seealso \code{\link[evmix:kernels]{kernels}}, \code{\link[evmix:kfun]{kfun}},
#' \code{\link[base:jitter]{jitter}}, \code{\link[stats:density]{density}} and
#' \code{\link[stats:bandwidth]{bw.nrd0}}
#' 
#' @aliases fbckden lbckden nlbckden
#' @family  kden kdengpd kdengpdcon bckden bckdengpd bckdengpdcon
#'          fkden fkdengpd fkdengpdcon fbckden fbckdengpd fbckdengpdcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(1, 1))
#' 
#' nk=500
#' x = rgamma(nk, shape = 1, scale = 2)
#' xx = seq(-1, 10, 0.01)
#' 
#' # cut and normalize is very quick 
#' fit = fbckden(x, linit = 0.2, bcmethod = "cutnorm")
#' hist(x, nk/5, freq = FALSE) 
#' rug(x)
#' lines(xx, dgamma(xx, shape = 1, scale = 2), col = "black")
#' # but cut and normalize does not always work well for boundary correction
#' lines(xx, dbckden(xx, x, lambda = fit$lambda, bcmethod = "cutnorm"), lwd = 2, col = "red")
#' # Handily, the bandwidth usually works well for other approaches as well
#' lines(xx, dbckden(xx, x, lambda = fit$lambda, bcmethod = "simple"), lwd = 2, col = "blue")
#' lines(density(x), lty = 2, lwd = 2, col = "green")
#' legend("topright", c("True Density", "BC KDE using cutnorm",
#'   "BC KDE using simple", "KDE Using density"),
#'   lty = c(1, 1, 1, 2), lwd = c(1, 2, 2, 2), col = c("black", "red", "blue", "green"))
#' 
#' # By contrast simple boundary correction is very slow
#' # a crude trick to speed it up is to ignore the normalisation and non-negative correction,
#' # which generally leads to bandwidth being biased high
#' fit = fbckden(x, linit = 0.2, bcmethod = "simple", proper = FALSE, nn = "none")
#' hist(x, nk/5, freq = FALSE) 
#' rug(x)
#' lines(xx, dgamma(xx, shape = 1, scale = 2), col = "black")
#' lines(xx, dbckden(xx, x, lambda = fit$lambda, bcmethod = "simple"), lwd = 2, col = "blue")
#' lines(density(x), lty = 2, lwd = 2, col = "green")
#' 
#' # but ignoring upper tail in likelihood works a lot better
#' q75 = qgamma(0.75, shape = 1, scale = 2)
#' fitnotail = fbckden(x[x <= q75], linit = 0.1, 
#'    bcmethod = "simple", proper = FALSE, nn = "none", extracentres = x[x > q75])
#' lines(xx, dbckden(xx, x, lambda = fitnotail$lambda, bcmethod = "simple"), lwd = 2, col = "red")
#' legend("topright", c("True Density", "BC KDE using simple", "BC KDE (upper tail ignored)",
#'    "KDE Using density"),
#'    lty = c(1, 1, 1, 2), lwd = c(1, 2, 2, 2), col = c("black", "blue", "red", "green"))
#' }
#' 

# maximum cross-validation likelihood fitting for boundary corrected KDE
fbckden <- function(x, linit = NULL, bwinit = NULL, kernel = "gaussian", extracentres = NULL,
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL, xmax = NULL,
  add.jitter = FALSE, factor = 0.1, amount = NULL,
  std.err = TRUE, method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.quant(extracentres, allownull = TRUE, allowna = TRUE, allowinf = TRUE)
  check.kbw(linit, bwinit)
  check.kernel(kernel)
  check.bcmethod(bcmethod)
  check.logic(proper)
  check.nn(nn)
  check.offset(offset, bcmethod, allowzero = TRUE)
  check.posparam(xmax, allownull = TRUE)  
  check.posparam(factor)
  check.posparam(amount, allownull = TRUE)
  check.logic(add.jitter)
  check.logic(std.err)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x) 

  if (any(x < 0)) stop("data must be non-negative")
  
  if (!is.null(extracentres)) {
    if (any(!is.finite(extracentres))) warning("non-finite extracentres have been removed")

    extracentres = extracentres[is.finite(extracentres)]
    
    if (any(extracentres < 0)) stop("negative kernel centres not permitted")

    check.quant(extracentres)
  }

  if ((method == "L-BFGS-B") | (method == "BFGS")) finitelik = TRUE

  if (add.jitter) x = pmax(jitter(x, factor, amount), 0)
  
  xuniq = unique(x)
  if (length(xuniq) < (0.95*n))
    warning("data may be rounded, as more than 5% are ties, so bandwidth could be biased to zero")
  
  # if bcmethod does not use standard kernels then lambda must be specified
  # then bw can be used, but lambda should be defaulted to if available
  kernelmethods = c("simple", "cutnorm", "renorm", "reflect", "logtrans")
  if (!(bcmethod %in% kernelmethods)) {
    if (is.null(linit))
      stop(paste("bandwidth bw only relevant for", kernelmethods, collapse = " "))
  } else {
    linit = klambda(bwinit, kernel, linit)    
  }

  if ((bcmethod == "copula") & (linit >= 1))
    stop("bandwidth must between (0, 1) for copula method")  
    
  upboundmethods = c("beta1", "beta2", "copula")
  if (!is.null(xmax) & !(bcmethod %in% upboundmethods))
    warning(paste("xmax only relevant for boundary correction methods", upboundmethods, collapse = " "))
  
  if (bcmethod %in% upboundmethods) {
    if (is.null(xmax)) stop("xmax is NULL")
    
    if (max(x) > xmax) stop("largest kernel centre must be below xmax")

    if (any(x == 0)) {
      warning("kernel centres of zero are invalid for gamma or beta method so ignored")
      x = x[x != 0]
    }

    if ((bcmethod != "gamma1") & (bcmethod != "gamma2")) {
      if (any(x == xmax)) {
        warning("kernel centres of xmax are invalid for beta or copula method so ignored")
        x = x[x != xmax]
      }
    }
    # need to recheck there are some valid kernel centres after these exclusions
    check.quant(x)
    n = length(x)

    if (!is.null(extracentres)) {
      if (max(extracentres) > xmax) stop("largest kernel centre must be below xmax")   
    }
  }
    
  # It is not always easy to choose a sensible good initial value for lambda
  # So try adjusting lambda up and down a little to find a valid one to start off
  llhinit = lbckden(x, lambda = linit, kernel = kernel, extracentres = extracentres,
    bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)

  tryi = 0
  lfirst = linit
  
  # try upto 2^5 larger than original
  while (is.infinite(llhinit) & (tryi < 5)) {
    linit = linit*2
    llhinit = lbckden(x, lambda = linit, kernel = kernel, extracentres = extracentres,
      bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
    tryi = tryi + 1
  }
  
  # try upto 2^-5 smaller than original
  if (is.infinite(llhinit)) {
    tryi = 0
    linit = lfirst
    while (is.infinite(llhinit) & (tryi < 5)) {
      linit = linit/2
      llhinit = lbckden(x, lambda = linit, kernel = kernel, extracentres = extracentres, 
        bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
      tryi = tryi + 1
    }
  }
  
  if (is.infinite(llhinit))
    stop("likelihood is undefined for initial bandwidth try another value")  

  if (tryi != 0)
    warning(paste("initial bandwidth was invalid, so linit=", linit, "is used instead"))

  fit = optim(par = linit, fn = nlbckden, x = x, kernel = kernel, extracentres = extracentres, 
    bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax,
    finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)

  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == linit) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }

  lambda = fit$par
  if (bcmethod %in% kernelmethods) {
    bw = kbw(lambda, kernel)
  } else {
    bw = NA
  }
  
  if (conv & std.err) {
    se = sqrt(1/fit$hessian)
  } else {
    se = NULL
  }

  list(call = call, x = as.vector(x), kerncentres = c(x, extracentres), 
    init = as.vector(linit), optim = fit, conv = conv, cov = 1/fit$hessian, mle = fit$par,
    se = se, nllh = fit$value, n = n, lambda = lambda, bw = bw, kernel = kernel,
    bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
}

#' @export
#' @aliases fbckden lbckden nlbckden
#' @rdname  fbckden

# cross-validation log-likelihood function for boundary corrected KDE
# will not stop evaluation unless it has to
lbckden <- function(x, lambda = NULL, bw = NULL, kernel = "gaussian", extracentres = NULL,
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL, xmax = NULL, log = TRUE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.quant(extracentres, allownull = TRUE, allowna = TRUE, allowinf = TRUE)
  check.param(lambda, allownull = TRUE)
  check.param(bw, allownull = TRUE)
  check.kernel(kernel)
  check.bcmethod(bcmethod)
  check.logic(proper)
  check.nn(nn)
  check.offset(offset, bcmethod, allowzero = TRUE)
  check.posparam(xmax, allownull = TRUE)  
  check.logic(log)
  
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)

  if (any(x < 0)) stop("data must be non-negative")
  
  if (!is.null(extracentres)) {
    if (any(!is.finite(extracentres))) warning("non-finite extracentres have been removed")

    extracentres = extracentres[is.finite(extracentres)]
    
    if (any(extracentres < 0)) stop("negative kernel centres not permitted")
  }
  kerncentres = c(x, extracentres)
  check.quant(kerncentres)

  nk = length(kerncentres)
  
  # if bcmethod does not use standard kernels then lambda must be specified
  # then bw can be used, but lambda should be defaulted to if available
  kernelmethods = c("simple", "cutnorm", "renorm", "reflect", "logtrans")
  if (!(bcmethod %in% kernelmethods)) {
    if (is.null(lambda))
      stop(paste("bandwidth bw only relevant for", kernelmethods, collapse = " "))
  } else {
    if (is.null(lambda) & is.null(bw)) stop("lambda and bw cannot both be NULL")
    
    if (is.null(lambda)) lambda = klambda(bw, kernel)
  }

  upboundmethods = c("beta1", "beta2", "copula")
  if (!is.null(xmax) & !(bcmethod %in% upboundmethods))
    warning(paste("xmax only relevant for boundary correction methods", upboundmethods, collapse = " "))
  
  if (bcmethod %in% upboundmethods) {
    if (is.null(xmax)) stop("xmax is NULL")
    
    if (max(kerncentres) > xmax) stop("largest kernel centre must be below xmax")

    if (any(kerncentres == 0)) {
      warning("kernel centres of zero are invalid for gamma or beta method so ignored")
      kerncentres = kerncentres[kerncentres != 0]
    }

    if ((bcmethod != "gamma1") & (bcmethod != "gamma2")) {
      if (any(kerncentres == xmax)) {
        warning("kernel centres of xmax are invalid for beta or copula method so ignored")
        kerncentres = kerncentres[kerncentres != xmax]
      }
    }
    # need to recheck there are some valid kernel centres after these exclusions
    check.quant(kerncentres)

    if (!is.null(extracentres)) {
      if (max(extracentres) > xmax) stop("largest kernel centre must be below xmax")   
    }
  }

  dbckdeni <- function(i, kerncentres, lambda, kernel, bcmethod, proper, nn, offset, xmax) {
    di = dbckden(kerncentres[i], kerncentres[-i], lambda, kernel = kernel, 
      bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax, log = TRUE)
  }

  if ((lambda <= 0) | ((bcmethod == "copula") & (lambda >= 1)) |
      ((bcmethod == "beta1") & (lambda >= 0.25*ifelse(is.null(xmax), Inf, xmax))) | 
      ((bcmethod == "beta2") & (lambda >= 0.25*ifelse(is.null(xmax), Inf, xmax)))) {
    l = -Inf
  } else {
    l = sum(sapply(1:n, FUN = dbckdeni, kerncentres = kerncentres, lambda = lambda, kernel = kernel,
      bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax))
  }
  
  if (!log) l = exp(l)

  l
}

#' @export
#' @aliases fbckden lbckden nlbckden
#' @rdname  fbckden

# negative cross-validation log-likelihood function for boundary corrected KDE
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlbckden <- function(lambda, x, bw = NULL, kernel = "gaussian", extracentres = NULL,
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL, xmax = NULL, finitelik = FALSE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.quant(extracentres, allownull = TRUE, allowna = TRUE, allowinf = TRUE)
  check.param(lambda, allownull = TRUE)
  check.param(bw, allownull = TRUE)
  check.kernel(kernel)
  check.bcmethod(bcmethod)
  check.logic(proper)
  check.nn(nn)
  check.offset(offset, bcmethod, allowzero = TRUE)
  check.posparam(xmax, allownull = TRUE)  
  check.logic(finitelik)
  
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  nllh = -lbckden(x, lambda, kernel = kernel, extracentres = extracentres,
    bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax) 

  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
