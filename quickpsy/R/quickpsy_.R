#' Fits psychometric functions
#'
#' \code{quickpsy_} is the standard evaluation SE function associated
#' to the non-standard evaluation NSE function \code{quickpsy}.
#' \href{http://adv-r.had.co.nz/Computing-on-the-language.html}{SE functions can be more easily called from other functions.}
#' In SE functions, you need to quote the names of the variables.
#' @param d Data frame with the results of a Yes-No experiment to fit.
#' It should have a
#' \href{http://vita.had.co.nz/papers/tidy-data.html}{tidy} form in which
#' each column corresponds to a variable and each row is an observation.
#' @param x Name of the explanatory variable.
#' @param k Name of the response variable. The response variable could be the
#' number of trials in which a yes-type response was given or a vector of 0s
#' (or -1s; no-type response) and 1s (yes-type response) indicating the
#' response on each trial.
#' @param n Only necessary if \code{k} refers to the number of trials
#' in which a yes-type response was given. It corresponds to the name of the
#' variable indicating the total number of trials.
#' @param grouping Name of the grouping variables. It should be specified as
#' \code{grouping = .(variable_name1, variable_name2)}.
#' @param random Name of the random variable. It should be specified as
#' \code{random = .(variable_name1, variable_name2)}. In the current version
#' of quickpsy, the random variable has not special treatment. It does the
#' same as \code{grouping}.
#' @param within Name of the within variable. It should be specified as
#' \code{within = .(variable_name1, variable_name2)}. In the current version
#' of quickpsy, the within variable has not special treatment. It does the
#' same as \code{grouping}.
#' @param between Name of the between variable.  It should be specified as
#' \code{between = .(variable_name1, variable_name2)}. In the current version
#' of quickpsy, the between variable has not special treatment. It does the
#' same as \code{grouping}.
#' @param xmin Minimum value of the explanatory variable for which the curves
#' should be calculated (the default is the minimum value of the explanatory
#' variable).
#' @param xmax Maximum value of the explanatory variable for which the curves
#' should be calculated (the default is the maximum value of the explanatory
#' variable).
#' @param log If \code{TRUE}, the logarithm of the explanatory variable is used
#' to fit the curves (default is \code{FALSE}).
#' @param fun Name of the shape of the curve to fit. It could be a predefined
#' shape (\code{cum_normal_fun}, \code{logistic_fun}, \code{weibull_fun})
#' or the name of a function introduced by the user
#' (default is \code{cum_normal_fun}).
#' @param parini Initial parameters. quickpsy calculates default
#' initial parameters using probit analysis, but it is also possible to
#' specify a vector of initial parameters or a list of the form
#' \code{list(c(par1min, par1max), c(par2min, par2max))} to
#' constraint the lower and upper bounds of the parameters (when
#' \code{optimization = 'DE'}, parini should be also a list).
#' @param guess Value indicating the guess rate \eqn{\gamma} (default is 0). If
#' \code{TRUE}, the guess rate is estimated as the i + 1 paramEter where
#' i corresponds to the number of parameters of \code{fun}. If, for
#' example, \code{fun} is a predefined shape with parameters p1 and p2,
#' then the guess rate corresponds to parameter p3.
#' @param lapses Value indicating the lapse rate \eqn{\lambda} (default is 0).
#'  If \code{TRUE}, the lapse rate is estimated as the i + 1 parameter where
#' i corresponds to the number of parameters of \code{fun} plus one if
#' the guess rate is estimated. If, for example, \code{fun} is a
#' predefined shape with parameters p1 and p2,
#' then the lapse rate corresponds to parameter p3. If the guess rate is also
#' estimated, p3 will be the guess rate and p4 the lapse rate.
#' @param prob Probability to calculate the threshold (default is
#' \code{guess + .5 * (1 - guess)}).
#' @param thresholds If \code{FALSE}, thresholds are not calculated
#' (default is \code{TRUE}).
#' @param logliks If \code{TRUE}, the loglikelihoods are calculated
#'  (default is \code{FALSE}).
#' @param bootstrap \code{'parametric'} performs parametric bootstrap;
#' \code{'nonparametric'} performs non-parametric bootstrap;
#' \code{'none'} does not perform bootstrap (default is \code{'parametric'}).
#' @param B number of bootstrap samples (default is 100 ONLY).
#' @param ci Confidence intervals level based on percentiles (default is .95).
#' @param optimization Method used for optimizization. The default is 'optim' which uses
#' the \code{optim} function. It can also be \code{'DE'} which uses de function
#' \code{DEoptim} from the package DEoptim, which performs differential
#' evolution optimization. By using \code{DEoptim}, it is less likely that the
#' optimization finishes in a local minimum, but the optimization is slow.
#' When \code{'DE'} is used, \code{parini} should be specified as a list with
#' lower and upper bounds.
#' @seealso \code{\link{quickpsy}}
#' @examples
#' library(MPDiR) # contains the Vernier data
#' fit <- quickpsy_(Vernier, 'Phaseshift', 'NumUpward', 'N',
#'                 grouping = c('Direction', 'WaveForm', 'TempFreq'), B = 20)
#' plotcurves(fit)
#' @export
quickpsy_ <- function(d, x = 'x', k = 'k', n = 'n', grouping, random, within,
                      between, xmin = NULL, xmax = NULL, log = FALSE,
                      fun = 'cum_normal_fun', parini = NULL, guess = 0,
                      lapses = 0, prob = NULL, thresholds = T,  logliks = FALSE,
                      bootstrap = 'parametric', B = 100, ci = .95,
                      optimization = 'optim') {

  options(dplyr.print_max = 1e9)

  if (!is.null(prob)) thresholds <- T

  if (missing(n)) n <- NULL
  if (is.null(parini)) pariniset <- FALSE
  else pariniset <- TRUE

  cat('Estimating parameters...\n')
  qp <- fitpsy(d, x, k, n, random, within, between, grouping, xmin, xmax, log,
               fun, parini, pariniset, guess, lapses, optimization)

  qp <- c(qp, list(pariniset = pariniset))

  qp <- c(qp, list(ypred = ypred(qp)))
  if (sum(qp$ypred$ypred < 0) + sum(qp$ypred$ypred > 1) > 0)
    if (bootstrap == 'parametric')
      stop ('As y-predictions are not within (0,1), bootstrap should be \'nonparametric\'', call.=FALSE)

  qp <- c(qp, list(curves = curves(qp, xmin, xmax, log)))

  qp <- c(qp, list(sse=sse(qp)))

  if (thresholds) {
    if (is.null(prob))
      if (is.logical(guess) && guess) prob <- .5
      else  prob <- guess + .5 * (1 - guess)
    qp <- c(qp, list(thresholds = thresholds(qp, prob, log)))
  }

  if (logliks) qp <- c(qp, list(logliks = logliks(qp)))

  if (bootstrap == 'parametric' || bootstrap == 'nonparametric') {
    cat('Performing bootstrap...\n')
    qp <- c(qp, list(parbootstrap = parbootstrap(qp, bootstrap, B)))
    qp <- c(qp, list(parci = parci(qp, ci)))
     if (!(
       (length(qp$groups)==0) ||
       (length(qp$groups)==1 && nrow(unique(qp$averages[qp$groups]))==1)
       )) {
       qp <- c(qp, list(parcomparisons = parcomparisons(qp, ci)))
     }
    qp <- c(qp, list(curvesbootstrap = curvesbootstrap(qp, log = log)))
    if (thresholds) {
      qp <- c(qp,
              list(thresholdsbootstrap = thresholdsbootstrap(qp, prob, log)))
      qp <- c(qp, list(thresholdsci = thresholdsci(qp, ci)))

      if (!(
        (length(qp$groups)==0) ||
        (length(qp$groups)==1 && nrow(unique(qp$averages[qp$groups]))==1)
      )) {
        qp <- c(qp, list(thresholdcomparisons = thresholdcomparisons(qp, ci)))
      }
    }
  }
  else if (bootstrap != 'none')
    stop('Bootstrap should be \'parametric\', \'nonparametric\' or \'none\'.', call. = FALSE)

  if (log) qp$averages[[x]] <- exp(qp$averages[[x]])

  class(qp) <- 'quickpsy'
  qp
}


