#' Initial covariance parameters (ICP)
#' 
#' Guess the initial values for the covariance parameters required to fit a
#' variogram model.
#' 
#' @inheritParams vgmLags
#' 
#' @param z Numeric vector with the values of the response variable for which 
#' the initial values for the covariance parameters should be guessed.
#' 
#' @param lags Numeric scalar defining the width of the lag-distance classes,
#' or a numeric vector with the lower and upper bounds of the lag-distance
#' classes. If missing, the lag-distance classes are computed using 
#' \code{\link[pedometrics]{vgmLags}}. See \sQuote{Details} for more 
#' information.
#' 
#' @param method Character keyword defining the method used for guessing the
#' initial covariance parameters. Defauls to \code{method = "a"}. See 
#' \sQuote{Details} for more information.
#' 
#' @param min.npairs Positive integer defining the minimum number of 
#' point-pairs required so that a lag-distance class is used for guessing the
#' initial covariance parameters. Defaults to \code{min.npairs = 30}.
#' 
#' @param model Character keyword defining the variogram model that will be
#' fitted to the data. Currently, most basic variogram models are accepted.
#' See \code{\link[geoR]{cov.spatial}} for more information. Defaults to
#' \code{model = "matern"}.
#' 
#' @param nu numerical value for the additional smoothness parameter \eqn{\nu} 
#' of the correlation function. See \code{\link[RandomFields]{RMmodel}} and
#' argument \code{kappa} of \code{\link[geoR]{cov.spatial}} for more 
#' information.
#' 
#' @param plotit Should the guessed initial covariance parameters be plotted
#' along with the sample variogram? Defaults to \code{plotit = FALSE}.
#' 
#' @param estimator Character keyword defining the estimator for computing the
#' sample variogram, with options \code{"qn"}, \code{"mad"}, \code{"matheron"},
#' and \code{"ch"}. Defaults to \code{estimator = "qn"}. See 
#' \code{\link[georob]{sample.variogram}} for more details.
#' 
#' @return A vector of numeric values: the guesses for the covariance parameters
#' nugget, partial sill, and range.
#' 
#' @details There are five methods two guess the initial covariance parameters
#' (ICP). Two of them (\code{"a"} and \code{"b"}) rely a sample variogram with
#' exponentially spaced lag-distance classes, while the other three (\code{"b"},
#' \code{"d"}, and \code{"e"}) use equidistant lag-distance classes (see
#' \code{\link[pedometrics]{vgmLags}}). All of them are 
#' \href{https://en.wikipedia.org/wiki/Heuristic}{heuristic}.
#' 
#' Method \code{"a"} was developed in-house, and is the most elaborated of them,
#' specially for guessing the nugget variance. Method \code{"c"} is implemented 
#' in the \pkg{automap}-package and was developed by 
#' \href{http://dx.doi.org/10.1016/j.cageo.2008.10.011}{Hiemstra et al. (2009)}.
#' 
#' Method \code{"b"} was proposed by 
#' \href{http://dx.doi.org/10.1016/0098-3004(95)00095-X}{Jian et al. (1996)} and
#' is implemented in \href{https://support.sas.com/documentation/cdl/en/statug/63347/HTML/default/viewer.htm#statug_variogram_a0000000593.htm}{SAS/STAT(R) 9.22}.
#' Method \code{"d"} was developed by 
#' \href{http://dx.doi.org/10.1007/s11004-012-9434-1}{Desassis & Renard (2012)}.
#' Method \code{"e"} was proposed by 
#' \href{http://www.ccgalberta.com/ccgresources/report05/2003-122-varfit.pdf}{Larrondo et al. (2003)} and is implemented in the VARFIT module of 
#' \href{http://www.gslib.com/}{GSLIB}.
#' 
#' @references 
#' 
#' Desassis, N. & Renard, D. Automatic variogram modelling by iterative least
#' squares: univariate and multivariate cases. \emph{Mathematical Geosciences}.
#' Springer Science \eqn{+} Business Media, v. 45, p. 453-470, 2012.
#' 
#' Hiemstra, P. H.; Pebesma, E. J.; Twenhöfel, C. J. & Heuvelink, G. B. 
#' Real-time automatic interpolation of ambient gamma dose rates from the Dutch 
#' radioactivity monitoring network. \emph{Computers & Geosciences}. Elsevier 
#' BV, v. 35, p. 1711-1721, 2009.
#' 
#' Jian, X.; Olea, R. A. & Yu, Y.-S. Semivariogram modelling by weighted least
#' squares. \emph{Computers & Geosciences}. Elsevier BV, v. 22, p. 387-397, 
#' 1996.
#' 
#' Larrondo, P. F.; Neufeld, C. T. & Deutsch, C. V. \emph{VARFIT: a program for 
#' semi-automatic variogram modelling}. Edmonton: Department of Civil and
#' Environmental Engineering, University of Alberta, p. 17, 2003.
#' 
#' @author Alessandro Samuel-Rosa <\email{alessandrosamuelrosa@@gmail.com}>
#' 
#' @seealso \code{\link[pedometrics]{vgmLags}}, 
#'          \code{\link[georob]{sample.variogram}}, 
#'          \code{\link[automap]{autofitVariogram}}
#' 
#' @concept variogram
#' @export
#' 
#' @examples
#' data(meuse, package = "sp")
#' icp <- vgmICP(z = log(meuse$copper), coords = meuse[, 1:2])
# FUNCTION - MAIN ##############################################################
vgmICP <- 
  function (z, coords, lags, cutoff = 0.5, method = "a", min.npairs = 30, 
            model = "matern", nu = 0.5, estimator = "qn", plotit = FALSE) {
    
    # Check arguments
    cov_models <- c(
      "matern", "exponential", "gaussian", "spherical", "circular", "cubic",
      "wave", "linear", "power", "powered.exponential", "stable", "cauchy",
      "gencauchy", "gneiting", "gneiting.matern", "pure.nugget")
    if (!model %in% cov_models) {
      stop (paste("model '", model, "' is not implemented", sep = ""))
    }
    
    # Check if suggested packages are installed
    pkg <- c("georob", "geoR")
    id <- !sapply(pkg, requireNamespace, quietly = TRUE)
    if (any(id)) {
      pkg <- paste(pkg[which(id)], collapse = " ")
      stop(paste("Package(s) needed for this function to work but not",
                 "installed: ", pkg, sep = ""), call. = FALSE)
    }
    
    # Check lags and max.dist
    if (missing(lags)) {
      if (method %in% c("a", "c")) {
        lags <- vgmLags(coords, cutoff = cutoff)
      } else {
        lags <- vgmLags(coords, n.lags = 15, type = "equi", cutoff = cutoff)
      }
    }
    cutoff <- sqrt(
      sum(apply(apply(coords[, 1:2], 2, range), 2, diff) ^ 2)) * cutoff
    
    # Compute variogram
    v <- georob::sample.variogram(
      object = z, locations = coords, lag.dist.def = lags, max.lag = cutoff,
      estimator = estimator)
    lags0 <- length(v$lag.dist)
    
    # Merge lag-distance classes that have too few point-pairs
    if (any(v$npairs < min.npairs)) {
      message("correcting lags for minimum number of point-pairs...")
      idx <- which(v$npairs < min.npairs) + 1
      while (length(idx) >= 1) {
        lags <- lags[-idx[1]]
        v <- georob::sample.variogram(
          object = z, locations = coords, lag.dist.def = lags, 
          max.lag = cutoff, estimator = estimator)
        idx <- which(v$npairs < min.npairs) + 1
      }
    }
    if (plotit) {
      graphics::plot(v)
    }
    
    # SILL
    # The initial guess for the sill should be the easiest among all
    # parameters needed to fit a covariance model. Several rules are used in the
    # literature:
    # 1) the average of the maximum semivariance and the median semivariance in 
    #    the sample variogram (Hiemstra et al., 2009);
    # 2) the average of all the experimental points (Larrondo et al., 2003);
    # 3) the average  semivariance of the three last lag-distance classes
    #    (Jian et al., 1996)
    # 4) the total variance (Desassis & Renard, 2012).
    # I define the partial sill as the difference between the variance
    # of the data minus the nugget variance.
    sill <- switch(
      method,
      a = { # Samuel-Rosa
        # The average of the variance of the data and gamma at the two last 
        # lag-distance classes.
        mean(c(stats::var(z), v$gamma[c(length(v$gamma), length(v$gamma) - 1)]))
      },
      b = { # JianEtAl1996
        mean(v$gamma[seq(length(v$gamma) - 2, length(v$gamma))])
      },
      c = { # HiemstraEtAl2009
        mean(c(max(v$gamma), stats::median(v$gamma)))
      },
      d = { # DesassisEtAl2012
        stats::var(z)
      },
      e = { # LarrondoEtAl2003
        mean(v$gamma)
      }
    )
    if (plotit) {
      graphics::abline(h = sill, lty = "dashed")
    }
    
    # RANGE
    # In general, the initial guess for the range (scale) parameter is made 
    # based on the lag-distance classes and on the dimensions of the study area.
    # The most common rule is to compute the initial range as half the 
    # maximum distance up to which lag-distance classes have been defined 
    # (Jian et al., 1996; Larrondo et al., 2003; Desassis & Renard, 2012). 
    # Others set the initial range to a proportion of the diagonal of the study 
    # area, say 0.1, 0.35 or 0.5 (Hiemstra et al., 2009).
    # I think that this rather arbitrary and, possibly, suboptimal because the 
    # lag-distance classes usually are defined by some automatic procedure
    # implemented in the software being used, which does not account for the
    # features of the data that is being analysed.
    # I propose using the estimate of the variance and semivariance of the data
    # to make an initial guess for the range parameter. The variance is used
    # here because it is the initial guess of the total sill (see bellow).
    # I start computing the absolute difference between the semivariogram and 
    # the variance in each lag-distance class (except for the first). Then, I
    # record the index of the lag-distance class where the semivariance is 
    # closest to the variance. The separation distance at the centre of this 
    # lag-distance class is used as the initial guess for the range parameter.
    range <- switch(
      method, 
      a = { # Samuel-Rosa2015
        v$lag.dist[which.min(abs(v$gamma[-c(1, length(v$gamma))] - sill)) + 1]
      }, 
      b = { # JianEtAl1996
        lags[length(lags)] * 0.5 
      },
      c = { # HiemstraEtAl2009
        sqrt(sum(apply(apply(coords, 2, range), 2, diff) ^ 2)) * 0.1
      },
      d = { # DesassisEtAl2012
        lags[length(lags)] * 0.5
      },
      e = { # LarrondoEtAl2003
        lags[length(lags)]
      }
    )
    
    # Correct initial guess of the range parameter
    range <- range / 
      geoR::practicalRange(cov.model = model, phi = 1, kappa = nu)
    
    if (plotit) {
      graphics::abline(v = range, lty = "dashed")
    }
    
    # NUGGET
    # The initial guess for the nugget variance is commonly made using one of
    # the following rules:
    # 1) use the minimum semivariance in the sample variogram (Hiemstra et al.,
    #    2009)
    # 2) set the initial nugget value to zero (Larrondo et al., 2003)
    # We can also find rules that take into account the difference in the
    # semivariance between the first and second lag-distance classes weighted
    # by the difference in the size of these lag-distance classes (Jian et al.,
    # 1996). The resulting initial guess for the nugget variance is always lower
    # than the minimum semivariance value.
    nugget <- switch(
      method,
      a = { # Samuel-Rosa
        
        # Is the minimum gamma in lags others than the first?
        # 
        if (which.min(v$gamma) > 1) {
          # 
          # If the minimum gamma is in lags other than in the first, then we 
          # may have one of two possibilities at hand. First, it may be that
          # the best covariance model is that of a pure nugget effect. Second,
          # the data at hand is not appropriate to estimate the behaviour of 
          # the variogram close to the origin. So, what is the best initial
          # guess?
          # 
          # The first task is to check which of the first three lags has the
          # lowest gamma.
          # 
          if (which.min(v$gamma[1:3]) == 2) {
            # 
            # If the second lag has the lowest gamma, we have to find out if 
            # gamma in the first lag is higher than in the third lag.
            # 
            if (v$gamma[1] > v$gamma[3]) {
              #
              # When gamma is the first lag is higher than gamma in the third 
              # lag, the best initial guess is the average of gamma in the 
              # second and third lags.
              # 
              mean(v$gamma[c(2, 2, 3)])
            } else {
              #
              # When gamma in the first lag is lower than in the third lag, then
              # the best initial guess is the average of gamma in the first and
              # second lags.
              # 
              mean(v$gamma[c(1, 2, 2)])
            }
          } else {
            # 
            # The third lag has the lowest gamma.
            # 
            if (v$gamma[1] > v$gamma[2]) {
              mean(v$gamma[c(2, 3, 3)])
            } else {
              mean(v$gamma[c(1, 3, 3)])
            }
          }
          # 
          # The first task is to check if gamma in the first lag is larger than
          # gamma in the second lag. This could easily happen if we have a 
          # pure nugget effect.
          # Is the semivariance of the fist lag larger than that of the second
          # lag? Or is the semivariance of the fist lag larger than the 
          # variance of z? Or is the semivariance of the fist and second lags
          # larger than that of the third? This looks like a messy sample
          # variogram.
          # if (v$gamma[1] >= v$gamma[2] || v$gamma[1] >= var(z) || 
          # which.min(v$gamma[1:3]) == 3) {
          # 
          # Well... Let us simply use the minimum gamma value.
          # min(v$gamma[1:3])
          # 
          # }
          # min(v$gamma)
          
        } else {
          # The first lag-distance class holds the minimum gamma!
          # 
          # The task now is to estimate the shape of the semivariogram near the
          # origin. The first thing to do is to test if gamma at the second lag
          # is closer to the estimate of the sill, or to gamma at the first lag.
          # If gamma at the second lag is closer to the estimated sill, then it 
          # means that the sample variogram is very steep near the origin.
          # Because the distance between these two lags and the difference in 
          # their sizes are very small, it may be that gamma in the first 
          # lag-distance class is spuriously low (or not).
          # 
          # if (diff(v$gamma[1:2]) > abs(diff(c(v$gamma[2], var(z))))) {
          if (diff(v$gamma[1:2]) > abs(diff(c(v$gamma[2], sill)))) {
            # Gamma in the second lag is closer to the estimated sill than to 
            # gamma in the first lag.
            # 
            # The task now is to evaluate the third lag in relation to the first
            # two. In the ideal situation, gamma will increase monotonically 
            # from the first to the third lag. If not, i.e. gamma in the third
            # lag is lower than in the second lag, then our data is not 
            # appropriate to estimate the behaviour of the variogram near the 
            # origin: gamma in the first lag is too low, and gamma at the second
            # lag is too high (higher than gamma at the third lag). It may be
            # better to use a conservative initial guess, such as the average 
            # of gamma at the first and third lags.
            # 
            if (!identical(order(v$gamma[1:3]), 1:3)) {
              # mean(c(v$gamma[1], min(v$gamma[-1])))
              mean(v$gamma[c(1, 3)])
              
            } else {
              # Gamma increases monotonically from the first to the third lag,
              # but gamma at the first lag is too low. It may be better to use 
              # a conservative initial guess, such as the average of gamma at 
              # the first and second lags.
              # 
              # mean(v$gamma[c(1, 2)])
              mean(v$gamma[c(1, 1, 2)])
              
            }
          } else {
            # Gamma in the second lag is closer to gamma in the first lag than
            # to the estimated sill.
            # 
            # This is the most wanted situation because it should enable 
            # estimating the behaviour of the variogram close to the origin with
            # great precision. The question is whether we should use gamma in 
            # the fist lag as the initial guess, or a value lower than that. If
            # the variogram is very steep, then it could be more appropriate to
            # use a value lower than gamma in the first lag as our initial 
            # guess. We could use the estimated range to decide on this matter.
            # a) If the estimated range is lower than the mean lag distance of 
            # the third lag class, then it is likely that the sample variogram 
            # is very steep near the origin (exponential model?). Thus, the 
            # initial guess for the nugget is half the difference between zero 
            # and gamma at the first lag.
            # b) Otherwise, it is likely that the sample variogram rises slowly
            # near the origin (Gaussian model?). The initial guess for the
            # nugget is gamma at the first lag.
            # 
            if (range < v$lag.dist[2]) {
              if (lags0 == length(v$lag.dist)) {
                diff(c(0, v$gamma[1])) * 0.5
              } else {
                v$gamma[1]
              }
            } else {
              v$gamma[1]
            }
          }
        }
      },
      b = { # JianEtAl1996
        max(0, c(v$gamma[1] - (lags[1] / diff(lags[1:2])) * diff(v$gamma[1:2])))
      },
      c = { # HiemstraEtAl2009
        min(v$gamma)
      },
      d = { # DesassisEtAl2012
        1e-12
      },
      e = { # LarrondoEtAl2003
        1e-12
      }
    )
    if (plotit) {
      graphics::abline(h = nugget, lty = "dashed")
    }
    
    # Guess the partial sill for a pure nugget effect model
    if (nugget < sill) {
      p_sill <- sill - nugget
    } else {
      p_sill <- 1e-12
    }
    
    # Prepare output
    res <- c(range = range, p_sill = p_sill, nugget = nugget)
    return (res)
  }
