#' Evaluation of geostatistical models of uncertainty
#' 
#' Evaluate the local quality of a geostatistical model of uncertainty 
#' (GMU) using summary measures and graphical displays.
#' 
#' @param observed Vector of observed values at the validation points. See 
#' \sQuote{Details} for more information.
#' 
#' @param simulated Data frame or matrix with simulated values (columns) for
#' each validation point (rows). See \sQuote{Details} for more information.
#' 
#' @param pi Vector defining the width of the series of probability intervals.
#' Defaults to \code{pi = seq(0.01, 0.99, 0.01)}. See \sQuote{Details} for more 
#' information.
#' 
#' @param symmetric Logical for choosing the type of probability interval.
#' Defaults to \code{symmetric = TRUE}. See \sQuote{Details} for more 
#' information.
#' 
#' @param plotit Logical for plotting the results. Defaults to 
#' \code{plotit = TRUE}.
#' 
#' @details 
#' There is no standard way of evaluating the local quality of a GMU. The 
#' collection of summary measures and graphical displays presented here is far
#' from being comprehensive. A few definitions are given bellow.
#' 
#' \subsection{Error statistics}{
#' Error statistics measure how well the GMU predicts the measured values at the
#' validation points. Four error statistics are presented:
#' 
#' \describe{
#' \item{Mean error (ME)}{
#' Measures the bias of the predictions of the GMU, being defined as the mean of
#' the differences between the average of the simulated values and the observed
#' values, i.e. the average of all simulations is taken as the predicted value.
#' }
#' \item{Mean squared error (MSE)}{
#' Measures the accuracy of the predictions of the GMU, being defined as the
#' mean of the squared differences between the average of the simulated values
#' and the observed values.
#' }
#' \item{Scaled root mean squared error (SRMSE)}{
#' Measures how well the GMU estimate of the prediction error variance (PEV)
#' approximates the observed prediction error variance, where the first is 
#' given by the variance of the simulated values, while the second is given by
#' the squared differences between the average of the simulated values, i.e. 
#' the squared error (SE). The SRMSE is computed as the average of SE / PEV, 
#' where SRMSE > 1 indicates underestimation, while SRMSE < 1 indicates
#' overestimation.
#' }
#' \item{Pearson correlation coefficient}{
#' Measures how close the GMU predictions are to the observed values. A scatter 
#' plot of the observed values versus the average of the simulated values 
#' can be used to check for possible unwanted outliers and non-linearities. 
#' The square of the Pearson correlation coefficient measures the fraction 
#' of the overall spread of observed values that is explained by the GMU, 
#' that is, the amount of variance explained (AVE), also known as coefficient 
#' of determination or ratio of scatter.
#' }
#' }
#' }
#' \subsection{Coverage probabilities}{
#' The coverage probability of an interval is given by the number of times that
#' that interval contains its parameter over several replications of an
#' experiment. For example, consider the interquartile range \eqn{IQR = Q3 - Q1}
#' of a Gaussian distributed variable with mean equal to zero and variance equal 
#' to one. The nominal coverage probability of the IQR is 0.5, i.e. two 
#' quarters of the data fall within the IQR. Suppose we generate a Gaussian 
#' distributed \emph{random} variable with the same mean and variance and count
#' the number of values that fall within the IQR defined above: about 0.5 of 
#' its values will fall within the IQR. If we continue generating Gaussian
#' distributed \emph{random} variables with the same mean and variance, on 
#' average, 0.5 of the values will fall in that interval.
#' 
#' Coverage probabilities are very useful to evaluate the local quality of a
#' GMU: the closer the observed coverage probabilities of a sequence of 
#' probability intervals (PI) are to the nominal coverage probabilities of 
#' those PIs, the better the modelling of the local uncertainty.
#' 
#' Two types of PIs can be used here: symmetric, median-centred PIs, and 
#' left-bounded PIs. Papritz & Dubois (1999) recommend using left-bounded PIs 
#' because they are better at evidencing deviations for both large and small 
#' PIs. The authors also point that the coverage probabilities of the symmetric,
#' median-centred PIs can be read from the coverage probability plots produced 
#' using left-bounded PIs.
#' 
#' In both cases, the PIs are computed at each validation location using the 
#' quantiles of the conditional cumulative distribution function (ccdf) defined
#' by the set of realizations at that validation location. For a sequence of 
#' PIs of increasing width, we check which of them contains the observed value 
#' at all validation locations. We then average the results over all validation 
#' locations to compute the proportion of PIs (with the same width) that 
#' contains the observed value: this gives the coverage probability of the PIs.
#' 
#' Deutsch (1997) proposed three summary measures of the coverage 
#' probabilities to assess the local \emph{goodness} of a GMU: accuracy ($A$),
#' precision ($P$), and goodness ($G$). According to Deutsch (1997), a GMU can 
#' be considered \dQuote{good} if it is both accurate and precise. Although 
#' easy to compute, these measures seem not to have been explored by many 
#' geostatisticians, except for the studies developed by Pierre Goovaerts and 
#' his later software implementation (Goovaerts, 2009). Richmond (2001) 
#' suggests that they should not be used as the only measures of the local 
#' quality of a GMU.
#' 
#' \describe{
#' \item{Accuracy}{
#' An accurate GMU is that for which the proportion \eqn{p^*} of true values 
#' falling within the $p$ PI is equal to or larger than the nominal probability
#' $p$, that is, when \eqn{p^* \geq p}. In the coverage probability plot, a GMU 
#' will be more accurate when all points are on or above the 1:1 line. The 
#' range of $A$ goes from 0 (lest accurate) to 1 (most accurate).
#' }
#' \item{Precision}{
#' The \emph{precision}, $P$, is defined only for an accurate GMU, and measures
#' how close \eqn{p^*} is to $p$. The range of $P$ goes from 0 (lest precise) to 1
#' (most precise). Thus, a GMU will be more accurate when all points in the 
#' PI-width plot are on or above the 1:1 line.
#' }
#' \item{Goodness}{
#' The \emph{goodness}, $G$, is a measure of the departure of the points from
#' the 1:1 line in the coverage probability plot. $G$ ranges from 0 (minimum
#' goodness) to 1 (maximum goodness), the maximum $G$ being achieved when 
#' \eqn{p^* = p}, that is, all points in both coverage probability and interval 
#' width plots are exactly on the 1:1 line.
#' }
#' }
#' It is worth noting that the coverage probability and PI-width plots are 
#' relevant mainly to GMU created using \emph{conditional simulations}, that is,
#' simulations that are locally conditioned to the data observed at the 
#' validation locations. Conditioning the simulations locally serves the 
#' purposes of honouring the available data and reducing the variance of the
#' output realizations. This is why one would like to find the points falling
#' above the 1:1 line in both coverage probability and PI-width plots. For
#' \emph{unconditional simulations}, that is, simulations that are only globally
#' conditioned to the histogram (and variogram) of the data observed at the
#' validation locations, one would expect to find that, over a large number 
#' of simulations, the whole set of possible values (i.e. the global histogram)
#' can be generated at any node of the simulation grid. In other words, it is 
#' expected to find all points on the 1:1 line in both coverage probability and 
#' PI-width plots. Deviations from the 1:1 line could then be used as evidence 
#' of problems in the simulation.
#' }
#' @references 
#' 
#' Deutsch, C. Direct assessment of local accuracy and precision. Baafi, E. Y. 
#' & Schofield, N. A. (Eds.) \emph{Geostatistics Wollongong '96}. Dordrecht:
#' Kinwer Academic Publishers, v. I, p. 115-125, 1997.
#' 
#' Papritz, A. & Dubois, J. R. Mapping heavy metals in soil by (non-)linear
#' kriging: an empirical validation. Gómez-Hernández, J.; Soares, A. & 
#' Froidevaux, R. (Eds.) \emph{geoENV II -- Geostatistics for Environmental
#' Applications}. Springer, p. 429-440, 1999.
#' 
#' Goovaerts, P. Geostatistical modelling of uncertainty in soil science.
#' \emph{Geoderma}. v. 103, p. 3 - 26, 2001.
#' 
#' Goovaerts, P. AUTO-IK: a 2D indicator kriging program for the automated 
#' non-parametric modeling of local uncertainty in earth sciences. 
#' \emph{Computers & Geosciences}. v. 35, p. 1255-1270, 2009.
#' 
#' Richmond, A. J. Maximum profitability with minimum risk and effort. Xie, H.;
#' Wang, Y. & Jiang, Y. (Eds.) \emph{Proceedings 29th APCOM}. Lisse: A. A.
#' Balkema, p. 45-50, 2001.
#' 
#' Ripley, B. D. \emph{Stochastic simulation}. New York: John Wiley & Sons, 
#' p. 237, 1987.
#' 
#' @note Comments by Pierre Goovaerts \email{pierre.goovaerts@@biomedware.com} 
#' were important to describe how to use the coverage probability and 
#' PI-width plots when a GMU is created using unconditional simulations.
#' 
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @export
#' @examples
#' \dontrun{
#' set.seed(2001)
#' observed <- round(rnorm(100), 3)
#' simulated <- t(
#'   sapply(1:length(observed), function (i) round(rnorm(100), 3)))
#' resa <- checkGMU(observed, simulated, symmetric = T)
#' resb <- checkGMU(observed, simulated, symmetric = F)
#' resa$error;resb$error
#' resa$goodness;resb$goodness
#' }
# FUNCTION #####################################################################
checkGMU <-
  function (observed, simulated, pi = seq(0.01, 0.99, 0.01),
            symmetric = TRUE, plotit = TRUE) {
    
    # Initial settings
    n_pts <- length(observed)
    n_pis <- length(pi)
    
    # If required, compute the symmetric probability intervals
    if (symmetric) {
      pi_bounds <- 
        sapply(1:length(pi), function (i) c(1 - pi[i], 1 + pi[i]) / 2)
      message("Processing ", n_pis, " symmetric probability intervals...")
    } else {
      message("Processing ", n_pis, " probability intervals...")
    }
    
    # Do true values fall into each of the (symmetric) probability intervals?
    fall <- matrix(nrow = n_pts, ncol = n_pis)
    width <- matrix(nrow = n_pts, ncol = n_pis)
    g_fall <- matrix(nrow = n_pts, ncol = n_pis)
    g_width <- matrix(nrow = n_pts, ncol = n_pis)
    if (symmetric) { # Deutsch (1997)
      for (i in 1:n_pts) {
        x <- simulated[i, ]
        y <- observed[i]
        for (j in 1:n_pis) {
          # Local
          bounds <- stats::quantile(x = x, probs = pi_bounds[, j])
          fall[i, j] <- bounds[1] < y & y <= bounds[2]
          width[i, j] <- as.numeric(bounds[2] - bounds[1])
          # Global
          g_bounds <- stats::quantile(x = observed, probs = pi_bounds[, j])
          g_fall[i, j] <- g_bounds[1] < y & y <= g_bounds[2]
          g_width[i, j] <- as.numeric(g_bounds[2] - g_bounds[1])
        }
      }
    } else { # Papritz & Dubois (1999)
      for (i in 1:n_pts) {
        x <- simulated[i, ]
        y <- observed[i]
        lowwer <- min(x)
        g_lowwer <- min(observed)
        for (j in 1:n_pis) {
          # Local
          upper <- stats::quantile(x = x, probs = pi[j])
          fall[i, j] <- y <= upper
          width[i, j] <- as.numeric(upper - lowwer)
          # Global
          g_upper <- stats::quantile(x = observed, probs = pi[j])
          g_fall[i, j] <- y <= g_upper
          g_width[i, j] <- as.numeric(g_upper - g_lowwer)
        }
      }
    }
    
    # Compute the proportion of true values that fall into each of the 
    # (symmetric) probability intervals
    count <- apply(fall, 2, sum)
    prop <- count / n_pts
    
    g_count <- apply(g_fall, 2, sum)
    # g_prop <- g_count / n_pts
    
    # Compute the average width of the (symmetric) probability intervals into
    # each the true values fall
    width <- width * fall
    width <- apply(width, 2, sum) / count
    
    g_width <- g_width * g_fall
    g_width <- apply(g_width, 2, sum) / g_count
    
    # Compute summary statistics
    accu <- prop >= pi
    pi_idx <- which(accu)
    accu <- sum(prop >= pi) / n_pis # accuracy
    prec <- 1 - 2 * sum(prop[pi_idx] - pi[pi_idx]) / n_pis # precision
    pi_w <- ifelse(1:n_pis %in% pi_idx, 1, 2)
    good <- 1 - (sum(pi_w * abs(prop - pi)) / n_pis) # goodness
    pred <- apply(simulated, 1, mean) # predicted value
    pred_var <- apply(simulated, 1, stats::var) # prediction variance
    err <- pred - observed # error
    me <- mean(err) # mean error
    serr <- err ^ 2 # squared error
    mse <- mean(serr) # mean squared error
    srmse <- mean(serr / pred_var) # scaled root mean squared error
    corr <- stats::cor(pred, observed) # linear correlation
    error_stats <- data.frame(me = me, mse = mse, srmse = srmse, cor = corr)
    good_meas <- data.frame(A = accu, P = prec, G = good, symmetric = symmetric)
    
    if (plotit) {
      on.exit(graphics::par())
      graphics::par(mfrow = c(2, 2))
      cex <- ifelse(n_pts > 10, 0.5, 1)
      
      # Coverage probability plot
      graphics::plot(
        0:1, 0:1, type = 'n', main = "Coverage probability",
        xlab = "Probability interval", ylab = "Proportion")
      graphics::abline(a = 0, b = 1)
      graphics::points(x = pi, y = prop, cex = cex)
      if (symmetric) {
        graphics::text(x = 1, y = 0, labels = "Symmetric PIs", pos = 2)
      }
      
      # PI-width plot
      lim <- range(c(width, g_width), na.rm = TRUE)
      graphics::plot(
        x = width, y = g_width, ylim = lim, xlab = "Local", 
        ylab = "Global", cex = cex, xlim = lim, main = "PI width")
      graphics::abline(a = 0, b = 1)
      if (symmetric) {
        graphics::text(x = lim[2], y = lim[1], labels = "Symmetric PIs", pos = 2)
      }
      
      # Plot observed vs simulated values
      lim <- range(c(observed, pred))
      graphics::plot(
        x = observed, pred, main = "Observed vs Simulated", xlab = "Observed",
        ylim = lim, xlim = lim, ylab = "Simulated (average)", cex = cex)
      graphics::abline(a = 0, b = 1)
      
      # Plot box plots
      idx <- 1:n_pts
      idx <- idx[order(rank(observed))]
      if (n_pts > 100) { 
        sub_idx <- round(seq(1, n_pts, length.out = 100))
        graphics::boxplot(
          t(simulated[idx[sub_idx], ]), col = "yellow", pars = list(cex = cex),
          names = idx[sub_idx])
        graphics::points(
          observed[idx[sub_idx]], col = "red", pch = 17, cex = cex)
        xlab <- "Validation point (max of 100)"
      } else {
        graphics::boxplot(
          t(simulated[idx, ]), col = "yellow", pars = list(cex = cex),
          names = idx, xlab = "Validation point")
        graphics::points(observed[idx], col = "red", pch = 17, cex = cex)
        xlab <- "Validation point"
      }
      graphics::title(
        main = "Distribution of values", xlab = xlab, ylab = "Distribution")
    }
    
    # Output
    res <- list(data = data.frame(pi = pi, prop = prop, width = width), 
                error = error_stats, goodness = good_meas)
    return (res)
    
  }
