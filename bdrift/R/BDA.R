#' @title Beta Drift Anaylsis
#'
#' @description
#' \code{BDA} performs a range of parameter instability diagnostics
#' for financial multi-factor models and returns data frames containing the
#' drifting parameters and their standard errors, a list of summary statistics
#' and an overview plot for each factor.
#'
#' @details
#' \code{BDA} performs a threefold analyis of a user-specified baseline model.
#' First, \code{BDA} performs a rolling regression across the entire data frame
#' where \code{horizon} determines the regression window size. The function includes
#' all rolling parameter estimates and standard errors in the output, so users
#' can access them using \code{$tdrift} and \code{$tdrift.se} respectively.
#'
#' Second, \code{BDA} estimates the baseline model parameters with estimation
#' windows of varying length from (\code{min.hor} to \code{max.hor}). Users can access
#' the resulting parameter estimates and standard errors using \code{$hdrift} and
#' \code{$hdrift.se} respectively.
#'
#' Third, \code{BDA} checks the baseline model for observations that have a noteworthy 
#' impact on the parameter estimate.
#'
#' For further details on the summary statistics output and plotting, please
#' reference \code{\link{summary.BDA}} and \code{\link{plot.BDA}} respectively.
#'
#' Although \code{BDA} was primarily developed to analyze financial multi-factor
#' models, it is capable to analyze any model fit, as long as the underlying data
#' is of class \code{xts}. However, \code{BDA} was developed with large datasets
#' in mind, so that very small datasets might produce errors or non-sensical results.
#'
#'
#' @param data an xts object containing all relevant time series. Please note that
#' \code{BDA} assumes that the data is ordered so that the most recent data is at the
#' tail of the matrix.
#' @param spec contains the formula for the baseline model.
#' @param horizon the time period for which the paramters should be estimated.
#' (e.g. 250 for a year, assuming daily data). By default, half of the data length is used.
#' @param min.hor the minimum horizon used in the analysis, by default one month,
#' assuming daily data (21 obs if available).
#' @param max.hor the maximum horizon used in the analysis, by default three years,
#' assuming daily data (750 obs if available)
#' @param family type of regression family passed to the \code{glm} function.
#' For further details on family types refer to \code{\link{family}}. Please note
#' that at this point the built-in plotting does not support all families.
#' @param doplot logical. If \code{TRUE}, the function returns diagnostic
#' plots for each parameter.
#' @param ... aditional commands passed to the \code{glm} function.
#' @export
#' @return a list with 8 elements:
#' \item{CALL}{function call}
#' \item{base.model}{baseline model}
#' \item{tdrift}{xts matrix containing historical estimates of baseline model}
#' \item{tdrift.se}{xts matrix containing historical standard errors of baseline model}
#' \item{hdrift}{matrix containing estimates of baseline model
#' with varying horizon lengths}
#' \item{hdrift.se}{matrix containing standard errors of
#' baseline model with varying horizon lengths}
#' \item{jackknife}{jackknife procedure of object class lm.influence}
#' \item{sumstats}{list containing various summary statistics}
#'
#' @author Markus Peter Auer <mp.auer@@meanerreversion.com>
#' @examples
#' \dontrun{
#' ###################################################
#' ####   3-Factor Stock Example: ExxonMobil      ####
#' ###################################################
#' 
#' results1 <- BDA(data = FFfactors, 
#'                 spec = (XOM~Mkt.RF + SMB + HML),
#'                 horizon = 250, doplot = TRUE)
#' 
#' ###################################################
#' #### 5-Factor Active Fund Example: BlackRock   ####
#' ###################################################
#' 
#' results2 <- BDA(data = FFfactors, 
#'                 spec = (MDLRX~Mkt.RF + SMB + HML + RMW + CMA),
#'                 horizon = 250, doplot = TRUE)
#' 
#' ###################################################
#' ####   1-Factor Index Fund Example: Vanguard   ####
#' ###################################################
#' 
#' results3 <- BDA(data = FFfactors, spec = (VOO~SP500),
#'                 horizon = 250, doplot = FALSE)
#' }
#' ###################################################
#' ####        CRAN-compatible example            ####
#' ###################################################
#' 
#' results <- BDA(data = FFfactors[nrow(FFfactors):(nrow(FFfactors)-300),], 
#'                spec = (VOO~SP500),horizon = 250, doplot = TRUE)
#' message("NOTE: This is a shortened example. Reference the manual for more complex examples")
#' 


BDA <- function(data, spec, horizon = round(nrow(data)*0.5)-1,
                min.hor = 21, max.hor = 750, family = gaussian,
                doplot = TRUE, ...)
{
  N <- nrow(data)-horizon+1
  ind <- zoo::index(data)[horizon:(horizon+N-1)]
  data <- as.data.frame(data)
  data <- data[rev(rownames(data)),]
  base.model <- glm(spec, data=data[1:horizon,], family = family, ...)
  k <- length(base.model$coefficients)
  tdrift <- matrix(0, nrow=N, ncol = k)
  tdrift.se <- matrix(0, nrow=N, ncol = k)
  min.hor <- ifelse(min.hor>nrow(data), nrow(data), min.hor)
  max.hor <- ifelse(max.hor<horizon, horizon+1, max.hor)
  max.hor <- ifelse(max.hor>nrow(data), nrow(data)-1, max.hor)
  hdrift <- matrix(0, nrow = max.hor - min.hor+1, ncol = k)
  hdrift.se <- matrix(0, nrow = max.hor - min.hor+1, ncol = k)
  coef <- base.model$coef[1:k]
  se <- sqrt(diag(vcov(base.model)))[1:k]
  header <- c(names(coef[1:k]))
  
  ###################################################
  ####              Beta Drift                   ####
  ###################################################
  for (i in 1:N)
  {
    tdrift.model <- glm(spec, data = data[i:(horizon+i-1),], family = family, ...)
    tdrift[i,] <- tdrift.model$coef[1:k]
    tdrift.se[i,] <- sqrt(diag(vcov(tdrift.model)))[1:k]
  }
  colnames(tdrift) <- header[1:k]
  tdrift <- as.data.frame(tdrift)
  tdrift <- tdrift[rev(rownames(tdrift)),]
  tdrift <- xts::as.xts(tdrift, order.by = ind)
  colnames(tdrift.se) <- header[1:k]
  tdrift.se <- as.data.frame(tdrift.se)
  tdrift.se <- tdrift.se[rev(rownames(tdrift.se)),]
  tdrift.se <- xts::as.xts(tdrift.se, order.by = ind)

  ###################################################
  ####            Testing Range                  ####
  ###################################################
  for (m in min.hor:max.hor)
  {
    hdrift[(m +1 - min.hor),] <- glm(spec, data = data[1:m,],
                                    family= family, ...)$coef[1:k]
    hdrift.se[(m +1 - min.hor),] <- sqrt(diag(vcov(glm(spec, data = data[1:m,],
                                                     family= family, ...))))[1:k]
  }
  rownames(hdrift) <- c(min.hor:max.hor)
  colnames(hdrift) <- header[1:k]
  rownames(hdrift.se) <- c(min.hor:max.hor)
  colnames(hdrift.se) <- header[1:k]

  ###################################################
  ####                Jackknife                  ####
  ###################################################
  jackknife <- lm.influence(base.model)
  
  ###################################################
  ####            Summary Statistics             ####
  ###################################################
  base.para <- matrix(0, nrow = 1, ncol = k)
  base.coef <- matrix(0, nrow = 1, ncol = k)
  base.se <- matrix(0, nrow = 1, ncol = k)
  base.p <- matrix(0, nrow = 1, ncol = k)
  base.lower <- matrix(0, nrow = 1, ncol = k)
  base.upper <- matrix(0, nrow = 1, ncol = k)
  tdrift.mean <- matrix(0, nrow = 1, ncol = k)
  tdrift.median <- matrix(0, nrow = 1, ncol = k)
  tdrift.sd <- matrix(0, nrow = 1, ncol = k)
  tdrift.lower <- matrix(0, nrow = 1, ncol = k)
  tdrift.upper <- matrix(0, nrow = 1, ncol = k)
  max.beta <- matrix(0, nrow = 1, ncol = k)
  pos1 <- matrix(0, nrow = 1, ncol = k)
  min.beta <- matrix(0, nrow = 1, ncol = k)
  pos2 <- matrix(0, nrow = 1, ncol = k)
  hor.mean <- matrix(0, nrow = 1, ncol = k)
  hor.sd <- matrix(0, nrow = 1, ncol = k)
  notice.nr <- matrix(NA, nrow = 1, ncol = k)
  notice.date <- matrix(NA, nrow = length(jackknife$coefficients[,1]), ncol = k)
  susp.nr <- matrix(NA, nrow = 1, ncol = k)
  susp.date <- matrix(NA, nrow = length(jackknife$coefficients[,1]), ncol = k)
  rownames(notice.date) <- rownames(jackknife$coefficients)
  rownames(susp.date) <- rownames(jackknife$coefficients)
  colnames(notice.date) <- header[1:k]
  colnames(susp.date) <- header[1:k]
  crit1 <- rep(0.75, k)
  crit2 <- rep(0.5, k)

  for (z in (1:k)) {
    base.para[,z] <- base.model$coef[z]
    base.se[,z] <- sqrt(diag(vcov(base.model)))[z]
    base.p[,z] <- dt(base.para[,z]/base.se[,z],
                     df = horizon -k)
    base.lower[,z] <- base.para[,z] -
      qt(0.975, df = df.residual(base.model))*base.se[,z]
    base.upper[,z] <- base.para[,z] +
      qt(0.975, df = df.residual(base.model))*base.se[,z]
    tdrift.mean[,z] <- mean(tdrift[,z])
    tdrift.median[,z] <- median(tdrift[,z])
    tdrift.sd[,z] <- sd(tdrift[,z])
    tdrift.lower[,z] <- tdrift.mean[,z] -
      qnorm(0.975)*tdrift.sd[,z]
    tdrift.upper[,z] <- tdrift.mean[,z] +
      qnorm(0.975)*tdrift.sd[,z]
    max.beta[, z] <- as.numeric(max(hdrift[,z]))
    min.beta[, z] <- as.numeric(min(hdrift[,z]))
    pos1[, z] <- names(which.max(hdrift[,z]))
    pos2[, z] <- names(which.min(hdrift[,z]))
    hor.mean[, z] <- as.numeric(mean(hdrift[,z]))
    hor.sd[, z] <- as.numeric(sd(hdrift[,z]))

    notice.date[,z] <- 2*pt(abs(jackknife$coef[,z] / base.se[,z]),
                                     df = horizon-k, lower.tail = FALSE)
    susp.date[,z] <- 2*pt(abs(jackknife$coef[,z] / base.se[,z]),
                                   df = horizon-k, lower.tail = FALSE)
    notice.nr[,z] <- sum(as.numeric(notice.date[,z] <0.75))
    susp.nr[,z] <- sum(as.numeric(susp.date[,z] <0.5))
  }
  notice.date <- as.data.frame(notice.date)
  susp.date <- as.data.frame(susp.date)
  colnames(notice.date) <- header[1:k]
  colnames(susp.date) <- header[1:k]
  notice.date <- notice.date[rowSums(mapply(`<`, notice.date, crit1)) > 0, ]
  susp.date <- susp.date[rowSums(mapply(`<`, susp.date, crit2)) > 0, ]

  
  starify <- function (y) {if (y < .001) {
    y <- "***"
  } else if (y < .01) {
    y <- "**"
  } else if (y < .05) {
    y <-  "*"
  } else if (y < .1) {
    y <- "."
  } else {
    y <- ""
  }}
  sign <- lapply(base.p, starify)

  sd.from <- (base.para-tdrift.mean)/tdrift.sd
  stats <- data.frame(rbind(base.para, sign, base.se, base.p, base.upper,
                            base.lower, sd.from, tdrift.mean, tdrift.median,
                            tdrift.sd, tdrift.upper, tdrift.lower, max.beta,
                            pos1, min.beta, pos2, hor.mean, notice.nr, susp.nr),
                      row.names = c("base model coef", "significance", "std. error",
                                    "p-value", "base coef 95CI upper", "base coef 95CI lower",
                                    "sd from historic", "historic mean coef",
                                    "historic median coef", "historic sd",
                                    "historic mean 95CI upper","historic mean 95CI lower",
                                    "max. parameter", "max. @ horizon", "min. parameter",
                                    "min. @ horizon", "horizon mean coef",
                                    "noteworthy obs.(pval < 0.75)",
                                    "suspicious obs. (pval < 0.5)"))
  colnames(stats) <- header[1:k]
  sumstats <- list(stats, notice.date, susp.date)

  ###################################################
  ####                Outputs                    ####
  ###################################################
  results <- list(CALL = match.call(),
                  base.model = base.model,
                  tdrift = tdrift,
                  tdrift.se = tdrift.se,
                  hdrift = hdrift,
                  hdrift.se = hdrift.se,
                  jackknife = jackknife,
                  sumstats = sumstats)
  class(results) <- c("BDA", "list")

  ###################################################
  ####                Plotting                   ####
  ###################################################
  if (doplot == TRUE) {
    opar <- par(no.readonly = TRUE)
    for (j in (1:k))
    {
      par(mfrow = c(1,3), oma=c(0,0,2,0))
      ###################################################
      ####              Drift Plot                   ####
      ###################################################
      plot(as.numeric(tdrift[,j]), type = "l",  xaxt = "n",
           ylab = "estimated parameter", xlab="estimation date",
           main="time drift", ylim = c(min(tdrift[,j]), max(tdrift[,j])))
      axis(1, at=c(1, length(zoo::index(as.numeric(tdrift[,j])))),
           labels=c(start(tdrift[,j]), end(tdrift[,j])))
      polygon(c(zoo::index(as.numeric(tdrift[,j])),
                rev(zoo::index(as.numeric(tdrift[,j])))),
              c(c(rep((coef[j]+qt(0.975, df = df.residual(base.model))*se[j]),
                      length(zoo::index(as.numeric(tdrift[,j]))))),
                c(rep((coef[j]-qt(0.975, df = df.residual(base.model))*se[j]),
                      length(zoo::index(as.numeric(tdrift[,j])))))),
              col = scales::alpha("red", 0.15), border = FALSE)

      abline(h=c(0,1))
      abline(h = coef[j], lwd = 2, col = scales::alpha("red", 0.7))
      abline(h = mean(tdrift[,j]), lwd = 2, col = scales::alpha("blue", 0.7))
      sp1 <- smooth.spline(tdrift[,j], nknots = 15)
      lines(sp1, lty = 2, col = scales::alpha("blue", 0.5), lwd = 3)

      ###################################################
      ####            Horizon Plot                   ####
      ###################################################
      plot(hdrift[,j], type = "l", xaxt="n", ylab = "estimated parameter",
           main="horizon drift", xlab="estimation window size ",
           ylim = c(min(hdrift[,j]-qt(0.975, df = (min.hor:max.hor - length(coef(base.model))))*hdrift.se[,j]),
                    max(hdrift[,j]+qt(0.975, df = (min.hor:max.hor - length(coef(base.model))))*hdrift.se[,j])))
      axis(1, at=seq(from = 1,
                     to = max.hor-min.hor,
                     by = round((max.hor-min.hor)/10)),
           labels=seq(from = min.hor, to = max.hor,
                      by = round((max.hor-min.hor)/10)), las = 2)
      abline(h=c(0,1))
      polygon(c(zoo::index(hdrift), rev(zoo::index(hdrift))),
              c((hdrift[,j]+qt(0.975, df = min.hor:max.hor)*hdrift.se[,j]),
                rev(hdrift[,j]-qt(0.975, df = min.hor:max.hor)*hdrift.se[,j])),
              col = scales::alpha("blue", 0.15), border = FALSE)
      abline(v=(horizon-min.hor), col = scales::alpha("red", 0.5), lwd = 2)
      sp2 <- smooth.spline(hdrift[,j], nknots =5)
      lines(sp2, lty = 2, col = scales::alpha("blue", 0.5), lwd = 3)

      ###################################################
      ####           Jackknife Plot                  ####
      ###################################################
      plot(jackknife$coef[,j], type = "h",
           xlab = "index of omitted observations", ylim = c(-se[j],se[j]),
           ylab = "estimated deviance", main = "jackknife")
      polygon(c(1:length(jackknife$coef[,j]), length(jackknife$coef[,j]):1),
              c(rep((+qt(0.75, df = horizon-k)*se[j]),
                    length(zoo::index(jackknife$coef[,j]))),
                rep((-qt(0.75, df = horizon-k)*se[j]),
                    length(zoo::index(jackknife$coef[,j])))),
              col = scales::alpha("red", 0.1), border = FALSE)
      polygon(c(1:length(jackknife$coef[,j]), length(jackknife$coef[,j]):1),
              c(rep((+qt(0.625, df = horizon-k)*se[j]),
                    length(zoo::index(jackknife$coef[,j]))),
                rep((-qt(0.625, df = horizon-k)*se[j]),
                    length(zoo::index(jackknife$coef[,j])))),
              col = scales::alpha("red", 0.1), border = FALSE)
      abline(h=0)
      title(paste("Estimation for Parameter",toString(header[j])), outer=TRUE)
    }
    par(opar)
  }
  else {}
  message("BDA was completed successfully!")
  return(results)
}