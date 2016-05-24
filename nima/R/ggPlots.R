utils::globalVariables(c("ord.x", "z", "lower", "upper", "label"))

#' Quantile-Quantile Plots with ggplot2
#'
#' Produce standard quantile-quantile plots for modeling using ggplot2.
#'
#' @param x A numeric vector of residuals from a generalized linear model.
#' @param distribution The reference probability distribution for residuals.
#' @param ... Any additional parameters to be passed to distribution functions.
#' @param conf The confidence level to be used with confidence intervals.
#' @param labels The names to be used when identifying points on the Q-Q plot.
#' @param line.estimate Should quantiles be estimated, if so which quantiles?
#'
#' @importFrom stats quantile ppoints na.omit qnorm
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab ggtitle theme
#' @importFrom ggplot2 geom_abline geom_text geom_ribbon scale_size_continuous
#'
#' @export qqPlot_gg
#'
#' @examples
#' n <- 100; x1 <- rnorm(n); y1 <- rnorm(n);
#' linmod <- lm(y1 ~ x1)
#' x <- linmod$residuals
#' qqPlot_gg(x)

qqPlot_gg <- function(x, distribution = "norm", ..., line.estimate = NULL,
                      conf = 0.95, labels = names(x)) {
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))
  x <- na.omit(x)
  ord <- order(x)
  n <- length(x)
  P <- ppoints(length(x))
  df <- data.frame(ord.x = x[ord], z = q.function(P, ...))

  if(is.null(line.estimate)) {
    Q.x <- quantile(df$ord.x, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75), ...)
    b <- (diff(Q.x) / diff(Q.z))
    coef <- c(Q.x[1] - b * Q.z[1], b)
  } else {
    coef <- coef(line.estimate(ord.x ~ z))
  }

  zz <- qnorm(1 - ( (1 - conf) / 2))
  SE <- (coef[2] / d.function(df$z)) * sqrt(P * ( (1 - P) / n))
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE

  if(!is.null(labels)) {
    df$label <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower,labels[ord],"")
  }

  p <- ggplot(df, aes(x = z, y = ord.x)) +
    geom_point() +
    geom_abline(intercept = coef[1], slope = coef[2]) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    ylab("Model Residual Quantiles") + xlab("Theoretical Quantiles") +
    ggtitle("Quantile Plot Residuals")
  if(!is.null(labels)) {
    p <- p + geom_text( aes(label = label))
  }
  return(p)
}


#####################################
########### NEXT FUNCTION ###########
#####################################

utils::globalVariables(c("aes", ".fitted", ".resid", "geom_point", ".stdresid",
                       "geom_hline", "xlab", "ylab", "ggtitle", "stat_smooth",
                       "geom_abline", ".cooksd", "geom_bar", ".hat", "theme",
                       "scale_size_continuous"))

#' Linear Model Diagnostic Plots with ggplot2
#'
#' Produce standard diagnostic plots for linear models using ggplot2.
#'
#' @param model A linear model object produced by \code{lm()}.
#'
#' @importFrom ggplot2 ggplot aes geom_point stat_smooth geom_hline geom_abline
#' @importFrom ggplot2 xlab ylab ggtitle geom_bar scale_size_continuous theme
#' @importFrom gridExtra grid.arrange
#'
#' @export lmPlots_gg
#'
#' @examples
#' n <- 100; x1 <- rnorm(n); y1 <- rnorm(n);
#' linmod <- lm(y1 ~ x1)
#' lmPlots_gg(linmod)

lmPlots_gg <- function(model) {
  p1 <- ggplot(model, aes(x = .fitted, y = .resid)) + geom_point() +
        stat_smooth(method = "loess") +
        geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
        xlab("Fitted vals") + ylab("Residuals") + ggtitle("Residual vs. Fitted")

  p2 <- ggplot(model, aes(sample = .stdresid)) + geom_point(stat = "qq") +
        geom_abline() + xlab("Theoretical Quantiles") + ylab("Std. Residuals") +
        ggtitle("Gausian Q-Q")

  p3 <- ggplot(model, aes(x = .fitted, y = sqrt(abs(.stdresid)))) +
        geom_point(na.rm = TRUE) + stat_smooth(method = "loess", na.rm = TRUE) +
        xlab("Fitted Value") + ylab(expression(sqrt("|Std. residuals|"))) +
        ggtitle("Scale-Location")

  p4 <- ggplot(model, aes(x = seq_along(.cooksd), y = .cooksd)) +
        geom_bar(stat = "identity", position = "identity") +
        xlab("Obs. No.") + ylab("Cook's distance") + ggtitle("Cook's distance")

  p5 <- ggplot(model, aes(x = .hat, y = .stdresid)) +
        geom_point(aes(size = .cooksd), na.rm = TRUE) +
        stat_smooth(method = "loess", na.rm = TRUE) + xlab("Leverage") +
        ggtitle("Resid. vs Leverage") + ylab("Std. Resid.") +
        scale_size_continuous("Cook's Distance", range = c(1,5)) +
        theme(legend.position = "none")

  p6 <- ggplot(model, aes(x = .hat, y = .cooksd)) + geom_point(na.rm = TRUE) +
        stat_smooth(method = "loess", na.rm = TRUE) + xlab("Leverage") +
        ylab("Cook's Distance") + ggtitle("Cook's dist. vs. Lev.") +
        geom_abline(slope = seq(0,3,0.5), color = "gray", linetype = "dashed")

  return(gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3))
}


#####################################
########### NEXT FUNCTION ###########
#####################################

utils::globalVariables(c("time", "cens", "surv", "low", "up", "group"))

#' Kaplan-Meier Plot with ggplot2
#'
#' Produce a survival plot of a Kaplan-Meier estimator using ggplot2.
#'
#' @param s Survival model object, generated by methods like Kaplan-Meier.
#' @param CI Type of confidence interval for the survival object.
#' @param survcl Color for observations with outcomes of "survived".
#' @param censcl Color for observations with outcomes of "censored".
#' @param bw Boolean for desired background color in plots (black/white).
#' @param xlab A label for the x-axis, defaults to "Time".
#' @param ylab A label for the y-axis, defaults to "Survival".
#' @param main A main label for the survival plot, no default.
#' @param pltcens Should the plot include the censored values?
#' @param lty The number of survival curves to be generated from the model. 
#' @param ltci The number of lines to be generated for confidence intervals. 
#' @param shape The shapes of points plotted, passed to \code{geom_point}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline xlab ylab ggtitle
#' @importFrom ggplot2 stat_smooth geom_abline geom_bar geom_step theme 
#' @importFrom ggplot2 scale_size_continuous scale_colour_discrete theme_bw 
#' @importFrom ggplot2 scale_linetype_manual scale_colour_manual
#' @importFrom survival survfit Surv
#'
#' @export survPlot_gg
#'
#' @examples
#' library(survival)
#' time <- c(5,6,8,3,22)
#' age <- c(46,35,30,30,36)
#' drug <- c(0,1,1,1,0)
#' censor <- c(1,0,1,1,1)
#' survdat <- as.data.frame(cbind(time, age, drug, censor))
#' s <- survfit(Surv(survdat$time, survdat$censor) ~ 1, conf.type="none")
#' survPlot_gg(s)

survPlot_gg <- function(s, CI = "def", pltcens = TRUE, survcl = "gg.def",
                        censcl = "red", lty = 1, ltci = 2, shape = 3, 
                        bw = FALSE, xlab = "Time", ylab = "Survival",
                        main = "") {

  strata <- ifelse(is.null(s$strata) == TRUE, 1, length(s$strata))
  stopifnot(length(survcl) == 1 | length(survcl) == strata)
  stopifnot(length(lty) == 1 | length(lty) == strata)

  ggsurv.s <- function(s, CI = "def", pltcens = TRUE, survcl = "gg.def",
                       censcl = "red", lty = 1, ltci = 2, shape = 3, bw = FALSE, 
                       xlab = "Time", ylab = "Survival", main = ""){

    dat <- data.frame(time = c(0, s$time), surv = c(1, s$surv),
                      up = c(1, s$upper), low = c(1, s$lower),
                      cens = c(0, s$n.censor))
    dat.cens <- subset(dat, cens != 0)

    col <- ifelse(survcl == "gg.def", "black", survcl)

    pl <- ggplot(dat, aes(x = time, y = surv)) + xlab(xlab) + ylab(ylab) +
          ggtitle(main) + geom_step(col = col, lty = lty)

    pl <- if(CI == TRUE | CI == "def") {
      pl + geom_step(aes(y = up), color = col, lty = ltci) +
        geom_step(aes(y = low), color = col, lty = ltci)
    } else (pl)

    pl <- if(pltcens == TRUE & length(dat.cens) > 0) {
      pl + geom_point(data = dat.cens, aes(y = surv), shape = shape,
                      col = censcl)
    } else if (pltcens == TRUE & length(dat.cens) == 0) {
      stop ("There are no censored observations")
    } else (pl)

    pl <- if(bw == TRUE) {
      pl + theme_bw()
    } else (pl)
    pl
  }

  ggsurv.m <- function(s, CI = "def", pltcens = TRUE, survcl = "gg.def",
                       censcl = "red", lty = 1, ltci = 2, shape = 3, bw = FALSE, 
                       xlab = "Time", ylab = "Survival Prob.", main = "") {
    n <- s$strata
    groups <- factor(unlist(strsplit(names(s$strata), "="))[seq(2, 2 * strata,
                                                                by = 2)])
    gr.name <-  unlist(strsplit(names(s$strata), "="))[1]
    gr.df <- vector("list", strata)
    ind <- vector("list", strata)
    n.ind <- c(0,n); n.ind <- cumsum(n.ind)
    for (i in 1:strata) ind[[i]] <- (n.ind[i] + 1):n.ind[i + 1]

    for (i in 1:strata) {
      gr.df[[i]] <- data.frame(
        time = c(0, s$time[ ind[[i]] ]),
        surv = c(1, s$surv[ ind[[i]] ]),
        up = c(1, s$upper[ ind[[i]] ]),
        low = c(1, s$lower[ ind[[i]] ]),
        cens = c(0, s$n.censor[ ind[[i]] ]),
        group = rep(groups[i], n[i] + 1))
    }

    dat <- do.call(rbind, gr.df)
    dat.cens <- subset(dat, cens != 0)

    pl <- ggplot(dat, aes(x = time, y = surv, group = group)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(aes(col = group, lty = group))

    col <- if(length(survcl == 1)) {
      scale_colour_manual(name = gr.name, values = rep(survcl, strata))
    } else{
        scale_colour_manual(name = gr.name, values = survcl)
    }

    pl <- if(survcl[1] != "gg.def") {
      pl + col
    } else {
        pl + scale_colour_discrete(name = gr.name)
    }

    line <- if(length(lty) == 1) {
      scale_linetype_manual(name = gr.name, values = rep(lty, strata))
    } else {
        scale_linetype_manual(name = gr.name, values = lty)
    }

    pl <- pl + line

    pl <- if(CI == TRUE) {
      if(length(survcl) > 1 && length(lty) > 1){
        stop("Either survcl or lty should be of length 1 in order to plot
             95% CI with multiple strata")
      } else if((length(survcl) > 1 | survcl == "gg.def")[1]) {
        pl + geom_step(aes(y = up, color = group), lty = ltci) +
          geom_step(aes(y = low, color = group), lty = ltci)
      } else {
          pl +  geom_step(aes(y = up, lty = group), col = survcl) +
          geom_step(aes(y = low,lty = group), col = survcl)
        }
    } else (pl)

    pl <- if(pltcens == TRUE & length(dat.cens) > 0) {
      pl + geom_point(data = dat.cens, aes(y = surv), shape = shape,
                      col = censcl)
    } else if (pltcens == TRUE & length(dat.cens) == 0) {
      stop ("There are no censored observations")
    } else (pl)

    pl <- if(bw == TRUE) {pl + theme_bw()
    } else (pl)
    pl
  }
  pl <- if(strata == 1) {
          ggsurv.s(s, CI , pltcens, survcl , censcl, lty, ltci,
                   shape, bw, xlab, ylab, main)
  } else {
      ggsurv.m(s, CI, pltcens, survcl , censcl, lty, ltci,
               shape, bw, xlab, ylab, main)
         }
  return(pl)
}
