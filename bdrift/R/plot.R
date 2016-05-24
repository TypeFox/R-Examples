#' @title Plot Beta Drift Analyses
#'
#' @description
#' \code{plot.BDA} plots results from beta drift analyses.
#'
#' @details
#' \code{plot.BDA} produces three plots for each parameter of the baseline
#' model of the corresponding \code{BDA} function. Unless \code{single} is
#' set to \code{TRUE}, all three plots for a parameter are displayed in a
#' single plot window. 
#'
#' The first plot, titled "time drift", displays the drift of the parameter
#' across time. In addition to the parameter itself, a 15-knot cubic smooth spline
#' is displayed as a light-blue dashed line. The solid horizontal red line
#' represents the parameter estimate of the baseline model with the 95%
#' confidence interval as a red-shaded area surrounding it. The solid
#' horizontal blue line represents the mean of the drift series.
#'
#' The second plot, titled "horizon drift", displays the drift of the parameter
#' with respect to the estimation window size. In addition to the parameter
#' itself, a 5-knot cubic smooth spline is displayed as a light-blue dashed
#' line. The blue-shaded area represents the 95% confidence interval for the
#' parameters. The solid red vertical line highlights the estimation window size
#' of the baseline model.
#'
#' The third plot, titled "jackknife", displays the outcome of the jackknife
#' procedure for the baseline model. The light red-shaded and dark
#' red-shaded areas represent p-values <0.5 and <0.75 respectively,
#' as implied by the baseline model.
#'
#' \code{IMPORTANT NOTE}: This package was developed with the GUI
#' of RStudio in mind. The plotting function creates a potentially large number
#' of plots which can be comfortably viewed in RStudio, but require some
#' preparations in the standard R GUI. Start by executing \code{dev.new()}, which
#' opens a graphical device. Next, click the "History" tab and then click
#' "Recording" in the drop-down menu. If you run the plotting function now,
#' you can jump through the plots using the PageUp and PageDown key on your
#' keyboard.
#'
#' @param x an object of class \code{BDA}.
#' @param single logical. If \code{TRUE}, grouping of plots by parameters is disabled.
#' @param ... additional parameters.
#' @method plot BDA
#' @export
#' @return NULL
#' @author Markus Peter Auer <mp.auer@@meanerreversion.com>
#' @examples
#' \dontrun{
#' ###################################################
#' ####             Full example                  ####
#' ###################################################
#' 
#' results <- BDA(data = FFfactors, spec = (VOO~SP500),
#'                horizon = 250, doplot = FALSE)
#' plot(results)
#' }
#' 
#' ###################################################
#' ####        CRAN-compatible example            ####
#' ###################################################
#' 
#' results <- BDA(data = FFfactors[nrow(FFfactors):(nrow(FFfactors)-300),], 
#'                spec = (VOO~SP500),horizon = 250, doplot = FALSE)
#' plot(results)
#' message("NOTE: This is a shortened example. Reference the manual for more complex examples")

plot.BDA <- function(x, single = FALSE, ...){
  if (!inherits(x, "BDA"))
    stop("Object must be of class 'BDA'")
  horizon <- length(x$base.model$fitted.values)
  min.hor <- as.numeric(rownames(x$hdrift)[1])
  max.hor <- as.numeric(rownames(x$hdrift)[length(x$hdrift[,1])])
  k <- length(x$base.model$coefficients)
  coefs <- x$base.model$coef[1:k]
  se <- sqrt(diag(vcov(x$base.model)))[1:k]
  header <- c(names(coefs[1:k]))
  opar <- par(no.readonly = TRUE)
  ifelse(single == FALSE,
         par(mfrow = c(1,3), oma=c(0,0,2,0)),
         par(mfrow = c(1,1), oma=c(0,0,2,0)))
  for (j in 1:k)
  {
    ###################################################
    ####              Drift Plot                   ####
    ###################################################
    plot(as.numeric(x$tdrift[,j]), type = "l",  xaxt = "n",
         ylab = "estimated parameter", xlab="estimation date",
         main="time drift", ylim = c(min(x$tdrift[,j]), max(x$tdrift[,j])),...)
    axis(1, at=c(1, length(zoo::index(as.numeric(x$tdrift[,j])))),
         labels=c(start(x$tdrift[,j]), end(x$tdrift[,j])))
    polygon(c(zoo::index(as.numeric(x$tdrift[,j])),
              rev(zoo::index(as.numeric(x$tdrift[,j])))),
            c(c(rep((coefs[j]+qt(0.975, df = df.residual(x$base.model))*se[j]),
                    length(zoo::index(as.numeric(x$tdrift[,j]))))),
              c(rep((coefs[j]-qt(0.975, df = df.residual(x$base.model))*se[j]),
                    length(zoo::index(as.numeric(x$tdrift[,j])))))),
            col = scales::alpha("red", 0.15), border = FALSE)

    abline(h=c(0,1))
    abline(h = coefs[j], lwd = 2, col = scales::alpha("red", 0.7))
    abline(h = mean(x$tdrift[,j]), lwd = 2, col = scales::alpha("blue", 0.7))
    sp1 <- smooth.spline(x$tdrift[,j], nknots = 15)
    lines(sp1, lty = 2, col = scales::alpha("blue", 0.5), lwd = 3)
    if (single == TRUE){
      title(paste("Estimation for Parameter",toString(header[j])), outer=TRUE)}

    ###################################################
    ####            Horizon Plot                   ####
    ###################################################
    plot(x$hdrift[,j], type = "l", xaxt="n", ylab = "estimated parameter",
         main="horizon drift", xlab="estimation window size ",
         ylim = c(min(x$hdrift[,j]-qt(0.975, df = (min.hor:max.hor - length(coef(x$base.model))))*x$hdrift.se[,j]),
                  max(x$hdrift[,j]+qt(0.975, df = (min.hor:max.hor - length(coef(x$base.model))))*x$hdrift.se[,j])), 
         ...)
    axis(1, at=seq(from = 1,
                   to = max.hor-min.hor,
                   by = round((max.hor-min.hor)/10)),
         labels=seq(from = min.hor, to = max.hor,
                    by = round((max.hor-min.hor)/10)), las = 2)
    abline(h=c(0,1))
    polygon(c(zoo::index(x$hdrift), rev(zoo::index(x$hdrift))),
            c((x$hdrift[,j]+qt(0.975, df = min.hor:max.hor)*x$hdrift.se[,j]),
              rev(x$hdrift[,j]-qt(0.975, df = min.hor:max.hor)*x$hdrift.se[,j])),
            col = scales::alpha("blue", 0.15), border = FALSE)
    abline(v=(horizon-min.hor), col = scales::alpha("red", 0.5), lwd = 2)
    sp2 <- smooth.spline(x$hdrift[,j], nknots =5)
    lines(sp2, lty = 2, col = scales::alpha("blue", 0.5), lwd = 3)
    if (single == TRUE){
      title(paste("Estimation for Parameter",toString(header[j])), outer=TRUE)}

    ###################################################
    ####           Jackknife Plot                  ####
    ###################################################
    plot(x$jackknife$coef[,j], type = "h",
         xlab = "index of omitted observations", ylim = c(-se[j],se[j]),
         ylab = "estimated deviance", main = "jackknife", ...)
    polygon(c(1:length(x$jackknife$coef[,j]), length(x$jackknife$coef[,j]):1),
            c(rep((+qt(0.75, df = horizon-k)*se[j]),
                  length(zoo::index(x$jackknife$coef[,j]))),
              rep((-qt(0.75, df = horizon-k)*se[j]),
                  length(zoo::index(x$jackknife$coef[,j])))),
            col = scales::alpha("red", 0.1), border = FALSE)
    polygon(c(1:length(x$jackknife$coef[,j]), length(x$jackknife$coef[,j]):1),
            c(rep((+qt(0.625, df = horizon-k)*se[j]),
                  length(zoo::index(x$jackknife$coef[,j]))),
              rep((-qt(0.625, df = horizon-k)*se[j]),
                  length(zoo::index(x$jackknife$coef[,j])))),
            col = scales::alpha("red", 0.1), border = FALSE)
    abline(h=0)
    title(paste("Estimation for Parameter",toString(header[j])), outer=TRUE)
  }
  par(opar)
}
