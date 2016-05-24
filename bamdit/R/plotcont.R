#' plotcont
#'
#' This function plots the observe data in the ROC (Receiving Operating Charachteristics) space with the
#' posterior predictive contours. The predictive curves are approximated using a non-parametric smoother or with a parametric model. For the parametric model the current implementation supports only a logistic link function.
#'
#' @param m The modelfit.
#' @param data The data frame used to fit the model.
#' @param parametric.smooth Indicates if the predictive curve is a parametric or non-parametric smoother
#' @param level credibility levels of the predictive curve
#' @param link The link function used to fit the model. Possible values are \emph{logit}, \emph{cloglog} \emph{probit}.
#' @param title Optional parameter for setting a title in the plot.
#' @seealso \code{\link{metadiag}}.
#' @keywords file
#' @examples
#'
#' ## execute analysis
#' \dontrun{
#' data(glas)
#' glas.t <- glas[glas$marker == "Telomerase", 1:4]
#' glas.m1 <- metadiag(glas.t,
#' re = "normal",
#' link = "logit",
#' nr.burnin = 1000,
#' nr.iterations = 10000,
#' nr.chains = 4,
#' r2jags = TRUE)
#'
#' plotcont(m = glas.m1, data = glas.t)
#' }
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import R2jags
#' @import rjags
#' @importFrom stats cor integrate median pnorm qchisq quantile sd

#' @export
plotcont <- function(m, data,
                     parametric.smooth = TRUE,
                     level = c(0.5, 0.75, 0.95),
                     link = "logit",
                     title = paste("Predictive Posterior Contours (50%, 75% and 95%)")
                      )
{

  if(class(m)!="rjags")stop("You have to provide a valid rjags object as fitted model.")

  link.test <- link %in% c("logit", "cloglog", "probit")
  if(!link.test)stop("This link function is not implemented.")

  for(i in 1:length(level))
    {if(!(0 < level[i] & level[i] < 1))
      stop("At least one contour level is not correct. Use values between 0 and 1.")
     }

  # Predictive curve for model m1 ............................................................
  # Patch for the "note" no-visible-binding-for-global-variable

  dat.pred <- data.frame(fpr.new = 1 - m$BUGSoutput$sims.list$sp.new,
                          tpr.new = m$BUGSoutput$sims.list$se.new)
  tp <- data[,1]
  n1 <- data[,2]
  fp <- data[,3]
  n2 <- data[,4]

  if(tp>n1 || fp>n2)stop("the data is inconsistent")

  if(!missing(data)){
    tpr <-  tp / n1
    fpr <-  fp / n2
      n <- n1 + n2
  }else
    stop("NAs are not alow in this plot function")

  dat.hat <- data.frame(tpr = tpr,
                        fpr = fpr,
                        n = n)

# Base plot ...............................................................................
baseplot  <- ggplot(dat.pred, aes_string(x = "fpr.new", y = "tpr.new"))+
  scale_x_continuous(name = "FPR (1 - Specificity)", limits = c(0, 1)) +
  scale_y_continuous(name = "TPR (Sensitivity)", limits = c(0, 1)) +
  geom_point(data = dat.hat, aes_string(x = "fpr", y = "tpr", size = "n"),
             shape = 21, alpha = 0.35, fill = "blue") +
  scale_size_area(max_size = 10) +
  ggtitle(title)

# Non-parametric ..........................................................................
  nparplot  <- baseplot + stat_density2d(aes_string(colour = "..level.."))

# Parametric ..............................................................................
 if(link == "logit"){
           x <- with(dat.pred, log(fpr.new/(1-fpr.new)))
           y <- with(dat.pred, log(tpr.new/(1-tpr.new)))
 } else if(link == "cloglog")
 {
   x <- with(dat.pred, log(-log(1 - fpr.new)))
   y <- with(dat.pred, log(-log(1 - tpr.new)))
 } else if(link == "probit")
 {
   x <- with(dat.pred, qnorm(fpr.new))
   y <- with(dat.pred, qnorm(tpr.new))
 }

  mean.invlink.x <- mean(x)
  mean.invlink.y <- mean(y)
    sd.invlink.x <- sd(x)
    sd.invlink.y <- sd(y)
 rho.invlink.x.y <- cor(x,y)

parplot <- baseplot
for(k in level)
{
  # Set up the curve 95%...
  cc <- sqrt(qchisq(k, 2))
  tt <- seq(0, 2*pi, len=100)

  mu.x <-  mean(mean.invlink.x) + mean(sd.invlink.x) * cc * cos(tt)
  mu.y <-  mean(mean.invlink.y) + mean(sd.invlink.y) * cc * cos( tt + acos(mean(rho.invlink.x.y)))

  # Pooled summaries
  if(link == "logit"){
    TPR1 <- 1 / ( 1 + exp(-1*mu.y) )  # with logit link
    FPR1 <- 1 / ( 1 + exp(-1*mu.x) )  # with logit link
  } else if(link == "cloglog")
  {
    TPR1 <- 1 - exp(-1*exp(mu.y))      # with cloglog link
    FPR1 <- 1 - exp(-1*exp(mu.x))      # with cloglog link

  } else if(link == "probit")
  {
    TPR1 <- pnorm(mu.y)               # with logit link
    FPR1 <- pnorm(mu.x)               # with logit link
  }

  curve.alpha <- data.frame(FPR1, TPR1)
      parplot <- parplot + geom_path(data = curve.alpha, aes_string(x = "FPR1", y = "TPR1"))
}


if(parametric.smooth == TRUE) finalplot <- parplot
  else finalplot <- nparplot

return(finalplot)
}

