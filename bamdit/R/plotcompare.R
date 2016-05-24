#' plotcompare
#'
#' This function compares the predictive posterior surfaces of two fitted models. The current implementation supports only a logistic link function.
#'
#' @param m1 A model fitted to the data
#' @param m2 A second model fitted to the data
#' @param data The data frame used to fit the model.
#' @param level credibility level of the predictive curves
#' @param link1 The link function used to fit the model m1. Possible values are \emph{logit}, \emph{cloglog} \emph{probit}.
#' @param link2 The link function used to fit the model m2. Possible values are \emph{logit}, \emph{cloglog} \emph{probit}.
#' @param title The title of the plot
#' @param m1.name Label of the model 1
#' @param m2.name Label of the model 2
#' @param group A factor variable to display data of different groups
#' @param limits.x A vector with the limits of the horizontal axis
#' @param limits.y A vector with the limits of the vertical axis
#' @param group.colors A character vector with two color names
#' @seealso \code{\link{metadiag}}.
#' @keywords file
#' @examples
#'
#' ## execute analysis
#' \dontrun{
#' # Comparing results from two models same data
#' data(glas)
#' glas.t <- glas[glas$marker == "Telomerase", 1:4]
#' glas.m1 <- metadiag(glas.t)
#' glas.m2 <- metadiag(glas.t, re = "sm")
#' plotcompare(m1 = glas.m1, m2 = glas.m2)
#'
#' # Comparing results from two models fitted to two groups of data
#' # Studies with retrospective design
#' m1.ct <- metadiag(ct[ct$design==1, 1:4])
#'
#' # Studies with prospective design
#' m1.ct <- metadiag(ct[ct$design==1, 1:4])
#'
#' plotcompare(m1.ct, m2.ct, data = ct,
#'            m1.name = "Retrospective design", m2.name = "Prospective design",
#'              group = factor(ct$design),
#'              limits.x = c(0, 0.75), limits.y = c(0.65, 1))
#' }
#'
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import R2jags
#' @import rjags
#' @importFrom stats cor integrate median pnorm qchisq quantile sd

#'@export
plotcompare <- function(m1, m2, data, level = 0.95,
                        link1 = "logit",
                        link2 = "logit",
                        title = paste("Comparative Predictive Posterior Contours") ,
                        m1.name = "Model.1",
                        m2.name = "Model.2",
                        group = NULL,
                        limits.x = c(0, 1),
                        limits.y = c(0, 1),
                        group.colors = c("blue", "red"))
{


  # Checks ..................................................................................
  link.test1 <- link1 %in% c("logit", "cloglog", "probit")
  if(!link.test1)stop("Problem in the link function of model 1. This link function is not implemented")

  link.test2 <- link2 %in% c("logit", "cloglog", "probit")
  if(!link.test2)stop("Problem in the link function of model 2. This link function is not implemented")

  if(length(level)>1)stop("Only one contour level is accepted to compare models in this function.")

  if(!(0 < level & level < 1))stop("At least one contour level is not correct. Use values between 0 and 1!")

  # Base plot ...............................................................................

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

  if(is.null(group)){
  # The  plot ...............................................................................
  baseplot  <- ggplot(dat.hat, aes_string(x = "fpr", y = "tpr"))+
                scale_x_continuous(name = "FPR (1 - Specificity)", limits = limits.x) +
                scale_y_continuous(name = "TPR (Sensitivity)", limits = limits.y) +
                geom_point(data=dat.hat, aes(x=fpr, y=tpr, size=n), shape = 21, alpha = 0.35, fill="blue") +
                scale_size_area(max_size = 10) +
                ggtitle(title)
  # .........................................................................................
  } else if (is.factor(group)) {

    dat.hat <- data.frame(tpr = tpr,
                           fpr = fpr,
                           n = n,
                           group = group)

    baseplot  <- ggplot(dat.hat, aes_string(x = "fpr", y = "tpr"))+
      scale_x_continuous(name = "FPR (1 - Specificity)", limits = limits.x) +
      scale_y_continuous(name = "TPR (Sensitivity)", limits = limits.y) +
      geom_point(data=dat.hat, aes_string(x = "fpr", y = "tpr", size = "n", fill = "group"),
                 shape = 21, alpha = 0.35) +
      scale_size_area(max_size = 10) +
      scale_fill_manual(values = group.colors) +
      ggtitle(title)

  }
 # Predictive curve for model m1 ............................................................
 # Patch for the "note" no-visible-binding-for-global-variable

 dat.pred1 <- data.frame(fpr.new = 1 - m1$BUGSoutput$sims.list$sp.new,
                         tpr.new = m1$BUGSoutput$sims.list$se.new)
 # Set up the curve
 cc <- sqrt(qchisq(level, 2))
 tt <- seq(0, 2*pi, len=100)

 if(link1 == "logit"){
   x1 <- with(dat.pred1, log(fpr.new/(1-fpr.new)))
   y1 <- with(dat.pred1, log(tpr.new/(1-tpr.new)))
 } else if(link1 == "cloglog")
 {
   x1 <- with(dat.pred1, log(-log(1 - fpr.new)))
   y1 <- with(dat.pred1, log(-log(1 - tpr.new)))
 } else if(link1 == "probit")
 {
   x1 <- with(dat.pred1, qnorm(fpr.new))
   y1 <- with(dat.pred1, qnorm(tpr.new))
 }

 mean.invlink.x1 <- mean(x1)
 mean.invlink.y1 <- mean(y1)
 sd.invlink.x1 <- sd(x1)
 sd.invlink.y1 <- sd(y1)
 rho.invlink.x.y1 <- cor(x1, y1)

 # Set up the curve ...
 mu.x1 <-  mean(mean.invlink.x1) + mean(sd.invlink.x1) * cc * cos(tt)
 mu.y1 <-  mean(mean.invlink.y1) + mean(sd.invlink.y1) * cc * cos( tt + acos(mean(rho.invlink.x.y1)))

 # Pooled summaries
 if(link1 == "logit"){
   TPR1 <- 1 / ( 1 + exp(-1*mu.y1) )   # with logit link
   FPR1 <- 1 / ( 1 + exp(-1*mu.x1) )   # with logit link
 } else if(link1 == "cloglog")
 {
   TPR1 <- 1 - exp(-1*exp(mu.y1))      # with cloglog link
   FPR1 <- 1 - exp(-1*exp(mu.x1))      # with cloglog link

 } else if(link1 == "probit")
 {
   TPR1 <- pnorm(mu.y1)                # with logit link
   FPR1 <- pnorm(mu.x1)                # with logit link
 }

 # Predictive curve for model m2 ............................................................

 dat.pred2 <- data.frame(fpr.new = 1 - m2$BUGSoutput$sims.list$sp.new,
                         tpr.new = m2$BUGSoutput$sims.list$se.new)


 if(link2 == "logit"){
   x2 <- with(dat.pred2, log(fpr.new/(1-fpr.new)))
   y2 <- with(dat.pred2, log(tpr.new/(1-tpr.new)))
 } else if(link2 == "cloglog")
 {
   x2 <- with(dat.pred2, log(-log(1 - fpr.new)))
   y2 <- with(dat.pred2, log(-log(1 - tpr.new)))
 } else if(link2 == "probit")
 {
   x2 <- with(dat.pred2, qnorm(fpr.new))
   y2 <- with(dat.pred2, qnorm(tpr.new))
 }

 mean.invlink.x2 <- mean(x2)
 mean.invlink.y2 <- mean(y2)
 sd.invlink.x2 <- sd(x2)
 sd.invlink.y2 <- sd(y2)
 rho.invlink.x.y2 <- cor(x2, y2)

 # Set up the curve ...
 mu.x2 <-  mean(mean.invlink.x2) + mean(sd.invlink.x2) * cc * cos(tt)
 mu.y2 <-  mean(mean.invlink.y2) + mean(sd.invlink.y2) * cc * cos( tt + acos(mean(rho.invlink.x.y2)))

 # Pooled summaries
 if(link2 == "logit"){
   TPR2 <- 1 / ( 1 + exp(-1*mu.y2) )   # with logit link
   FPR2 <- 1 / ( 1 + exp(-1*mu.x2) )   # with logit link
 } else if(link2 == "cloglog")
 {
   TPR2 <- 1 - exp(-1*exp(mu.y2))      # with cloglog link
   FPR2 <- 1 - exp(-1*exp(mu.x2))      # with cloglog link

 } else if(link2 == "probit")
 {
   TPR2 <- pnorm(mu.y2)                # with logit link
   FPR2 <- pnorm(mu.x2)                # with logit link
 }

# data frame to plot
curve.alpha <- data.frame(FPR.pred = c(FPR1, FPR2),
                          TPR.pred = c(TPR1, TPR2),
                          model = gl(2, 100, labels = c(m1.name, m2.name))
                          )

# Final plot ..............................................................................
finalplot <- baseplot +
              geom_path(data = curve.alpha, aes_string(x = "FPR.pred", y = "TPR.pred", linetype = "model"))+
              theme(legend.position="top")

return(finalplot)
}



