#' bsroc
#'
#' This function plots the observe data in the ROC (Receiving Operating Charachteristics) space with the
#' Bayesian SROC (Summery ROC) curve. The predictive curves are approximated using a parametric model.
#'
#' @param m                    The model fitted.
#' @param data                 The data frame used to fit the model.
#' @param level                Credibility levels of the predictive curve
#' @param link                 The link function used to fit the model. Possible values are \emph{logit}, \emph{cloglog} \emph{probit}.
#' @param title                Optional parameter for setting a title in the plot.
#' @param fpr.x                Grid of values where the conditionlal distribution is calculated.
#' @param xlim.bsroc           Limits of the x-axis for the BSROC curve plot.
#' @param ylim.bsroc           Limits of the y-axis for the BSROC curve plot.
#' @param lower.auc            Lower limit of the AUC.
#' @param upper.auc            Upper limit of the AUC.
#' @param col.fill.points      Color used to fill points, default is blue.
#' @param results.bauc         Print results of the Bayesian Area Under the Curve, default value is TRUE.
#' @param results.bsroc        Print results of the Bayesian SROC curve, default value is FALSE.
#' @param plot.post.bauc       The BSROC and the posterior of the BAUC are ploted in the same page, default is FALSE.
#' @param binwidth.p           Histograms' binwidth, default is 1/30 range of the data.
#' @param scale.size.area      Scale area for the ploted points, default = 10.
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
#' bsroc(m = glas.m1, data = glas.t)
#' }
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import R2jags
#' @import rjags
#' @importFrom stats cor integrate median pnorm qchisq quantile sd
#' @export
bsroc <- function(m,
                  data,
                  level = c(0.05, 0.5, 0.95),
                  link = "logit",
                  title = "Bayesian SROC Curve",
                  fpr.x =  seq(0.01, 0.99, 0.01),
                  xlim.bsroc = c(0,1),
                  ylim.bsroc = c(0,1),
                  lower.auc = 0,
                  upper.auc = 0.99,
                  col.fill.points = "blue",
                  results.bauc = TRUE,
                  results.bsroc = FALSE,
                  plot.post.bauc = FALSE,
                  binwidth.p = 1/30,
                  scale.size.area = 10 )
{

  link.test <- link == "logit"      #%in%  c("logit", "cloglog", "probit")
  if(!link.test)stop("This link function is not implemented")

  for(i in 1:length(level))
  {if(!(0 < level[i] & level[i] < 1))
    stop("At least one credibility level is not correct. Use values between 0 and 1!")
  }

  tp <- data[,1]
  n1 <- data[,2]
  fp <- data[,3]
  n2 <- data[,4]

  if(tp>n1 || fp>n2)stop("the data is inconsistent")

  if(!missing(data)){
    tpr <- tp / n1
    fpr <- fp / n2
      n <- n1 + n2
  }else
    stop("NAs are not alow in this plot function")

  dat.hat <- data.frame(tpr = tpr,
                        fpr = fpr,
                          n = n)
  # Parametric ..............................................................................

  # Patch for the "note" no-visible-binding-for-global-variable
  A <- m$BUGSoutput$sims.list$mu.D
  rho <- m$BUGSoutput$sims.list$rho
  sigma.D <- m$BUGSoutput$sims.list$sigma.D
  sigma.S <- m$BUGSoutput$sims.list$sigma.S
  B <- rho*sigma.D/sigma.S

  ###############################################################################
  # BSROC
  f <- function(FPR, A, B){
    x <- A/(1-B) + (B+1)/(1-B)*log(FPR/(1-FPR))
    f.value <- exp(x)/(1+exp(x))
    return(f.value)
  }

  ###############################################################################
  # Simple graph for the BSROC

  BSROC  <- Vectorize(f, vectorize.args = c("A", "B"))
  BSROC.out <- BSROC(A, B, FPR = fpr.x)
  bsrocCI <- apply(BSROC.out, 1, quantile, prob = sort(level), na.rm =T)
  bsrocCI <- t(bsrocCI)

  dat.bsroc <- data.frame(x = fpr.x,
                          y.low = bsrocCI[ ,1],
                          y.med = bsrocCI[ ,2],
                           y.up = bsrocCI[ ,3])

  bsroc.plot <-  ggplot(data = dat.bsroc, aes_string(x = "x", y = "y.med")) +
                   geom_line(aes_string(x = "x", y = "y.up")) +
                   geom_line(aes_string(x = "x", y = "y.med")) +
                   geom_line(aes_string(x = "x", y = "y.low")) +
                   geom_point(data = dat.hat, aes_string(size = "n", x = "fpr", y = "tpr"),
                              shape = 21, alpha = 0.35, fill = col.fill.points)+
                   scale_x_continuous(name = "FPR (1 - Specificity)", limits = xlim.bsroc ) +
                   scale_y_continuous(name = "TPR (Sensitivity)", limits = ylim.bsroc) +
                   scale_size_area(max_size = scale.size.area) +
                   ggtitle(title)

  # Bayesian Area Under the Curve (BAUC)...........................................................

  # BAUC function
  area <- function(A, B)integrate(f, lower = lower.auc,
                                     upper = upper.auc,
                                         A = A,
                                         B = B)$value
  v.area <- Vectorize(area)

  # Calculations restricted for the slope -1 < B < 1 (Verde 2008 pag 11.)..........................
  #
  index.B.range <- B < 0.9 & B > -0.9
  A.auc <- A[index.B.range]
  B.auc <- B[index.B.range]

  bauc <- v.area(A.auc, B.auc)

  # Print results .........................................................

  if(results.bauc==TRUE){
  cat("\n")
  cat("Summary results for the Bayesian Area Under the Curve (BAUC)","\n")
  cat("------------------------------------------------------------","\n")

  bauc.ci <- quantile(bauc, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
  print(bauc.ci)
  cat("------------------------------------------------------------","\n")
  cat("\n")
  }

  if(results.bsroc == TRUE)
  {
    level.names <- paste(paste("TPR",level*100, sep=" "),"%", sep="")
    names(dat.lin) <- c("FPR", level.names)

    cat("Summary results for the BSROC curve","\n")
    cat("------------------------------------------------------------","\n")
    print(dat.bsroc, "\n")
  }

  post.bauc <- ggplot(data.frame(bauc), aes(x = bauc)) +
               geom_histogram(colour = "black", fill = "dodgerblue", binwidth = binwidth.p) +
               xlab("Bayesian Area Under the Curve") +
               geom_vline(xintercept = median(bauc)) +
               geom_vline(xintercept = quantile(bauc, prob = c(0.025, 0.975)), linetype = "longdash")+
               ggtitle("Posterior Distribution")

  if(plot.post.bauc == TRUE)
  {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2)))

  vplayout <- function(x, y)
    viewport(layout.pos.row = x,
             layout.pos.col = y)

  print(bsroc.plot, vp = vplayout(1, 1))
  print(post.bauc, vp = vplayout(1, 2))
  } else
  {
    print(bsroc.plot)
  }

  return()
}


