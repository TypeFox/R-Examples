riskchart <- function(pr, ac, ri=NULL,
                      title=NULL,
                      title.size=10,
                      show.legend=TRUE,
                      optimal=NULL, optimal.label="",
                      chosen=NULL, chosen.label="",
                      include.baseline=TRUE,
                      dev="", filename="",
                      show.knots=NULL,
                      show.lift=TRUE,
                      show.precision=TRUE,
                      show.maximal=TRUE,
                      risk.name="Risk",
                      recall.name="Recall",
                      precision.name="Precision",
                      thresholds=NULL,
                      legend.horiz=TRUE)
{

  # Riskcharts are only implemented for binary classificaton.
  
  if (length(unique(ac)) > 2)
    stop("Risk charts are for binary classification models only.")
  
  # 121209 Initialise variables otherwise appearing to be unintialised
  # to R CMD check.

  x <- y <- NULL

  # ggplot2 version of the risk chart.

  if (! requireNamespace("ggplot2")) stop(Rtxt("riskchart requires the ggplot2 package"))

  data <- evaluateRisk(pr, ac, ri)
  
  stopifnot(c("Caseload", "Recall", "Precision") %in% colnames(data)) 

  # Is Risk included in the supplied data?
  
  include.risk <- 'Risk' %in% colnames(data)
  
  score <- as.numeric(row.names(data))
  locateScores <- function(s) 
  { 
    # Convert score to caseload (percentile) for plotting
    idx <- findInterval(s, score)
    idx[idx==0] <- 1
    data$Caseload[idx]
  } 
  scores <- data.frame(ticks=seq(0.1, 0.9, by = 0.1))
  scores$pos <- locateScores(scores$ticks)

  # Calculate AUC for each curve. The calculateAUC is for the whole
  # plot, rather than just under the maximal possible curve. So we
  # need to refine it to be relative to the maximal curve, as
  # suggested by Hari Koesmarno, 9 Aug 2013. I suspect I should always
  # show the AUC wrt the maximal obtainable rather than wrt to the
  # whole plot but for now I'll show the old if show.maximal is
  # selected, as a method of backward compatibility.
  
  recall.auc <- with(data, calculateAUC(Caseload, Recall))
  recall.auc <- round(100*ifelse(show.maximal,
                                 (2*recall.auc)/(2-data$Precision[1]),
                                 recall.auc))
  recall.name <- sprintf("%s (%s%%)", recall.name, recall.auc)

  if (include.risk)
  {
    # Should this be adjusted too wrt the maximal obtainable? I think
    # so. By the time we get to maximal (x=sr and y=100), then we have
    # the maximal risk covered. However, given each case is not equal
    # in terms of risk, we potentially get odd numbers >100% - maybe
    # not, actually.
    risk.auc <- with(data, calculateAUC(Caseload, Risk))
    risk.auc <- round(100*ifelse(show.maximal,
                                 (2*risk.auc)/(2-data$Precision[1]),
                                 risk.auc))
    risk.name <- sprintf("%s (%s%%)", risk.name, risk.auc)
  }
  
  # How to get recall.name etc in the following - by default it is not
  # defined?

  p <- ggplot2::ggplot(data, ggplot2::aes(x=100*Caseload))
  p <- p + ggplot2::scale_colour_manual(breaks=c("Recall", "Risk", "Precision"),
                                        labels=c(recall.name, risk.name, precision.name),
                                        values=c("blue", "forestgreen", "red"))
  p <- p + ggplot2::scale_linetype_manual(breaks=c("Recall", "Risk", "Precision"),
                                          labels=c(recall.name, risk.name, precision.name),
                                          values=c(3, 1, 2))

  if (show.maximal)
  {
    # This is the maximal performance possible from any model - a
    # straight line 100% performance at a caseload including only
    # productive cases.
    
    mx <- data.frame(x=c(0, 100*data$Precision[1], 100), y=c(0, 100, 100))
    p <- p + ggplot2::geom_line(data=mx, ggplot2::aes(x, y), colour="grey")
  }
  
  p <- p + ggplot2::geom_line(ggplot2::aes(y=100*Recall, colour="Recall", linetype="Recall"))
  p <- p + ggplot2::geom_line(ggplot2::aes(y=100*Precision, colour="Precision", linetype="Precision"))
  p <- p + ggplot2::geom_text(data=scores, ggplot2::aes(x=100*pos, y=102, label=ticks), size=2)
  p <- p + ggplot2::annotate("text", x=0, y=110, label="Risk Scores", hjust=0, size=3)
  p <- p + ggplot2::geom_line(ggplot2::aes(y=100*Caseload))
  p <- p + ggplot2::ggtitle(title) + ggplot2::xlab("Caseload (%)") + ggplot2::ylab("Performance (%)")
  if (legend.horiz)
    p <- p + ggplot2::theme(legend.title=ggplot2::element_blank(),
                            plot.title=ggplot2::element_text(size=title.size),
                            legend.justification=c(0, 0),
                            legend.position=c(.10, 0.01),
                            legend.direction="horizontal")
  else
    p <- p + ggplot2::theme(legend.title=ggplot2::element_blank(),
                            plot.title=ggplot2::element_text(size=title.size),
                            legend.justification=c(0, 0),
                            legend.position=c(.30, .02))
  p <- p + ggplot2::guides(colour=ggplot2::guide_legend(keywidth=3, labels=1:3))
  p <- p + ggplot2::scale_x_continuous(breaks=seq(0, 100, 20))
  p <- p + ggplot2::scale_y_continuous(breaks=seq(0, 100, 20))
  
  if (include.risk)
    p <- p + ggplot2::geom_line(ggplot2::aes(y=100*Risk, colour="Risk", linetype="Risk"))

  if (include.baseline && show.precision)
    p <- p + ggplot2::geom_text(data=data.frame(x=100, y=100*data$Precision[1]+4,
                           text=sprintf("%0.0f%%", 100*data$Precision[1])),
                         ggplot2::aes(x, y, label=text), size=3)

  if (show.lift)
  {
    base.lift <- 100*data$Precision[1]
    lift.seq <- seq(base.lift, 100, base.lift)
    lifts <- data.frame(x=110, y=lift.seq, label=lift.seq/base.lift)
    # 140906 R-Devel CMD check complains label not defined in:
    # p <- p + ggplot2::geom_text(data=lifts, ggplot2::aes(x, y, label=label), size=3)
    # So use aes_string
    p <- p + ggplot2::geom_text(data=lifts,
                                ggplot2::aes_string("x", "y", label="label"), size=3)
    p <- p + ggplot2::annotate("text", x=110, y=101.5, label="Lift", size=3)
  }

  if (! is.null(thresholds))
    p <- p + ggplot2::geom_vline(xintercept=100*locateScores(thresholds),
                        linetype="twodash", color="grey")

  p <- p + ggplot2::theme(legend.background=ggplot2::element_rect(fill=NA),
                 legend.key=ggplot2::element_rect(fill=NA, linetype=0))
  
  return(p)
}
