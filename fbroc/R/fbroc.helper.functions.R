# transforms list to roc data.frame
c.list.to.roc <- function(input) {
  roc.frame <- as.data.frame(input[1:3])
  names(roc.frame) <- c("TPR", "FPR", "threshold")
  return(roc.frame)
}

# Check if number is a single numeric between 0 and 1
validate.single.numeric <- function(number, var.name) {
  if (is.null(number)) stop(paste("Please pass ", var.name, " to perf.roc!", sep = ""))
  if (class(number) != "numeric") stop(paste(var.name, " must be numeric!", sep = ""))
  if (length(number) != 1) stop(paste(var.name, " must have length 1!", sep = ""))
  if (is.na(number)) stop(paste(var.name, " is NA!", sep = ""))
  if ((number < 0) | (number > 1)) stop(paste(var.name, " must be in [0, 1]!", sep = ""))
  return(number)
}


fbroc.plot.base <- function(plot.frame) {
  roc.plot <- ggplot(data = plot.frame, aes(x = FPR, y = TPR)) +               
              ggtitle("ROC Curve") + xlab("False Positive Rate") +
              ylab("True Positive Rate") + theme_bw() +
              theme(title = element_text(size = 22),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                axis.text.x = element_text(size = 16),
                axis.text.y = element_text(size = 16))
  return(roc.plot)
}

# adds confidence region to roc curve plot
fbroc.plot.add.conf <- function(roc1, conf.level = conf.level, steps = steps, fill = fill) {
  conf.frame <- conf(roc1, conf.level = conf.level, steps = steps)
  conf.frame$Segment <- 1
 
  geom_ribbon(data = conf.frame, fill = fill, alpha = 0.5,
              aes(y = NULL, ymin = Lower.TPR, ymax = Upper.TPR))
}

# add performance metric visualization to roc plot (paired roc curve)
fbroc.plot.add.metric.paired <- function(roc.plot,show.metric, perf, col1, col2) {
  if (show.metric == "tpr") {
    extra.frame <- data.frame(FPR = perf$params, 
                              TPR = c(perf$Observed.Performance.Predictor1, 
                                      perf$Observed.Performance.Predictor2),
                              Segment = 1,
                              lower = c(perf$CI.Performance.Predictor1[1], 
                                        perf$CI.Performance.Predictor2[1]),
                              upper = c(perf$CI.Performance.Predictor1[2],
                                        perf$CI.Performance.Predictor2[2]))
  
    roc.plot <- roc.plot + geom_errorbar(data = extra.frame, width = 0.02, size = 1.25,
                              aes(ymin = lower, ymax = upper), col = c(col1, col2), alpha = 0.7) + 
      geom_point(data = extra.frame, size = 4, col = c(col1, col2))
  }
  if (show.metric == "fpr") {
    extra.frame <- data.frame(TPR = perf$params, 
                              FPR = c(perf$Observed.Performance.Predictor1, 
                                      perf$Observed.Performance.Predictor2),
                              Segment = 1,
                              lower = c(perf$CI.Performance.Predictor1[1], 
                                        perf$CI.Performance.Predictor2[1]),
                              upper = c(perf$CI.Performance.Predictor1[2],
                                        perf$CI.Performance.Predictor2[2]))
    roc.plot <- roc.plot + geom_errorbarh(data = extra.frame, height = 0.02, size = 1.25,
                              aes(xmin = lower, xmax = upper), col = c(col1, col2), alpha = 0.7) +
      geom_point(data = extra.frame, size = 4)
  }
  return(roc.plot)
}

# add performance metric visualization to roc plot (roc curve)
fbroc.plot.add.metric <- function(roc.plot, show.metric, perf, col) {
  if (show.metric == "tpr") {
    extra.frame <- data.frame(FPR = perf$params, TPR = perf$Observed.Performance, Segment = 1,
                              lower = perf$CI.Performance[1], upper = perf$CI.Performance[2])
    roc.plot <- roc.plot + geom_errorbar(data = extra.frame, width = 0.02, size = 1.25,
                                         aes(ymin = lower, ymax = upper)) + 
      geom_point(data = extra.frame, size = 4)
  }
  if (show.metric == "fpr") {
    extra.frame <- data.frame(TPR = perf$params, FPR = perf$Observed.Performance, Segment = 1,
                              lower = perf$CI.Performance[1], upper = perf$CI.Performance[2])
    roc.plot <- roc.plot + geom_errorbarh(data = extra.frame, height = 0.02, size = 1.25,
                                          aes(xmin = lower, xmax = upper)) +
      geom_point(data = extra.frame, size = 4)
  }
  return(roc.plot)
}

