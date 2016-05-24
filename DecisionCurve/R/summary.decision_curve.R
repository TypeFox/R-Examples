#' Displays a useful description of a decision_curve object
#'
#' @param object decision_curve object to summarise
#' @param ... other arguments ignored (for compatibility with generic)
#' @param measure name of summary measure to print out. For standardized net benefit: "sNB" (default), net benefit: "NB", true positive rate: "TPR", false positive rate: "FPR".
#' @param nround number of decimal places to round (default 3).
#' @method summary decision_curve
#' @import pander
#' @import reshape
#' @examples
#'#helper function
#'
#'#load simulated data
#'data(dcaData)
#'
#'full.model <- decision_curve(Cancer~Age + Female + Smokes + Marker1 + Marker2,
#'data = dcaData,
#'thresholds = seq(0, .4, by = .05),
#'bootstraps = 25)
#'
#'summary(full.model) #outputs standardized net benefit by default
#'
#'summary(full.model, nround = 2, measure = "TPR")
#'
#' @export

summary.decision_curve <- function(object, ..., measure = c("sNB", "NB", "TPR", "FPR"), nround = 3){
  x <- object
  #get measure name for printing
  measure <- match.arg(measure)
  measure.names.df <- data.frame(measure = c("sNB", "NB", "TPR", "FPR"), measure.names = c("Standardized Net Benefit",
                                                                                      "Net Benefit",
                                                                                      "Sensitivity (TPR)",
                                                                                      "1-Specificity (FPR)"))
  measure.name <- as.character(measure.names.df[match(measure, measure.names.df$measure), "measure.names"])

  #if this is true, confidence intervals have been calculated
  conf.int <- ncol(x$derived.data) > 10
  model <- NULL #appease check
  xx.wide <- cast(x$derived.data, thresholds+cost.benefit.ratio~model, value = measure)
  #need to add prob.high risk from the formula and convert to percent
  xx.wide$prob.high.risk <- subset(x$derived.data, !is.element(model, c("All", "None")))$prob.high.risk*100

  #rearrange terms to make sure we have the right ordering
  formula.name <- unique(x$derived.data$model)
  formula.name <- formula.name[!is.element(formula.name, c("None", "All"))]

  xx.wide <- xx.wide[, c("thresholds", "cost.benefit.ratio", "prob.high.risk", "All", formula.name, "None")]



  if( conf.int){

    xx.lower <- cast(x$derived.data, thresholds+cost.benefit.ratio~model, value = paste(measure, "_lower",sep = ""))
    xx.lower$prob.high.risk <- subset(x$derived.data, !is.element(model, c("All", "None")))$prob.high.risk_lower*100
    xx.lower <- xx.lower[, c("thresholds", "cost.benefit.ratio", "prob.high.risk", "All", formula.name, "None")]

    xx.upper <- cast(x$derived.data, thresholds+cost.benefit.ratio~model, value = paste(measure, "_upper",sep = ""))
    xx.upper$prob.high.risk <- subset(x$derived.data, !is.element(model, c("All", "None")))$prob.high.risk_upper*100
    xx.upper <- xx.upper[, c("thresholds", "cost.benefit.ratio", "prob.high.risk", "All", formula.name, "None")]

  }else{
    xx.lower <- NULL
    xx.upper <- NULL
  }

  out <- xx.wide
  out[,-c(1:2)] <- round(out[,-c(1:2)], nround)
  names(out)[1] <- "risk\nthreshold"
  names(out)[2] <- c("cost:benefit\n ratio")
  names(out)[3] <- c("percent\n high risk")



  if(conf.int){
    cat(paste0("\n", measure.name, " (", round(100*x$confidence.intervals),  "% Confidence Intervals):"))

    not.preds <- match(c("risk\nthreshold", "cost:benefit\n ratio",  "None"),
                       names(out))
    n.preds <- ncol(out) - length(not.preds)
    for( i in 1:n.preds){
       out[,-c(not.preds)][i] <- paste0(out[,-c(not.preds)][[i]],
                                 "\n(",
                                 round(c(unlist(xx.lower[,-c(not.preds)][[i]])), nround),
                                 ", ",
                                 round(c(unlist(xx.upper[,-c(not.preds)][[i]])), nround),
                                 ")")
    }
  }else{
    cat(paste0("\n", measure.name, ":"))

  }

  out.table <- pandoc.table(out, split.table = Inf, keep.line.breaks = TRUE)
  out.table
  cat('\n')

  invisible(out.table)
}
