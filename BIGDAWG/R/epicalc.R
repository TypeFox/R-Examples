#' Case-Control Odds ratio calculation and graphing
#'
#' cci function port epicalc version 2.15.1.0 (Virasakdi Chongsuvivatwong, 2012)
#' @param caseexp Number of cases exposed
#' @param controlex Number of controls exposed
#' @param casenonex Number of cases not exosed
#' @param controlnonex Number of controls not exposed
#' @param cctable A 2-by-2 table. If specified, will supercede the outcome and exposure variables
#' @param graph If TRUE (default), produces an odds ratio plot
#' @param design Specification for graph; can be "case control","case-control", "cohort" or "prospective"
#' @param main main title of the graph
#' @param xlab label on X axis
#' @param ylab label on Y axis
#' @param xaxis two categories of exposure in graph
#' @param yaxis two categories of outcome in graph
#' @param alpha level of significance
#' @param fisher.or whether odds ratio should be computed by the exact method
#' @param exact.ci.or whether confidence limite of the odds ratio should be computed by the exact method
#' @param decimal number of decimal places displayed
#' @note This function is for internal BIGDAWG use only.
cci <- function (caseexp, controlex, casenonex, controlnonex, cctable = NULL, graph = TRUE, design = "cohort", main, xlab, ylab, xaxis, yaxis, alpha = 0.05, fisher.or = FALSE, exact.ci.or = TRUE, decimal = 2) {

  if (is.null(cctable)) {
    frame <- cbind(Outcome <- c(1, 0, 1, 0), Exposure <- c(1, 
                                                           1, 0, 0), Freq <- c(caseexp, controlex, casenonex, 
                                                                               controlnonex))
    Exposure <- factor(Exposure)
    expgrouplab <- c("Non-exposed", "Exposed")
    levels(Exposure) <- expgrouplab
    Outcome <- factor(Outcome)
    outcomelab <- c("Negative", "Positive")
    levels(Outcome) <- outcomelab
    table1 <- xtabs(Freq ~ Outcome + Exposure, data = frame)
  }
  else {
    table1 <- as.table(get("cctable"))
  }
  fisher <- fisher.test(table1)
  caseexp <- table1[2, 2]
  controlex <- table1[1, 2]
  casenonex <- table1[2, 1]
  controlnonex <- table1[1, 1]
  se.ln.or <- sqrt(1/caseexp + 1/controlex + 1/casenonex + 
                     1/controlnonex)
  if (!fisher.or) {
    or <- caseexp/controlex/casenonex * controlnonex
    p.value <- chisq.test(table1, correct = FALSE)$p.value
  }
  else {
    or <- fisher$estimate
    p.value <- fisher$p.value
  }
  if (exact.ci.or) {
    ci.or <- as.numeric(fisher$conf.int)
  }
  else {
    ci.or <- or * exp(c(-1, 1) * qnorm(1 - alpha/2) * se.ln.or)
  }
  if (graph == TRUE) {
    caseexp <- table1[2, 2]
    controlex <- table1[1, 2]
    casenonex <- table1[2, 1]
    controlnonex <- table1[1, 1]
    if (!any(c(caseexp, controlex, casenonex, controlnonex) < 
               5)) {
      if (design == "prospective" || design == "cohort" || 
            design == "cross-sectional") {

        if (missing(main)) 
          main <- "Odds ratio from prospective/X-sectional study"
        if (missing(xlab)) 
          xlab <- ""
        if (missing(ylab)) 
          ylab <- paste("Odds of being", ifelse(missing(yaxis), 
                                                "a case", yaxis[2]))
        if (missing(xaxis)) 
          xaxis <- c("non-exposed", "exposed")
        axis(1, at = c(0, 1), labels = xaxis)
      }
      else {

        if (missing(main)) 
          main <- "Odds ratio from case control study"
        if (missing(ylab)) 
          ylab <- "Outcome category"
        if (missing(xlab)) 
          xlab <- ""
        if (missing(yaxis)) 
          yaxis <- c("Control", "Case")
        axis(2, at = c(0, 1), labels = yaxis, las = 2)
        mtext(paste("Odds of ", ifelse(xlab == "", "being exposed", 
                                       paste("exposure being", xaxis[2]))), side = 1, 
              line = ifelse(xlab == "", 2.5, 1.8))
      }
      title(main = main, xlab = xlab, ylab = ylab)
    }
  }
  if (!fisher.or) {
    results <- list(or.method = "Asymptotic", or = or, se.ln.or = se.ln.or, 
                    alpha = alpha, exact.ci.or = exact.ci.or, ci.or = ci.or, 
                    table = table1, decimal = decimal)
  }
  else {
    results <- list(or.method = "Fisher's", or = or, alpha = alpha, 
                    exact.ci.or = exact.ci.or, ci.or = ci.or, table = table1, 
                    decimal = decimal)
  }
  class(results) <- c("cci", "cc")
  return(results)
}

#' Creation of a 2x2 table using the indicated orientation.
#'
#' make2x2 function port epicalc version 2.15.1.0 (Virasakdi Chongsuvivatwong, 2012)
#' @param caseexp Number of cases exposed
#' @param controlex Number of controls exposed
#' @param casenonex Number of cases not exosed
#' @param controlnonex Number of controls not exposed
#' @note This function is for internal BIGDAWG use only.
make2x2 <- function (caseexp, controlex, casenonex, controlnonex)  {
  
  table1 <- c(controlnonex, casenonex, controlex, caseexp)
  dim(table1) <- c(2, 2)
  rownames(table1) <- c("Non-diseased", "Diseased")
  colnames(table1) <- c("Non-exposed", "Exposed")
  attr(attr(table1, "dimnames"), "names") <- c("Outcome", "Exposure")
  table1

}