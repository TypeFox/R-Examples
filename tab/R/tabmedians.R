tabmedians <- function(x, y, latex = FALSE, xlevels = NULL, yname = NULL, quantiles = NULL,
                       quantile.vals = FALSE, parenth = "iqr", text.label = NULL, parenth.sep = "-",
                       decimals = NULL, p.include = TRUE, p.decimals = c(2, 3), p.cuts = 0.01, 
                       p.lowerbound = 0.001, p.leading0 = TRUE, p.avoid1 = FALSE, 
                       overall.column = TRUE, n.column = FALSE, n.headings = TRUE, 
                       bold.colnames = TRUE, bold.varnames = FALSE, variable.colname = "Variable") {
  
  # If yname unspecified, use variable name
  if (is.null(yname)) {
    yname <- deparse(substitute(y))
    if (grepl("\\$", yname)) {
      yname <- strsplit(yname, "\\$")[[1]][2]
    }
  }
  
  # If any inputs are not correct class, return error
  if (!is.logical(latex)) {
    stop("For latex input, please enter TRUE or FALSE")
  }
  if (!is.null(xlevels) && !is.character(xlevels)) {
    stop("For xlevels input, please enter vector of character strings")
  }
  if (!is.character(yname)) {
    stop("For yname input, please enter character string")
  }
  if (!is.null(quantiles) && ! (is.numeric(quantiles) & length(quantiles) == 1 && quantiles > 1 & quantiles == round(quantiles))) {
    stop("For quantiles input, please enter a whole number greater than 1")
  }
  if (!is.logical(quantile.vals)) {
    stop("For quantile.vals input, please enter TRUE or FALSE")
  }
  if (! parenth %in% c("none", "iqr", "range", "minmax", "q1q3")) {
    stop("For parenth input, please enter 'none', 'iqr', 'range', 'minmax', or 'q1q3'")
  }
  if (!is.null(text.label) && !is.character(text.label)) {
    stop("For text.label input, please enter a character string or just leave it unspecified. Use 'none' to request no label")
  }
  if (!is.character(parenth.sep)) {
    stop("For parenth.sep input, please enter a character string (only used if parenth is set to 'minmax' or 'q1q3'; usually '-' or ', ')")
  }
  if (!is.null(decimals) && !is.numeric(decimals)) {
    stop("For decimals input, please enter numeric value")
  }
  if (!is.logical(p.include)) {
    stop("For p.include input, please enter TRUE or FALSE")
  }
  if (!is.numeric(p.decimals)) {
    stop("For p.decimals input, please enter numeric value or vector")
  }
  if (!is.numeric(p.cuts)) {  
    stop("For p.cuts input, please enter numeric value or vector")
  }
  if (!is.numeric(p.lowerbound)) {
    stop("For p.lowerbound input, please enter numeric value")
  }
  if (!is.logical(p.leading0)) {
    stop("For p.leading0 input, please enter TRUE or FALSE")
  }
  if (!is.logical(p.avoid1)) {
    stop("For p.avoid1 input, please enter TRUE or FALSE")
  }
  if (!is.logical(overall.column)) {
    stop("For overall.column input, please enter TRUE or FALSE")
  }
  if (!is.logical(n.column)) {
    stop("For n.column input, please enter TRUE or FALSE")
  }
  if (!is.logical(n.headings)) {
    stop("For n.headings input, please enter TRUE or FALSE")
  }
  if (!is.logical(bold.colnames)) {
    stop("For bold.colnames input, please enter TRUE or FALSE")
  }
  if (!is.logical(bold.varnames)) {
    stop("For bold.varnames input, please enter TRUE or FALSE")
  }
  if (!is.character(variable.colname)) {
    stop("For variable.colname input, please enter a character string")
  }
  
  # Drop missing values
  locs.complete <- which(!is.na(x) & !is.na(y))
  x <- x[locs.complete]
  y <- y[locs.complete]
  
  # Create quantiles if necessary
  if (!is.null(quantiles)) {
    x <- cut(x = x, breaks = quantile(x, probs = seq(0, 1, 1/quantiles)), include.lowest = TRUE, right = TRUE, dig.lab = 3)
  }
  
  # Get medians and sample sizes
  medians <- tapply(X = y, INDEX = x, FUN = median)
  ns <- tapply(X = y, INDEX = x, FUN = length)
  
  # If decimals is unspecified, set to appropriate value
  if (is.null(decimals)) {
    max.median <- max(abs(medians))
    if (max.median >= 1000) {
      decimals <- 0
    } else if (max.median < 1000 & max.median >= 10) {
      decimals <- 1
    } else if (max.median < 10 & max.median >= 0.1) {
      decimals <- 2
    } else if (max.median < 0.1 & max.median >= 0.01) {
      decimals <- 3
    } else if (max.median < 0.01 & max.median >= 0.001) {
      decimals <- 4
    } else if (max.median < 0.001 & max.median >= 0.0001) {
      decimals <- 5
    } else if (max.median < 0.0001) {
      decimals <- 6
    }
  }
  
  # Convert decimals to variable for sprintf
  spf <- paste("%0.", decimals, "f", sep = "")
  
  # Get unique values of x
  xvals <- sort(unique(x))
  
  # If xlevels unspecified, set to actual values
  if (is.null(xlevels)) {
    if (!is.null(quantiles)) {
      if (quantile.vals == TRUE) {
        xlevels <- paste("Q", 1:length(xvals), " ", as.character(xvals), sep = "")
      } else {
        xlevels <- paste("Q", 1:length(xvals), sep = "")
      }
    } else {
      xlevels <- as.character(xvals)
    }
  }
  
  # Initialize table
  tbl <- matrix("", nrow = 1, ncol = length(xlevels) + 4)
  
  # Figure out text.label
  if (is.null(text.label)) {
    if (parenth == "iqr") {
      text.label <- ", Median (IQR)"
    } else if (parenth == "Range") {
      text.label <- ", Median (range)"
    } else if (parenth == "minmax") {
      text.label <- ", Median (Min-Max)"
    } else if (parenth == "q1q3") {
      text.label <- ", Median (Q1-Q3)"
    } else if (parenth == "none") {
      text.label <- ", Median"
    }
  } else if (text.label == "none") {
    text.label <- NULL
  }
  
  # Add variable name and n column to table
  tbl[1, 1] <- paste(yname, text.label, sep = "")
  tbl[1, 2] <- sprintf("%.0f", length(y))
  
  # Add Median (parenth) overall and in each x group
  if (parenth == "minmax") {
    parent1 <- tapply(X = y, INDEX = x, FUN = min)
    parent2 <- tapply(X = y, INDEX = x, FUN = max)
    parent <- paste(sprintf(spf, parent1), parenth.sep, sprintf(spf, parent2), sep = "")
    tbl[1, 3] <- paste(sprintf(spf, median(y)), " (", sprintf(spf, min(y)), parenth.sep, sprintf(spf, max(y)), ")", sep = "")
    tbl[1, 4:(ncol(tbl)-1)] <- paste(sprintf(spf, medians), " (", parent, ")", sep = "")
  } else if (parenth == "range") {
    parent1 <- tapply(X = y, INDEX = x, FUN = min)
    parent2 <- tapply(X = y, INDEX = x, FUN = max)
    parent <- paste(sprintf(spf, parent2 - parent1))
    tbl[1, 3] <- paste(sprintf(spf, median(y)), " (", sprintf(spf, max(y) - min(y)), ")", sep = "")
    tbl[1, 4:(ncol(tbl)-1)] <- paste(sprintf(spf, medians), " (", parent, ")", sep = "")
  } else if (parenth == "q1q3") {
    parent1 <- tapply(X = y, INDEX = x, FUN = function(x) quantile(x, probs = 0.25))
    parent2 <- tapply(X = y, INDEX = x, FUN = function(x) quantile(x, probs = 0.75))
    parent <- paste(sprintf(spf, parent1), parenth.sep, sprintf(spf, parent2), sep = "")
    tbl[1, 3] <- paste(sprintf(spf, median(y)), " (", sprintf(spf, quantile(y, probs = 0.25)), parenth.sep, sprintf(spf, quantile(y, 0.75)), ")", sep = "")
    tbl[1, 4:(ncol(tbl)-1)] <- paste(sprintf(spf, medians), " (", parent, ")", sep = "")
  } else if (parenth == "iqr") {
    parent1 <- tapply(X = y, INDEX = x, FUN = function(x) quantile(x, probs = 0.25))
    parent2 <- tapply(X = y, INDEX = x, FUN = function(x) quantile(x, probs = 0.75))
    parent <- paste(sprintf(spf, parent2 - parent1))
    tbl[1, 3] <- paste(sprintf(spf, median(y)), " (", sprintf(spf, quantile(y, probs = 0.75) - quantile(y, probs = 0.25)), ")", sep = "")
    tbl[1, 4:(ncol(tbl)-1)] <- paste(sprintf(spf, medians), " (", parent, ")", sep = "")
  } else if (parenth == "none") {
    tbl[1, 3] <- sprintf(spf, median(y))
    tbl[1, 4:(ncol(tbl)-1)] <- sprintf(spf, medians)
  }
  
  # Add p-value from statistical test depending on number of levels of x
  if (p.include == TRUE) {
    
    if (length(xlevels) == 2) {
      
      # Mann-Whitney U test a.k.a. Wilcoxon rank-sum test
      p <- wilcox.test(y ~ x)$p.value
      
    } else {
      
      # Kruskal-Wallis rank-sum test
      p <- kruskal.test(y ~ as.factor(x))$p.value
      
    }
    
  } else {
    p <- NA
  }
  
  # Add p-value from t-test
  if (p.include == TRUE) {
    tbl[1, ncol(tbl)] <- formatp(p = p, cuts = p.cuts, decimals = p.decimals, lowerbound = p.lowerbound, leading0 = p.leading0, avoid1 = p.avoid1)
  }
  
  # Add column names, with sample sizes for each group if requested
  if (n.headings == FALSE) {
    colnames(tbl) <- c(variable.colname, "N", "Overall", xlevels, "P")
  } else {
    colnames(tbl) <- c(variable.colname, "N", paste(c("Overall", xlevels), " (n = ", c(sum(ns), ns), ")", sep = ""), "P")
  }
  
  # Drop N column if requested
  if (n.column == FALSE) {
    tbl <- tbl[, -which(colnames(tbl) == "N"), drop = FALSE]
  }
  
  # Drop overall column if requested
  if (overall.column == FALSE) {
    tbl <- tbl[, -grep("^Overall", colnames(tbl)), drop = FALSE]
  }
  
  # Drop p column if requested
  if (p.include == FALSE) {
    tbl <- tbl[, -which(colnames(tbl) == "P"), drop = FALSE]
  }
  
  # If latex is TRUE, do some re-formatting
  if (latex == TRUE) {
    if (p.include == TRUE) {
      plocs <- which(substr(tbl[, "P"], 1, 1) == "<")
      if (length(plocs) > 0) {
        tbl[plocs, "P"] <- paste("$<$", substring(tbl[plocs, "P"], 2), sep = "")
      }
    }
    if (bold.colnames == TRUE) {
      colnames(tbl) <- paste("$\\textbf{", colnames(tbl), "}$", sep = "")
    }
    if (bold.varnames == TRUE) {
      tbl[1, 1] <- paste("$\\textbf{", tbl[1, 1], "}$")
    }
  }
  
  # Return table
  return(tbl)
  
}