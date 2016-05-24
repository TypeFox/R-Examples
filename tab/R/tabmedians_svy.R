tabmedians.svy <- function(svy, x, y, latex = FALSE, xlevels = NULL, yname = "Y variable",
                           test = "wilcoxon", decimals = 1, p.include = TRUE, p.decimals = c(2, 3),
                           p.cuts = 0.01, p.lowerbound = 0.001, p.leading0 = TRUE, p.avoid1 = FALSE,
                           n.column = FALSE, n.headings = TRUE, parenth = "iqr", text.label = NULL,
                           parenth.sep = "-", bold.colnames = TRUE, bold.varnames = FALSE,
                           variable.colname = "Variable") {
  
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
  if (! test %in% c("wilcoxon", "vanderWaerden", "median", "KruskalWallis")) {
    stop("For test input, please enter a possible value for the 'test' input of the 
         svyranktest function in the survey package: 'wilcoxon', 'vanderWaerden', 
         'median', or 'KruskalWallis'. See documentation for tabmedians.svy and 
         svyranktest for details.")
  }
  if (!is.numeric(decimals)) {
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
  if (!is.logical(n.column)) {
    stop("For n.column input, please enter TRUE or FALSE")
  }
  if (!is.logical(n.headings)) {
    stop("For n.headings input, please enter TRUE or FALSE")
  }
  if (! parenth %in% c("minmax", "range", "q1q3", "iqr", "none")) {
    stop("For parenth input, please enter one of the following: 'minmax', 'range', 'q1q3', 'iqr', 'none'")
  }
  if (!is.null(text.label) && !is.character(text.label)) {
    stop("For text.label input, please enter something like 'Median (IQR)' or just leave it unspecified")
  }
  if (!is.character(parenth.sep)) {
    stop("For parenth.sep input, please enter a character string")
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
  
  # Convert decimals to variable for sprintf
  spf <- paste("%0.", decimals, "f", sep = "")
  
  # Save x and y as character strings
  xstring <- x
  ystring <- y
  
  # Extract vectors x and y
  x <- svy$variables[, xstring]
  y <- svy$variables[, ystring]
  
  # Update survey object to include y and x explicitly
  svy2 <- update(svy, y = y, x = x)
  
  # Drop missing values if present
  locs <- which(!is.na(x) & !is.na(y))
  if (length(locs) < nrow(svy2)) {
    svy2 <- subset(svy2, !is.na(x) & !is.na(y))
    x <- svy2$variables[, xstring]
    y <- svy2$variables[, ystring]
  }
  
  # Get unique values of x
  xvals <- sort(unique(x))
  
  # If xlevels unspecified, set to actual values
  if (is.null(xlevels)) {
    xlevels <- xvals
  }
  
  # Initialize table
  tbl <- matrix("", nrow = 1, ncol = length(xlevels) + 4)
  
  # Get medians and values for parentheses, and add variable name to 1st cell entry to table
  #medians <- c(svyquantile(~y, svy2, 0.5), svyby(~y, ~x, svy2, svyquantile, quantile = 0.5, keep.var = FALSE)[, 2])
  overallmedian <- svyquantile(~y, svy2, 0.5)
  medians <- svyby(~y, ~x, svy2, svyquantile, quantile = 0.5, keep.var = FALSE)[, 2]
  ns <- tapply(X = y, INDEX = x, FUN = length)
  if (parenth == "minmax") {
    parent1 <- tapply(X = y, INDEX = x, FUN = min)
    parent2 <- tapply(X = y, INDEX = x, FUN = max)
    parent <- paste(sprintf(spf, parent1), parenth.sep, sprintf(spf, parent2), sep = "")
    if (is.null(text.label)) {
      text.label <- "Median (Min-Max)"
    }
    tbl[1, 1] <- paste(yname, ", ", text.label, sep = "")
    tbl[1, 2] <- sprintf("%.0f", sum(ns))
    tbl[1, 3] <- paste(sprintf(spf, overallmedian), " (", sprintf(spf, min(y)), parenth.sep, sprintf(spf, max(y)), ")", sep = "")
    tbl[1, 4:(ncol(tbl)-1)] <- paste(sprintf(spf, medians), " (", parent, ")", sep = "")
  } else if (parenth == "range") {
    parent1 <- tapply(X = y, INDEX = x, FUN = min)
    parent2 <- tapply(X = y, INDEX = x, FUN = max)
    parent <- paste(sprintf(spf, parent2 - parent1))
    if (is.null(text.label)) {
      text.label <- "Median (Range)"
    }
    tbl[1, 1] <- paste(yname, ", ", text.label, sep = "")
    tbl[1, 2] <- sprintf("%.0f", sum(ns))
    tbl[1, 3] <- paste(sprintf(spf, overallmedian), " (", sprintf(spf, max(y) - min(y)), ")", sep = "")
    tbl[1, 4:(ncol(tbl)-1)] <- paste(sprintf(spf, medians), " (", parent, ")", sep = "")
  } else if (parenth == "q1q3") {
    parent1 <- svyby(~y, ~x, svy2, svyquantile, quantile = 0.25, keep.var = FALSE)[, 2]
    parent2 <- svyby(~y, ~x, svy2, svyquantile, quantile = 0.75, keep.var = FALSE)[, 2]
    parent <- paste(sprintf(spf, parent1), parenth.sep, sprintf(spf, parent2), sep = "")
    if (is.null(text.label)) {
      text.label <- "Median (Q1-Q3)"
    }
    tbl[1, 1] <- paste(yname, ", ", text.label, sep = "")
    tbl[1, 2] <- sprintf("%.0f", sum(ns))
    tbl[1, 3] <- paste(sprintf(spf, overallmedian), " (", sprintf(spf, svyquantile(~y, svy2, 0.25)), 
                       parenth.sep, sprintf(spf, svyquantile(~y, svy2, 0.75)), ")", sep = "")
    tbl[1, 4:(ncol(tbl)-1)] <- paste(sprintf(spf, medians), " (", parent, ")", sep = "")
  } else if (parenth == "iqr") {
    parent1 <- svyby(~y, ~x, svy2, svyquantile, quantile = 0.25, keep.var = FALSE)[, 2]
    parent2 <- svyby(~y, ~x, svy2, svyquantile, quantile = 0.75, keep.var = FALSE)[, 2]
    parent <- paste(sprintf(spf, parent2 - parent1))
    if (is.null(text.label)) {
      text.label <- "Median (IQR)"
    }
    tbl[1, 1] <- paste(yname, ", ", text.label, sep = "")
    tbl[1, 2] <- sprintf("%.0f", sum(ns))
    tbl[1, 3] <- paste(sprintf(spf, overallmedian), " (", sprintf(spf, svyquantile(~y, svy2, 0.75) - svyquantile(~y, svy2, 0.25)), ")", sep = "")
    tbl[1, 4:(ncol(tbl)-1)] <- paste(sprintf(spf, medians), " (", parent, ")", sep = "")
  } else if (parenth == "none") {
    if (is.null(text.label)) {
      text.label <- "Median"
    }
    tbl[1, 1] <- paste(yname, ", ", text.label, sep = "")
    tbl[1, 2] <- sprintf("%.0f", sum(ns))
    tbl[1, 3] <- sprintf(spf, overallmedian)
    tbl[1, 4:(ncol(tbl)-1)] <- sprintf(spf, medians)
  }
  
  # Add p-value from statistical test
  if (p.include == TRUE) {
    p <- svyranktest(y ~ x, svy2, test)$p.value
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
  
  # Drop p column if requested
  if (p.include == FALSE) {
    tbl <- tbl[, -which(colnames(tbl) == "P")]
  }
  
  # If latex is TRUE, do some re-formatting
  if (latex == TRUE) {
    plocs <- which(substr(tbl[, "P"], 1, 1) == "<")
    if (length(plocs) > 0) {
      tbl[plocs, "P"] <- paste("$<$", substring(tbl[plocs, "P"], 2), sep = "")
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