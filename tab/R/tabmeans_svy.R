tabmeans.svy <- function(x, y, svy, latex = FALSE, xlevels = NULL, yname = "Y variable", 
                         test = "Wald", decimals = 1, p.decimals = c(2, 3), p.cuts = 0.01, 
                         p.lowerbound = 0.001, p.leading0 = TRUE, p.avoid1 = FALSE, n.column = FALSE,
                         n.headings = TRUE, bold.colnames = TRUE, bold.varnames = FALSE, 
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
  if (! test %in% c("Wald", "LRT")) {
    stop("For test input, please enter 'Wald' or 'LRT'")
  }
  if (!is.numeric(decimals)) {
    stop("For decimals input, please enter numeric value")
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
  
  # Get means and SE's overall and by levels of x, and get sample size in each x
  totmean <- svymean(y, design = svy2)
  means <- svyby(~y, by = ~x, FUN = svymean, design = svy2)
  ns <- tapply(X = y, INDEX = x, FUN = length)
  
  # Add mean (SE) values to table
  tbl[1, 1] <- paste(yname, ", M (SE)", sep = "")
  tbl[1, 2] <- sprintf("%.0f", sum(ns))
  tbl[1, 3] <- paste(sprintf(spf, totmean), " (", sprintf(spf, sqrt(attr(totmean, "var"))), ")", sep = "")
  tbl[1, 4:(ncol(tbl)-1)] <- paste(sprintf(spf, means$"y"), " (", sprintf(spf, means$"se"), ")", sep = "")
  
  # ANOVA
  fit1 <- svyglm(y ~ 1, design = svy2)
  fit2 <- svyglm(y ~ as.factor(x), design = svy2)
  pval <- anova(fit1, fit2, method = test)$p
  tbl[1, ncol(tbl)] <- formatp(p = pval, cuts = p.cuts, decimals = p.decimals, lowerbound = p.lowerbound, leading0 = p.leading0, avoid1 = p.avoid1)
  
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