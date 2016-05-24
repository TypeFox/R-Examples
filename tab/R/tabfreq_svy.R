tabfreq.svy <- function(x, y, svy, latex = FALSE, xlevels = NULL, yname = "Y variable",
                        ylevels = NULL, test = "F", decimals = 1, p.decimals = c(2,3), p.cuts = 0.01,
                        p.lowerbound = 0.001, p.leading0 = TRUE, p.avoid1 = FALSE, n.column = FALSE,
                        n.headings = TRUE, compress = FALSE, compress.val = NULL, 
                        bold.colnames = TRUE, bold.varnames = FALSE, bold.varlevels = FALSE, 
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
  if (!is.null(ylevels) && !is.character(ylevels)) {
    stop("For ylevels input, please enter vector of character strings")
  }
  if (! test %in% c("F", "Chisq", "Wald", "adjWald", "lincom", "saddlepoint")) {
    stop("For test input, please enter a possible value for the 'statistic' input of the 
         svychisq function in the survey package: 'F', 'Chisq', 'Wald', 'adjWald', 'lincom', 
         or 'saddlepoint'. See svychisq documentation for details.")
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
  if (!is.logical(compress)) {
    stop("For compress input, please enter TRUE or FALSE")
  }
  if (!is.logical(bold.colnames)) {
    stop("For bold.colnames input, please enter TRUE or FALSE")
  }
  if (!is.logical(bold.varnames)) {
    stop("For bold.varnames input, please enter TRUE or FALSE")
  }
  if (!is.logical(bold.varlevels)) {
    stop("For bold.varlevels input, please enter TRUE or FALSE")
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
  
  # Basic table to get levels of x and y
  counts <- table(y, x)
  
  # If ylevels unspecified, set to actual values
  if (is.null(ylevels)) {
    ylevels <- rownames(counts)
  }
  
  # Initialize table
  tbl <- matrix("", nrow = nrow(counts)+1, ncol = ncol(counts) + 4) 
  
  # Add variable name and levels of Y to first row
  tbl[, 1] <- c(paste(yname, ", n (%)", sep = ""), paste("  ", ylevels, sep = ""))
  
  # Add N column
  tbl[1, 2] = sprintf("%.0f", sum(counts))
  
  # n (%) for each cell
  for (ii in 1:nrow(counts)) {
    yval <- rownames(counts)[ii]
    totmean <- svymean(y == yval, design = svy2)
    tbl[ii+1, 3] <- paste(sprintf("%.0f", rowSums(counts)[ii]), " (", sprintf(spf, totmean*100), ")", sep = "")
    yval <- rownames(counts)[ii]
    levmeans <- svyby(~y == yval, by = ~x, FUN = svymean, design = svy2)
    tbl[ii+1, 4:(ncol(tbl)-1)] <- paste(sprintf("%.0f", counts[ii, ]), " (", sprintf(spf, levmeans$"y == yvalTRUE"*100), ")", sep = "")
  }
  
  # Statistical test
  if (nrow(counts) == 1) {
    pval <- "-"
  } else {
    pval <- svychisq(~y + x, design = svy2, statistic = test)$p.value
  }
  tbl[1, ncol(tbl)] <- formatp(p = pval, cuts = p.cuts, decimals = p.decimals, lowerbound = p.lowerbound, leading0 = p.leading0, avoid1 = p.avoid1)
  
  #   # If y binary and compress is TRUE, compress table to a single row
  #   if (nrow(counts) <= 2 & compress == TRUE) {
  #     tbl <- matrix(c(tbl[1, 1:2], tbl[nrow(tbl), 3:(ncol(tbl)-1)], tbl[1, ncol(tbl)]), nrow = 1)
  #   }
  
  # If xlevels unspecified, set to actual values
  if (is.null(xlevels)) {
    xlevels <- colnames(counts)
  }
  
  # If y binary and compress is TRUE, compress table to a single row
  if (nrow(counts) <= 2 & compress == TRUE) {
    if (is.null(compress.val)) {
      tbl <- matrix(c(tbl[1, 1:2], tbl[nrow(tbl), 3:(ncol(tbl)-1)], tbl[1, ncol(tbl)]), nrow = 1)
    } else {
      whichrow <- which(rownames(counts) == as.character(compress.val)) + 1
      tbl <- matrix(c(tbl[1, 1:2], tbl[whichrow, 3:(ncol(tbl)-1)], tbl[1, ncol(tbl)]), nrow = 1)
    }
  }
  
  # Add column names, with sample sizes for each group if requested
  if (n.headings == FALSE) {
    colnames(tbl) <- c(variable.colname, "N", "Overall", xlevels, "P")
  } else {
    colnames(tbl) <- c(variable.colname, "N", paste(c("Overall", xlevels), " (n = ", c(sum(counts), apply(counts, 2, sum)), ")", sep = ""), "P")
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
    spacelocs <- which(substr(tbl[, variable.colname], 1, 2) == "  ")
    if (length(spacelocs) > 0) {
      tbl[spacelocs, variable.colname] <- paste("\\hskip .3cm ", substring(tbl[spacelocs, variable.colname], 3), sep = "")
    }
    chars <- strsplit(tbl[, variable.colname], "")
    for (ii in 1:length(chars)) {
      percentlocs <- which(chars[[ii]] == "%")
      if (length(percentlocs) > 0) {
        chars[[ii]][percentlocs] <- "\\%"
      }
    }
    tbl[, variable.colname] <- sapply(chars, function(x) paste(x, sep = "", collapse = ""))
    if (bold.colnames == TRUE) {
      colnames(tbl) <- paste("$\\textbf{", colnames(tbl), "}$", sep = "")
    }
    if (bold.varnames == TRUE) {
      tbl[1, 1] <- paste("$\\textbf{", tbl[1, 1], "}$")
    }
    if (bold.varlevels == TRUE) {
      tbl[2:nrow(tbl), 1] <- paste("$\\textbf{", tbl[2:nrow(tbl), 1], "}$", sep = "")
    }
  }
  
  # Return table
  return(tbl)
  
}