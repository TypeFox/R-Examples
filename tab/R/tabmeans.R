tabmeans <- function(x, y, latex = FALSE, variance = "unequal", xname = NULL, xlevels = NULL, 
                     yname = NULL, quantiles = NULL, quantile.vals = FALSE, parenth = "sd", 
                     text.label = NULL, parenth.sep = "-", decimals = NULL, p.include = TRUE, 
                     p.decimals = c(2, 3), p.cuts = 0.01, p.lowerbound = 0.001, p.leading0 = TRUE,
                     p.avoid1 = FALSE, overall.column = TRUE, n.column = FALSE, n.headings = TRUE,
                     bold.colnames = TRUE, bold.varnames = FALSE, variable.colname = "Variable", 
                     fig = FALSE, fig.errorbars = "z.ci") {
  
  # If yname/xname unspecified, use variable names
  if (is.null(yname)) {
    yname <- deparse(substitute(y))
    if (grepl("\\$", yname)) {
      yname <- strsplit(yname, "\\$")[[1]][2]
    }
  }
  if (is.null(xname)) {
    xname <- deparse(substitute(x))
    if (grepl("\\$", xname)) {
      xname <- strsplit(xname, "\\$")[[1]][2]
    }
  }
  
  # If any inputs are not correct class, return error
  if (!is.logical(latex)) {
    stop("For latex input, please enter TRUE or FALSE")
  }
  if (! variance %in% c("equal", "unequal", "ftest")) {
    stop("For variance input, please enter 'equal', 'unequal', or 'ftest'")
  }
  if (!is.character(xname)) {
    stop("For xname input, please enter character string")
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
  if (! parenth %in% c("none", "sd", "se", "t.ci", "z.ci", "none")) {
    stop("For parenth input, please enter 'none', 'sd', 'se', 't.ci', or 'z.ci'")
  }
  if (!is.null(text.label) && !is.character(text.label)) {
    stop("For text.label input, please enter a character string or just leave it unspecified. Use 'none' to request no label")
  }
  if (!is.character(parenth.sep)) {
    stop("For parenth.sep input, please enter a character string (only used if parenth is set to 't.ci' or 'z.ci'; usually '-' or ', ')")
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
  if (!is.logical(fig)) {
    stop("For fig input, please enter TRUE or FALSE")
  }
  
  # Drop missing values
  locs.complete <- which(!is.na(x) & !is.na(y))
  x <- x[locs.complete]
  y <- y[locs.complete]
  
  # Create quantiles if necessary
  if (!is.null(quantiles)) {
    x <- cut(x = x, breaks = quantile(x, probs = seq(0, 1, 1/quantiles)), include.lowest = TRUE, right = TRUE, dig.lab = 3)
  }
  
  # Get means and sample sizes
  means <- tapply(X = y, INDEX = x, FUN = mean)
  ns <- tapply(X = y, INDEX = x, FUN = length)
  
  # If decimals is unspecified, set to appropriate value
  if (is.null(decimals)) {
    max.mean <- max(abs(means))
    if (max.mean >= 1000) {
      decimals <- 0
    } else if (max.mean < 1000 & max.mean >= 10) {
      decimals <- 1
    } else if (max.mean < 10 & max.mean >= 0.1) {
      decimals <- 2
    } else if (max.mean < 0.1 & max.mean >= 0.01) {
      decimals <- 3
    } else if (max.mean < 0.01 & max.mean >= 0.001) {
      decimals <- 4
    } else if (max.mean < 0.001 & max.mean >= 0.0001) {
      decimals <- 5
    } else if (max.mean < 0.0001) {
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
  
  # Calculate p-value based on ANOVA or t-test depending on number of levels of x
  if (p.include == TRUE) {
    
    if (length(xlevels) == 2) {
      
      if (variance == "equal") {
        p <- t.test(x = y[x == xvals[1]], y = y[x == xvals[2]], var.equal = TRUE)$p.value
        message(paste("Equal variance t-test was used to compare mean ", yname, " in the two groups.", sep = ""))
      } else if (variance == "unequal") {
        p <- t.test(x = y[x == xvals[1]], y = y[x == xvals[2]], var.equal = FALSE)$p.value
        message(paste("Unequal variance t-test was used to compare mean ", yname, " in the two groups.", sep = ""))
      } else if (variance == "ftest") {
        f <- var.test(x = y[x == xvals[1]], y = y[x == xvals[2]])
        if (f$p.value < 0.05) {
          p <- t.test(x = y[x == xvals[1]], y = y[x == xvals[2]], var.equal = FALSE)$p.value
          message(paste("Unequal variance t-test was used to compare mean ", yname, " in the two groups.", sep = ""))
        } else {
          p <- t.test(x = y[x == xvals[1]], y = y[x == xvals[2]], var.equal = TRUE)$p.value
          message(paste("Equal variance t-test was used to compare mean ", yname, " in the two groups.", sep = ""))
        }
      }
      
    } else {
      
      # ANOVA
      p <- anova(lm(y ~ as.factor(x)))$"Pr(>F)"[1]
      message(paste("ANOVA was used to compare means for ", yname, sep = ""))
      
    }
    
  } else {
    p <- NA
  }
  
  if (fig == FALSE) {
    
    # Initialize table
    tbl <- matrix("", nrow = 1, ncol = length(xlevels)+4)
    
    # Figure out text.label
    if (is.null(text.label)) {
      if (parenth == "sd") {
        text.label <- ", M (SD)"
      } else if (parenth == "se") {
        text.label <- ", M (SE)"
      } else if (parenth %in% c("t.ci", "z.ci")) {
        text.label <- ", M (95% CI)"
      } else if (parenth == "none") {
        text.label <- ", M"
      }
    } else if (text.label == "none") {
      text.label <- NULL
    }
    
    # Add variable name and n column to table
    tbl[1, 1] <- paste(yname, text.label, sep = "")
    tbl[1, 2] <- sprintf("%.0f", length(y))
    
    # Add M (parenth) overall and in each x group
    if (parenth == "sd") {
      sds <- tapply(X = y, INDEX = x, FUN = sd)
      tbl[1, 3] <- paste(sprintf(spf, mean(y)), " (", sprintf(spf, sd(y)), ")", sep = "")
      tbl[1, 4:(ncol(tbl)-1)] <- paste(sprintf(spf, means), " (", sprintf(spf, sds), ")", sep = "") # Move to end
    } else if (parenth == "se") {
      ses <- tapply(X = y, INDEX = x, FUN = function(x) sd(x)/sqrt(length(x)))
      tbl[1, 3] <- paste(sprintf(spf, mean(y)), " (", sprintf(spf, sd(y)/sqrt(length(y))), ")", sep = "")
      tbl[1, 4:(ncol(tbl)-1)] <- paste(sprintf(spf, means), " (", sprintf(spf, ses), ")", sep = "")
    } else if (parenth == "t.ci") {
      ses <- tapply(X = y, INDEX = x, FUN = function(x) sd(x)/sqrt(length(x)))
      ci.lower <- means - qt(p = 0.975, df = ns - 1) * ses
      ci.upper <- means + qt(p = 0.975, df = ns - 1) * ses
      tbl[1, 3] <- paste(sprintf(spf, mean(y)), " (", sprintf(spf, mean(y) - qt(p = 0.975, df = length(y) - 1) * sd(y)/sqrt(length(y))), 
                         parenth.sep, sprintf(spf, mean(y) + qt(p = 0.975, df = length(y) - 1) * sd(y)/sqrt(length(y))), ")", sep = "")
      tbl[1, 4:(ncol(tbl)-1)] <- paste(sprintf(spf, means), " (", sprintf(spf, ci.lower), parenth.sep, sprintf(spf, ci.upper), ")", sep = "")
    } else if (parenth == "z.ci") {
      ses <- tapply(X = y, INDEX = x, FUN = function(x) sd(x)/sqrt(length(x)))
      ci.lower <- means - qnorm(0.975) * ses
      ci.upper <- means + qnorm(0.975) * ses
      tbl[1, 3] <- paste(sprintf(spf, mean(y)), " (", sprintf(spf, mean(y) - qnorm(0.975) * sd(y)/sqrt(length(y))), 
                         parenth.sep, sprintf(spf, mean(y) + qnorm(0.975) * sd(y)/sqrt(length(y))), ")", sep = "")
      tbl[1, 4:(ncol(tbl)-1)] <- paste(sprintf(spf, means), " (", sprintf(spf, ci.lower), parenth.sep, sprintf(spf, ci.upper), ")", sep = "")
    } else if (parenth == "none") {
      tbl[1, 3] <- sprintf(spf, mean(y))
      tbl[1, 4:(ncol(tbl)-1)] <- sprintf(spf, means)
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
  
  } else {
    
    if (fig.errorbars == "sd") {
      
      lowerbars <- means - tapply(X = y, INDEX = x, FUN = sd)
      upperbars <- means + tapply(X = y, INDEX = x, FUN = sd)
      ylabel <- paste(yname, " (Mean +/- 1 SD)", sep = "")
      
    } else if (fig.errorbars == "se") {
      
      lowerbars <- means - tapply(X = y, INDEX = x, FUN = function(x) sd(x)/sqrt(length(x)))
      upperbars <- means + tapply(X = y, INDEX = x, FUN = function(x) sd(x)/sqrt(length(x)))
      ylabel <- paste(yname, " (Mean +/- 1 SE)", sep = "")
      
    } else if (fig.errorbars == "t.ci") {
      
      lowerbars <- means - qt(p = 0.975, df = ns - 1) * tapply(X = y, INDEX = x, FUN = function(x) sd(x)/sqrt(length(x)))
      upperbars <- means + qt(p = 0.975, df = ns - 1) * tapply(X = y, INDEX = x, FUN = function(x) sd(x)/sqrt(length(x)))
      ylabel <- paste(yname, " (Mean +/- 95% CI)", sep = "")
      
    } else if (fig.errorbars == "z.ci") {
      
      lowerbars <- means - qnorm(p = 0.975) * tapply(X = y, INDEX = x, FUN = function(x) sd(x)/sqrt(length(x)))
      upperbars <- means + qnorm(p = 0.975) * tapply(X = y, INDEX = x, FUN = function(x) sd(x)/sqrt(length(x)))
      ylabel <- paste(yname, " (Mean +/- 95% CI)", sep = "")
      
    }
    
    if (fig.errorbars == "none") {
      
      bar.range <- max(means) - min(means)
      ylim1 <- min(means) - 0.2*bar.range
      ylim2 <- max(means) + 0.2*bar.range  
      
      tbl <- plot(x = 1:length(means), y = means, main = paste("Mean ", yname, " by ", xname, sep = ""), 
                  xlim = c(0.5, (length(xlevels)+0.5)), ylim = c(ylim1, ylim2), 
                  ylab = paste(yname, " (Mean)", sep = ""), xlab = xname, 
                  xaxt = "n", cex.lab = 1.1, pch = 19)
      axis(side = 1, at = 1:length(xlevels), labels = xlevels)
      
      
    } else {
      
      bar.range <- max(c(lowerbars, upperbars)) - min(c(lowerbars, upperbars))
      ylim1 <- min(c(lowerbars, upperbars) - 0.1*bar.range)
      ylim2 <- max(c(lowerbars, upperbars) + 0.1*bar.range)
      
      tbl <- plot(x = NULL, y = NULL, main = paste("Mean ", yname, "by ", xname), 
                  xlim = c(0.5, (length(xlevels)+0.5)), ylim = c(ylim1, ylim2), 
                  ylab = ylabel, xlab = xname, xaxt = "n", cex.lab = 1.1)
      
      for (ii in 1:length(lowerbars)) {
        
        endpoints <- c(lowerbars[ii], upperbars[ii])
        points(x = ii, y = means[ii], pch = 19)
        lines(x = rep(ii, 2), y = endpoints)
        lines(x = c((ii - 0.03), (ii + 0.03)), y = rep(endpoints[1], 2))
        lines(x = c((ii - 0.03), (ii + 0.03)), y = rep(endpoints[2], 2))
        
      }
      axis(side = 1, at = 1:length(xlevels), labels = xlevels)
      
    }
    
    tbl <- recordPlot()
    
  }
  
  # Return table
  return(tbl)
  
}