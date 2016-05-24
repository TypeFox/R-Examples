tabgee <- function(geefit, latex = FALSE, xlabels = NULL, ci.beta = TRUE, decimals = 2, 
                   p.decimals = c(2, 3), p.cuts = 0.01, p.lowerbound = 0.001, p.leading0 = TRUE, 
                   p.avoid1 = FALSE, basic.form = FALSE, intercept = TRUE, n.id = FALSE, 
                   n.total = FALSE, or = TRUE, robust = TRUE, data = NULL, greek.beta = FALSE,
                   binary.compress = TRUE, bold.colnames = TRUE, bold.varnames = FALSE, 
                   bold.varlevels = FALSE, predictor.colname = "Variable") {
  
  # If any inputs are not correct class, return error
  if (!all(class(geefit) == c("gee", "glm"))) {
    stop("For geefit input, please enter an object returned from the gee function")
  }
  if (!is.logical(latex)) {
    stop("For latex input, please enter TRUE or FALSE")
  }
  if (!is.logical(ci.beta)) {
    stop("For ci.beta input, please enter TRUE or FALSE")
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
  if (!is.logical(basic.form)) {
    stop("For basic.form input, please enter TRUE or FALSE")
  }
  if (!is.logical(intercept)) {
    stop("For intercept input, please enter TRUE or FALSE")
  }
  if (!is.logical(n.id)) {
    stop("For n.id input, please enter TRUE or FALSE")
  }
  if (!is.logical(n.total)) {
    stop("For n.total input, please enter TRUE or FALSE")
  }
  if (!is.logical(or)) {
    stop("For or input, please enter TRUE or FALSE")
  }
  if (!is.logical(robust)) {
    stop("For robust input, please enter TRUE or FALSE")
  }
  if (!is.null(data)) {
    if (!(is.data.frame(data) | is.matrix(data))) {
      stop("For data input, please enter data frame or matrix.")
    }
  }
  if (!is.logical(greek.beta)) {
    stop("For greek.beta input, please enter TRUE or FALSE")
  }
  if (!is.logical(binary.compress)) {
    stop("For binary.compress input, please enter TRUE or FALSE")
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
  if (!is.character(predictor.colname)) {
    stop("For predictor.colname input, please enter a character string")
  }
  
  # Convert decimals to variable for sprintf
  spf <- paste("%0.", decimals, "f", sep = "")
  
  # Store gee.fit results
  coef <- summary(geefit)$coefficients
  xnames <- rownames(coef)
  model <- attr(geefit$terms, "dataClasses")
  
  # Column for SE's and Z's based on robust input
  secol <- ifelse(robust == FALSE, 2, 4)
  zcol <- ifelse(robust == FALSE, 3, 5)
  
  # Initialized vectors for formatting factor variables in table
  spaces <- c()
  refs <- c()
  
  # Get indices of variable names
  predcounter <- 0
  pred <- c()
  for (ii in 1:(length(model)-1)) {
    pred[ii] <- predcounter + 1
    if (! model[ii] == "factor" | (model[ii] == "factor" & length(unique(data[, names(model)[ii]])) == 2 & binary.compress == TRUE)) {
      predcounter <- predcounter + 1
    } else {
      predcounter <- predcounter + length(unique(data[, names(model)[ii]]))
    }
  }
  
  # Initialize table
  tbl <- matrix("", nrow = 100, ncol = 8)
  tbl[1, 2] <- length(unique(geefit$id))
  tbl[1, 3] <- sprintf("%0.0f", geefit$nobs)
  
  # Create index variables for table and gee coefficients
  tabindex <- 1
  coefindex <- 1
  
  # Enter intercept information if available
  if (xnames[1] == "(Intercept)" & intercept == TRUE) {
    beta <- coef[1, 1]
    se <- coef[1, secol]
    z <- coef[1, zcol]
    p <- pnorm(-abs(z))*2
    tbl[1, 1] <- "Intercept"
    tbl[1, 4] <- paste(sprintf(spf, coef[1, 1]), " (", sprintf(spf, se), ")", sep = "")
    tbl[1, 5] <- paste("(", sprintf(spf, beta-1.96*se), ", ",sprintf(spf, beta+1.96*se), ")", sep = "")
    tbl[1, 6:7] <- "-"
    tbl[1, 8] <- formatp(p = p, cuts = p.cuts, decimals = p.decimals, lowerbound = p.lowerbound,
                         leading0 = p.leading0, avoid1 = p.avoid1)
    tabindex <- tabindex + 1
    coefindex <- coefindex + 1
    
  } else {
    coefindex <- coefindex + 1
  }
  
  # If there are one or more interaction terms OR basic.form is TRUE, just do basic formatting straight
  # from the table of coefficients
  if ((":" %in% unlist(strsplit(rownames(coef), ""))) | basic.form == TRUE | (sum(model == "factor") > 0 & is.null(data))) {
    beta <- coef[2:nrow(coef), 1]
    se <- coef[2:nrow(coef), secol]
    or <- exp(beta)
    z <- coef[tabindex:nrow(coef), zcol]
    p <- pnorm(-abs(z))*2
    tbl[2:nrow(coef), 1] <- rownames(coef)[-1]
    tbl[2:nrow(coef), 4] <- paste(sprintf(spf, beta), " (", sprintf(spf, se), ")", sep = "")
    tbl[2:nrow(coef), 5] <- paste("(", sprintf(spf, beta-1.96*se), ", ", sprintf(spf, beta+1.96*se), ")", sep = "")
    tbl[2:nrow(coef), 6] <- sprintf(spf, exp(beta))
    tbl[2:nrow(coef), 7] <- paste("(", sprintf(spf, exp(beta-1.96*se)), ", ", sprintf(spf, exp(beta+1.96*se)), ")", sep = "")
    tbl[tabindex:nrow(coef), 8] <- formatp(p = p, cuts = p.cuts, decimals = p.decimals, lowerbound = p.lowerbound,
                                           leading0 = p.leading0, avoid1 = p.avoid1)
    tabindex <- nrow(coef)+1
  } else {
  
    # Otherwise format factors neatly
    for (ii in 2:(length(model)-1)) {
      if (model[ii] != "factor" | (model[ii] == "factor" & length(unique(data[, names(model)[ii]])) == 2 & binary.compress == TRUE))  {
        beta <- coef[coefindex, 1]
        se <- coef[coefindex, secol]
        or <- exp(beta)
        z <- coef[coefindex, zcol]
        p <- pnorm(-abs(z))*2
        tbl[tabindex, 1] <- names(model)[ii]
        tbl[tabindex, 4] <- paste(sprintf(spf, beta), " (", sprintf(spf, se), ")", sep = "")
        tbl[tabindex, 5] <- paste("(", sprintf(spf, beta-1.96*se), ", ", sprintf(spf, beta+1.96*se), ")", sep = "")
        tbl[tabindex, 6] <- sprintf(spf, exp(beta))
        tbl[tabindex, 7] <- paste("(", sprintf(spf, exp(beta-1.96*se)), ", ", sprintf(spf, exp(beta+1.96*se)), ")", sep = "")
        tbl[tabindex, 8] <- formatp(p = p, cuts = p.cuts, decimals = p.decimals, lowerbound = p.lowerbound,
                                    leading0 = p.leading0, avoid1 = p.avoid1)
        tabindex <- tabindex + 1
        coefindex <- coefindex + 1
        
      } else {
        levels <- sort(unique(data[, names(model)[ii]]))
        if (length(levels) == 2 & binary.compress == TRUE) {
          beta <- coef[coefindex, 1]
          se <- coef[coefindex, secol]
          or <- exp(beta)
          z <- coef[coefindex, zcol]
          p <- pnorm(-abs(z))*2
          tbl[tabindex, 1] <- xnames[coefindex]
          tbl[tabindex, 4] <- paste(sprintf(spf, beta), " (", sprintf(spf, se), ")", sep = "")
          tbl[tabindex, 5] <- paste("(", sprintf(spf, beta-1.96*se), ", ", sprintf(spf, beta+1.96*se), ")", sep = "")
          tbl[tabindex, 6] <- sprintf(spf, exp(beta))
          tbl[tabindex, 7] <- paste("(", sprintf(spf, exp(beta-1.96*se)), ", ", sprintf(spf, exp(beta+1.96*se)), ")", sep = "")
          tbl[tabindex, 8] <- formatp(p = p, cuts = p.cuts, decimals = p.decimals, lowerbound = p.lowerbound,
                                      leading0 = p.leading0, avoid1 = p.avoid1)
          tabindex <- tabindex + 1
          coefindex <- coefindex + 1
        } else {
          tbl[tabindex, 1] <- names(model)[ii]
          tabindex <- tabindex + 1
          tbl[tabindex, 1] <- paste("  ", levels[1], " (ref)", sep = "")
          tbl[tabindex, 4:8] <- "-"
          spaces <- c(spaces, tabindex)
          refs <- c(refs, tabindex)
          tabindex <- tabindex + 1
          
          for (jj in 2:length(levels)) {
            tbl[tabindex, 1] <- paste("  ", levels[jj], sep = "")
            beta <- coef[coefindex, 1]
            se <- coef[coefindex, secol]
            or <- exp(beta)
            z <- coef[coefindex, zcol]
            p <- pnorm(-abs(z))*2
            tbl[tabindex, 4] <- paste(sprintf(spf, beta), " (", sprintf(spf, se), ")", sep = "")
            tbl[tabindex, 5] <- paste("(", sprintf(spf, beta-1.96*se), ", ", sprintf(spf, beta+1.96*se), ")", sep = "")
            tbl[tabindex, 6] <- sprintf(spf, exp(beta))
            tbl[tabindex, 7] <- paste("(", sprintf(spf, exp(beta-1.96*se)), ", ", sprintf(spf, exp(beta+1.96*se)), ")", sep = "")
            tbl[tabindex, 8] <- formatp(p = p, cuts = p.cuts, decimals = p.decimals, lowerbound = p.lowerbound,
                                        leading0 = p.leading0, avoid1 = p.avoid1)
            spaces <- c(spaces, tabindex)
            tabindex <- tabindex + 1
            coefindex <- coefindex + 1
          }  
        }
      }
    }
  }
  
  # Truncate table at correct number of rows
  tbl <- tbl[1:(tabindex-1),, drop = FALSE]
  
  # Add column names
  colnames(tbl) <- c(predictor.colname, "Clusters", "Observations", "Beta (SE)", "95% CI for Beta", "OR", "95% CI for OR", "P")
  
  # Remove n columns if requested
  if (n.id == FALSE) {
    tbl <- tbl[, -which(colnames(tbl) == "Clusters"), drop = FALSE]
  }
  if (n.total == FALSE) {
    tbl <- tbl[, -which(colnames(tbl) == "Observations"), drop = FALSE]
  }
  
  # If ci.beta is FALSE, remove column
  if (ci.beta == FALSE) {
    tbl <- tbl[, -which(colnames(tbl) == "95% CI for Beta"), drop = FALSE]
  }
  
  # Adjust OR columns if necessary
  if (geefit$family$link == "log") {
    colnames(tbl)[colnames(tbl) == "OR"] <- "exp(Beta)"
    colnames(tbl)[colnames(tbl) == "95% CI for OR"] <- "95% CI for exp(Beta)"
  } else if (!(geefit$family$link == "logit" & geefit$family$family %in% c("binomial", "quasi", "quasibibinomial"))) {
    tbl <- tbl[, -which(colnames(tbl) %in% c("OR", "95% CI for OR")), drop = FALSE]
  }
  
  # Add variable labels if possible
  if (!is.null(xlabels)) {
    xlabels[spaces] <- paste("  ", xlabels[spaces], sep = "")
    xlabels[refs] <- paste(xlabels[refs], "(ref)")
    tbl[1:nrow(tbl), 1] <- xlabels
  }
  
  # If latex is TRUE, do some re-formatting
  if (latex == TRUE) {
    if (greek.beta == TRUE) {
      colnames(tbl)[which(colnames(tbl) == "Beta (SE)")] <- "$\\hat{\\beta}$ (SE)"
      colnames(tbl)[which(colnames(tbl) == "95% CI for Beta")] <- "95% CI for $\\beta$"
      colnames(tbl)[which(colnames(tbl) == "exp(Beta)")] <- "exp($\\beta)$"
      colnames(tbl)[which(colnames(tbl) == "95% CI for exp(Beta)")] <- "95\\% CI for exp($\\beta$)"
    }
    plocs <- which(substr(tbl[, "P"], 1, 1) == "<")
    if (length(plocs) > 0) {
      tbl[plocs, "P"] <- paste("$<$", substring(tbl[plocs, "P"], 2), sep = "")
    }
    spacelocs <- which(substr(tbl[, predictor.colname], 1, 2) == "  ")
    if (length(spacelocs) > 0) {
      tbl[spacelocs, predictor.colname] <- paste("\\hskip .4cm ", substring(tbl[spacelocs, predictor.colname], 3), sep = "")
    }
    chars <- strsplit(colnames(tbl), "")
    for (ii in 1:length(chars)) {
      percentlocs <- which(chars[[ii]] == "%")
      if (length(percentlocs) > 0) {
        chars[[ii]][percentlocs] <- "\\%"
      }
    }
    colnames(tbl) <- sapply(chars, function(x) paste(x, sep = "", collapse = ""))
    if (bold.colnames == TRUE) {
      colnames(tbl) <- paste("$\\textbf{", colnames(tbl), "}$", sep = "")
    }
    if (bold.varnames == TRUE) {
      tbl[pred, 1] <- paste("$\\textbf{", tbl[pred, 1], "}$")
    }
    if (bold.varlevels == TRUE) {
      tbl[c(1:nrow(tbl))[! c(1:nrow(tbl)) %in% pred], 1] <- paste("$\\textbf{", tbl[c(1:nrow(tbl))[! c(1:nrow(tbl)) %in% pred], 1], "}$", sep = "")
    }
  }
  
  # Return tbl
  return(tbl)
  
}