tabglm <- function(glmfit, latex = FALSE, xlabels = NULL, ci.beta = TRUE, inference = "wald", 
                   decimals = 2, p.decimals = c(2, 3), p.cuts = 0.01, p.lowerbound = 0.001, 
                   p.leading0 = TRUE, p.avoid1 = FALSE, basic.form = FALSE, intercept = TRUE, 
                   n = FALSE, events = FALSE, greek.beta = FALSE, binary.compress = TRUE, 
                   bold.colnames = TRUE, bold.varnames = FALSE, bold.varlevels = FALSE, 
                   predictor.colname = "Variable") {
  
  # If glmfit is not correct class, return error
  if (!all(class(glmfit) == c("glm", "lm"))) {
    stop("For glmfit input, please enter an object returned from the glm function")
  }
  # If any predictors are ordered factors, return error
  if (any(class(glmfit$model[, -1])[1] == "ordered") & basic.form == FALSE) {
    stop("tabglm does not work with ordered factors because dummie coding can be unpredictable. 
         Please re-run the glm after converting ordered factors to regular factors")
  }
  # If any other inputs are not correct class, return error
  if (!is.logical(latex)) {
    stop("For latex input, please enter TRUE or FALSE")
  }
  if (!is.logical(ci.beta)) {
    stop("For ci.beta input, please enter TRUE or FALSE")
  }
  if (! inference %in% c("wald", "wald.z", "profile", "profile.z")) {
    stop("For inference input, please enter 'wald', 'wald.z', 'profile', or 'profile.z'")
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
  if (!is.logical(n)) {
    stop("For n input, please enter TRUE or FALSE")
  }
  if (!is.logical(events)) {
    stop("For events input, please enter TRUE or FALSE")
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
  
  # Store glm.fit results
  coef <- summary(glmfit)$coefficients
  xnames <- rownames(coef)
  model <- glmfit$model
  
  # Assign perc variable for constructing confidence intervals
  if (inference %in% c("wald", "profile")) {
    if (colnames(summary(glmfit)$coefficients)[3] == "t value") {
      degfree <- glmfit$df.residual
      perc <- qt(0.975, degfree)
    } else {
      degfree <- Inf
      perc <- qnorm(0.975)
    }
  } else if (inference %in% c("wald.z", "profile.z")) {
    degfree <- Inf
    perc <- qnorm(0.975)
  }
  if (inference %in% c("profile", "profile.z")) {
    ci <- confint(glmfit)
  }
  
  # Initialized vectors for formatting factor variables in table
  spaces <- c()
  refs <- c()
  
  # Get indices of variable names
  predcounter <- 0
  pred <- c()
  for (ii in 2:ncol(model)) {
    pred[(ii-1)] <- predcounter + 1
    if (! class(model[, ii]) == "factor" | (class(model[, ii]) == "factor" & length(unique(model[, ii])) == 2 & binary.compress == TRUE)) {
      predcounter <- predcounter + 1
    } else {
      predcounter <- predcounter + length(unique(model[, ii])) + 1
    }
  }
  if (intercept == TRUE) {
    pred <- c(1, pred + 1)
  }
  
  # Initialize table
  tbl <- matrix("", nrow = 100, ncol = 8)
  tbl[1, 2] <- nrow(model)
  tbl[1, 3] <- sprintf("%0.0f", sum(model[, 1]))
  
  # Create index variables for table and glm coefficients
  tabindex <- 1
  coefindex <- 1
  
  # Enter intercept information if available
  if (xnames[1] == "(Intercept)" & intercept == TRUE) {
    beta <- coef[1, 1]
    se <- coef[1, 2]
    stat <- coef[1, 3]
    p <- pt(-abs(stat), degfree)*2
    tbl[1, 1] <- "Intercept"
    tbl[1, 4] <- paste(sprintf(spf, coef[1, 1]), " (", sprintf(spf, se), ")", sep = "")
    if (inference %in% c("wald", "wald.z")) {
      tbl[1, 5] <- paste("(", sprintf(spf, beta-perc*se), ", ", sprintf(spf, beta+perc*se), ")", sep = "")
    } else {
      tbl[1, 5] <- paste("(", sprintf(spf, ci[1, 1]), ", ", sprintf(spf, ci[1, 2]), ")", sep = "")
    }
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
  if ((":" %in% unlist(strsplit(rownames(coef), ""))) | basic.form == TRUE) {
    beta <- coef[2:nrow(coef), 1]
    se <- coef[2:nrow(coef), 2]
    stat <- coef[2:nrow(coef), 3]
    or <- exp(beta)
    p = pt(-abs(stat), degfree)*2
    tbl[2:nrow(coef), 1] <- rownames(coef)[-1]
    tbl[2:nrow(coef), 4] <- paste(sprintf(spf, beta), " (", sprintf(spf, se), ")", sep = "")
    if (inference %in% c("wald", "wald.z")) {
      tbl[2:nrow(coef), 5] <- paste("(", sprintf(spf, beta-perc*se), ", ", sprintf(spf, beta+perc*se), ")", sep = "")
      tbl[2:nrow(coef), 7] <- paste("(", sprintf(spf, exp(beta-perc*se)), ", ", sprintf(spf, exp(beta+perc*se)), ")", sep = "")
    } else {
      tbl[2:nrow(coef), 5] <- paste("(", sprintf(spf, ci[2:nrow(ci), 1]), ", ", sprintf(spf, ci[2:nrow(ci), 2]))
      tbl[2:nrow(coef), 7] <- paste("(", sprintf(spf, exp(ci[2:nrow(ci), 1])), ", ", sprintf(spf, exp(ci[2:nrow(ci), 2])), ")", sep = "")
    }
    tbl[2:nrow(coef), 6] <- sprintf(spf, exp(beta))
    tbl[2:nrow(coef), 8] <- formatp(p = p, cuts = p.cuts, decimals = p.decimals, lowerbound = p.lowerbound,
                                    leading0 = p.leading0, avoid1 = p.avoid1)
    tabindex <- nrow(coef)+1
    
  } else {
  
    # Otherwise format factors neatly
    for (ii in 2:ncol(model)) {
      if (class(model[, ii])[1] != "factor" | (class(model[, ii]) == "factor" & length(unique(model[, ii])) == 2 & binary.compress == TRUE)) {
        beta <- coef[coefindex, 1]
        se <- coef[coefindex, 2]
        stat <- coef[coefindex, 3]
        or <- exp(beta)
        p <- pt(-abs(stat), degfree)*2
        tbl[tabindex, 1] <- colnames(model)[ii]
        tbl[tabindex, 4] <- paste(sprintf(spf, beta), " (", sprintf(spf, se), ")", sep = "")
        if (inference %in% c("wald", "wald.z")) {
          tbl[tabindex, 5] <- paste("(", sprintf(spf, beta-perc*se), ", ", sprintf(spf, beta+perc*se), ")", sep = "")
          tbl[tabindex, 7] <- paste("(", sprintf(spf, exp(beta-perc*se)), ", ", sprintf(spf, exp(beta+perc*se)), ")", sep = "")
        } else {
          tbl[tabindex, 5] <- paste("(", sprintf(spf, ci[coefindex, 1]), ", ", sprintf(spf, ci[coefindex, 2]), ")", sep = "")
          tbl[tabindex, 7] <- paste("(", sprintf(spf, exp(ci[coefindex, 1])), ", ", sprintf(spf, exp(ci[coefindex, 2])), ")", sep = "")
        }
        tbl[tabindex, 6] <- sprintf(spf, exp(beta))
        tbl[tabindex, 8] <- formatp(p = p, cuts = p.cuts, decimals = p.decimals, lowerbound = p.lowerbound,
                                    leading0 = p.leading0, avoid1 = p.avoid1)
        tabindex <- tabindex + 1
        coefindex <- coefindex + 1
        
      } else {
        levels <- sort(unique(model[, ii]))
        if (length(levels) == 2 & binary.compress == TRUE) {
          beta <- coef[coefindex, 1]
          se <- coef[coefindex, 2]
          stat <- coef[coefindex, 3]
          or <- exp(beta)
          p <- pt(-abs(stat), degfree)*2
          tbl[tabindex, 1] <- xnames[coefindex]
          tbl[tabindex, 4] <- paste(sprintf(spf, beta), " (", sprintf(spf, se), ")", sep = "")
          if (inference %in% c("wald", "wald.z")) {
            tbl[tabindex, 5] <- paste("(", sprintf(spf, beta-perc*se), ", ", sprintf(spf, beta+perc*se), ")", sep = "")
            tbl[tabindex, 7] <- paste("(", sprintf(spf, exp(beta-perc*se)), ", ", sprintf(spf, exp(beta+perc*se)), ")", sep = "")
          } else {
            tbl[tabindex, 5] <- paste("(", sprintf(spf, ci[coefindex, 1]), ", ", sprintf(spf, ci[coefindex, 2]), ")", sep = "")
            tbl[tabindex, 7] <- paste("(", sprintf(spf, exp(ci[coefindex, 1])), ", ", sprintf(spf, exp(ci[coefindex, 2])), ")", sep = "")
          }
          tbl[tabindex, 6] <- sprintf(spf, exp(beta))
          tbl[tabindex, 8] <- formatp(p = p, cuts = p.cuts, decimals = p.decimals, lowerbound = p.lowerbound,
                                      leading0 = p.leading0, avoid1 = p.avoid1)
          tabindex <- tabindex + 1
          coefindex <- coefindex + 1
        } else {
          tbl[tabindex, 1] <- colnames(model)[ii]
          tabindex <- tabindex + 1
          tbl[tabindex, 1] <- paste("  ", levels[1], " (ref)", sep = "")
          tbl[tabindex, 4:8] <- "-"
          spaces <- c(spaces, tabindex)
          refs <- c(refs, tabindex)
          tabindex <- tabindex + 1
          for (jj in 2:length(levels)) {
            beta <- coef[coefindex, 1]
            se <- coef[coefindex, 2]
            stat <- coef[coefindex, 3]
            or <- exp(beta)
            p <- pt(-abs(stat), degfree)*2
            tbl[tabindex, 1] <- paste("  ", levels[jj], sep = "")
            tbl[tabindex, 4] <- paste(sprintf(spf, beta), " (", sprintf(spf, se), ")", sep = "")
            if (inference %in% c("wald", "wald.z")) {
              tbl[tabindex, 5] <- paste("(", sprintf(spf, beta-perc*se), ", ", sprintf(spf, beta+perc*se), ")", sep = "")
              tbl[tabindex, 7] <- paste("(", sprintf(spf, exp(beta-perc*se)), ", ", sprintf(spf, exp(beta+perc*se)), ")", sep = "")
            } else {
              tbl[tabindex, 5] <- paste("(", sprintf(spf, ci[coefindex, 1]), ", ", sprintf(spf, ci[coefindex, 2]), ")", sep = "")
              tbl[tabindex, 7] <- paste("(", sprintf(spf, exp(ci[coefindex, 1])), ", ", sprintf(spf, exp(ci[coefindex, 2])), ")", sep = "")
            }
            tbl[tabindex, 6] <- sprintf(spf, exp(beta))
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
  colnames(tbl) <- c(predictor.colname, "N", "Events", "Beta (SE)", "95% CI for Beta", "OR", "95% CI for OR", "P")

  # If not a binary response or events is FALSE, remove events column
  if (glmfit$family$family != "binomial" | events == FALSE) {
    tbl <- tbl[, -which(colnames(tbl) == "Events"), drop = FALSE]
  }
  
  # If n is FALSE, remove column
  if (n == FALSE) {
    tbl <- tbl[, -which(colnames(tbl) == "N"), drop = FALSE]
  }

  # If ci.beta is FALSE, remove column
  if (ci.beta == FALSE) {
    tbl <- tbl[, -which(colnames(tbl) == "95% CI for Beta"), drop = FALSE]
  }
  
  # Adjust OR columns if necessary
  if (glmfit$family$link == "log") {
    colnames(tbl)[colnames(tbl) == "OR"] <- "exp(Beta)"
    colnames(tbl)[colnames(tbl) == "95% CI for OR"] <- "95% CI for exp(Beta)"
  } else if (!(glmfit$family$link == "logit" & glmfit$family$family %in% c("binomial", "quasi", "quasibibinomial"))) {
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