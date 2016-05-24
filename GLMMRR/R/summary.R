#' Summarizing GLMMRR Fits for fixed-effect models
#'
#' @param object
#' an object of class RRglm.
#' @param printResiduals
#' print scaled Pearson residuals (default: false).
#' @param limitRRparameters
#' set limit for list of Randomized Response parameters (default: 10).
#' @param digits
#' minimal number of \emph{significant digits} (default: 5).
#' @param ...
#' further arguments passed to or from other methods.
#'
#' @method summary RRglm
#' @export
summary.RRglm <- function(object, printResiduals = FALSE, limitRRparameters = 10, digits = 5, ...)
{
  # Grab summary.glm output
  output.glm <- summary.glm(object, ...)

  # Create a dataframe of RR paramters and glm data to work with
  df.work <- data.frame("RRmodel" = object$RRmodel, "p1" = object$RRp1, "p2" = object$RRp2, "residuals" = as.numeric(na.omit(resid(object, type = "dev"))))

  # Combine p1 and p2 into a single column
  df.work$level <- paste("(", df.work$p1, " | ", df.work$p2, ")", sep = "")

  # Create a subset for each RR model
  rrmodels <- levels(as.factor(df.work$RRmodel))
  ls.perModelData <- list()
  for(ii in 1:length(rrmodels))
  {
    ls.perModelData[[ii]] = df.work[df.work$RRmodel == rrmodels[ii], ]
  }

  # For each RR model, which unique combinations of p1 | p2 are used
  ls.perModelp1p2 <- list()
  for (ii in 1:length(ls.perModelData))
  {
    ls.perModelp1p2[[ii]] <- levels(as.factor(ls.perModelData[[ii]]$level))
  }

  # Generate output data frames per RR model
  ls.perModelOutput <- list()
  for(ii in 1:length(rrmodels))
  {
    # Limit number of columns if desired
    # One column for each level of p1 | p2
    if(length(ls.perModelp1p2[[ii]]) < limitRRparameters)
      print.limit <- length(ls.perModelp1p2[[ii]])
    else
      print.limit <- limitRRparameters

    if (printResiduals)
    {
      # Obtain information about residuals for the model
      tmp <- summary(ls.perModelData[[ii]]$residuals)

      # Create a data frame for each RR model containing the desired information by level
      ls.perModelOutput[ii] <- data.frame("Total" = c(tmp[1], tmp[2], tmp[3], tmp[5], tmp[6]))
      for(jj in 1:print.limit)
      {
        # The subset for this model
        tmp.modelData <- ls.perModelData[[ii]]

        # Obtain information about residuals for the model and the given level
        tmp <- summary(tmp.modelData[tmp.modelData$level == ls.perModelp1p2[[ii]][[jj]], ]$residuals)

        # Attach information to the data frame
        ls.perModelOutput[[ii]] <- cbind(ls.perModelOutput[[ii]], e = c(tmp[1], tmp[2], tmp[3], tmp[5], tmp[6]))
      }

      # Set column names
      colnames(ls.perModelOutput[[ii]]) <- c("Total", ls.perModelp1p2[[ii]])
    }
  }

  # Print the output given by summary.glm()
  print(output.glm)

  cat("\n")
  cat("### GLMMRR - Binary Randomized Response Data ###")
  cat("\n")
  cat("Generalized linear fixed-effects model")
  cat("\n\n")
  cat("Family:\t\t\t", object$family$family, "\n")
  cat("Link function:\t\t", object$family$link, "\n")
  cat("Model(s):\t\t", rrmodels[1], ls.perModelp1p2[[1]], "\n")
  if(length(rrmodels) > 1)
  {
    for (ii in 2:length(rrmodels))
    {
      cat("\t\t\t", rrmodels[ii], ls.perModelp1p2[[ii]], "\n")
    }
  }
  if(printResiduals)
  {
    for (ii in 1:length(rrmodels))
    {
      cat("\n")
      cat("---------------------------------------------------------")
      cat("\n")
      cat("Deviance residuals for", rrmodels[ii], "(p1 | p2) \n")
      print(ls.perModelOutput[[ii]], digits = digits)
      cat("\n")
    }
  }
  cat("\n")
}

#' Summarizing GLMMRR Fits for mixed-effect models
#'
#' @param object
#' an object of class RRglmerMod.
#' @param printResiduals
#' print scaled Pearson residuals (default: false).
#' @param limitRRparameters
#' set limit for list of Randomized Response parameters (default: 10).
#' @param digits
#' minimal number of \emph{significant digits} (default: 5).
#' @param ...
#' further arguments passed to or from other methods.
#'
#' @method summary RRglmerMod
#' @export
summary.RRglmerMod <- function(object, printResiduals = FALSE, limitRRparameters = 10, digits = 5, ...)
{
  # Grab summary.merMod output
  output.merMod <- summary(as(object, "glmerMod"), ...)

  # Create a dataframe of RR paramters
  df.work <- data.frame("RRmodel" = object@RRparam$RRmodel, "p1" = object@RRparam$p1, "p2" = object@RRparam$p2, "residuals" = as.numeric(na.omit(resid(object, type = "pearson", scale = TRUE))))

  # Combine p1 and p2 into a single column
  df.work$level <- paste("(", df.work$p1, " | ", df.work$p2, ")", sep = "")

  # Create a subset for each RR model
  rrmodels <- levels(as.factor(df.work$RRmodel))
  ls.perModelData <- list()
  for(ii in 1:length(rrmodels))
  {
    ls.perModelData[[ii]] = df.work[df.work$RRmodel == rrmodels[ii], ]
  }

  # For each RR model, which unique combinations of p1 | p2 are used
  ls.perModelp1p2 <- list()
  for (ii in 1:length(ls.perModelData))
  {
    ls.perModelp1p2[[ii]] <- levels(as.factor(ls.perModelData[[ii]]$level))
  }

  # Generate output data frames per RR model
  ls.perModelOutput <- list()
  for(ii in 1:length(rrmodels))
  {
    # Limit number of columns if desired
    # One column for each level of p1 | p2
    if(length(ls.perModelp1p2[[ii]]) < limitRRparameters)
      print.limit <- length(ls.perModelp1p2[[ii]])
    else
      print.limit <- limitRRparameters

    if (printResiduals)
    {
      # Obtain information about residuals for the model
      tmp <- summary(ls.perModelData[[ii]]$residuals)

      # Create a data frame for each RR model containing the desired information by level
      ls.perModelOutput[ii] <- data.frame("Total" = c(tmp[1], tmp[2], tmp[3], tmp[5], tmp[6]))
      for(jj in 1:print.limit)
      {
        # The subset for this model
        tmp.modelData <- ls.perModelData[[ii]]

        # Obtain information about residuals for the model and the given level
        tmp <- summary(tmp.modelData[tmp.modelData$level == ls.perModelp1p2[[ii]][[jj]], ]$residuals)

        # Attach information to the data frame
        ls.perModelOutput[[ii]] <- cbind(ls.perModelOutput[[ii]], e = c(tmp[1], tmp[2], tmp[3], tmp[5], tmp[6]))
      }

      # Set column names
      colnames(ls.perModelOutput[[ii]]) <- c("Total", ls.perModelp1p2[[ii]])
    }
  }

  print(output.merMod)

  cat("\n")
  cat("### GLMMRR - Binary Randomized Response Data ###")
  cat("\n")
  cat("Generalized linear mixed-effects model")
  cat("\n\n")
  cat("Model(s):\t", rrmodels[1], ls.perModelp1p2[[1]], "\n")
  if(length(rrmodels) > 1)
  {
    for (ii in 2:length(rrmodels))
    {
      cat("\t\t", rrmodels[ii], ls.perModelp1p2[[ii]], "\n")
    }
  }
  if(printResiduals)
  {
    for (ii in 1:length(rrmodels))
    {
      cat("\n")
      cat("---------------------------------------------------------")
      cat("\n")
      cat("Scaled Pearson residuals for", rrmodels[ii], "(p1 | p2) \n")
      print(ls.perModelOutput[[ii]], digits = digits)
      cat("\n")
    }
  }
  cat("\n")
}
