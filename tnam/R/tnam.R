# This file contains the data preparation and estimation functions for temporal 
# or cross-sectional network autocorrelation models. Written by Philip Leifeld.


# display version number and date when the package is loaded
.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage(
    'Package:  tnam\n', 
    'Version:  ', desc$Version, '\n', 
    'Date:     ', desc$Date, '\n', 
    'Authors:  Philip Leifeld (Eawag and University of Bern)\n',
    '          Skyler J. Cranmer (The Ohio State University)\n'
  )
}


# function which aggregates data for glm analysis
tnamdata <- function(formula, center.y = FALSE) {
  
  # parse the formula
  if (class(formula) != "formula") {
    stop("'formula' must be a formula object.")
  }
  lhs <- deparse(formula[[2]])  # name of the response variable
  lhs <- eval(parse(text = lhs))  # get the actual response data
  rhs <- paste0(deparse(formula[[3]]), collapse = "")  # rhs of formula
  rhs <- gsub("\\s+", " ", rhs)  # get rid of redundant spaces
  rhs <- strsplit(rhs, " \\+ ")[[1]]  # parse separate formula elements
  
  # create data frame with response variable, time, and nodes
  time <- numeric()
  node <- character()
  response <- numeric()
  if (class(lhs) == "list") {
    for (i in 1:length(lhs)) {
      if (!is.numeric(lhs[[i]])) {
        stop(paste("The response variable should be numeric or a list of", 
            "numerics or a data frame with one time point per column."))
      }
      if (is.null(names(lhs[[i]])) || length(names(lhs[[i]])) != 
          length(lhs[[i]])) {
        stop(paste("The outcome variable must have node labels if multiple", 
            "time points are present."))
      }
      node <- c(node, names(lhs[[i]]))
      time <- c(time, rep(i, length(lhs[[i]])))
      if (center.y == TRUE) {
        lhs[[i]] <- lhs[[i]] - mean(lhs[[i]], na.rm = TRUE)
      }
      response <- c(response, lhs[[i]])
    }
  } else if (class(lhs) == "data.frame") {
    for (i in 1:ncol(lhs)) {
      if (!is.numeric(lhs[, i])) {
        stop(paste("The response variable should be numeric or a list of", 
            "numerics or a data frame with one time point per column."))
      }
      if (is.null(rownames(lhs)) || length(rownames(lhs)) != nrow(lhs)) {
        stop(paste("The outcome variable must have node labels if multiple", 
            "time points are present."))
      }
      node <- c(node, rownames(lhs))
      time <- c(time, rep(i, nrow(lhs)))
      if (center.y == TRUE) {
        lhs[, i] <- lhs[, i] - mean(lhs[, i], na.rm = TRUE)
      }
      response <- c(response, lhs[, i])
    }
  } else if (!is.numeric(lhs)) {
    stop("Data type of the response variable could not be recognized.")
  } else {
    response <- lhs
    if (center.y == TRUE) {
      response <- response - mean(response, na.rm = TRUE)
    }
    time <- rep(1, length(lhs))
    node <- as.character(1:length(lhs))
  }
  dat <- data.frame(response = response, time = time, node = node)
  
  # compute results according to rhs
  resultlist <- list()
  for (i in 1:length(rhs)) {
    result <- eval(parse(text = rhs[i]))
    resultlist[[i]] <- result
  }
  
  # check compatibility of labels
  for (i in 1:length(resultlist)) {
    for (j in 1:length(resultlist)) {
      itime <- length(unique(resultlist[[i]]$time))
      jtime <- length(unique(resultlist[[j]]$time))
      if ((itime > 1 || jtime > 1) && i < j) {
        inters <- length(intersect(resultlist[[i]]$node, resultlist[[j]]$node))
        if (inters == 0) {
          stop(paste("Model terms", i, "and", j, "do not have any", 
              "intersecting node labels. Please attach names, row names, or", 
              "vertex names to the 'y' or 'networks' argument."))
        }
      }
    }
  }
  
  # take care of the lags
  for (i in 1:length(resultlist)) {
    lag.i <- attributes(resultlist[[i]])$lag
    if (is.null(lag.i) || length(lag.i) == 0) {
      lag.i <- 0
    }
    resultlist[[i]]$time <- resultlist[[i]]$time + lag.i
  }
  
  # merge results with response variable and take care of lags
  for (i in 1:length(resultlist)) {
    dat <- merge(dat, resultlist[[i]], by = c("time", "node"), all.x = TRUE, 
        all.y = FALSE)
    colnames(dat)[3] <- "response"
    dat$node <- as.character(dat$node)
    if (ncol(resultlist[[i]]) == 4) {
      dat <- dat[, -ncol(dat)]
    }
  }
  dat <- dat[, c(3, 1, 2, 4:ncol(dat))]
  
  return(dat)
}


# temporal network autocorrelation model
tnam <- function(formula, family = gaussian, re.node = FALSE, 
    re.time = FALSE, time.linear = FALSE, time.quadratic = FALSE, 
    center.y = FALSE, na.action = na.omit, ...) {
  
  # prepare the data frame
  dat <- tnamdata(formula, center.y = center.y)
  
  # check if GLM is appropriate
  if (re.node == FALSE && re.time == FALSE && length(unique(dat$time)) > 1) {
    warning(paste("Different time points are available. You might want to use", 
        "a mixed effects model using arguments 're.time' and/or 're.node'."))
  }
  
  # take care of the node variable: keep as random effect or remove
  if (re.node == TRUE && length(unique(dat$time)) > 1) {
    glmest.node <- FALSE
  } else {
    dat <- dat[, -3]
    glmest.node <- TRUE
  }
  
  # take care of the time variable
  if (time.linear == TRUE && time.quadratic == TRUE && re.time == TRUE) {
    # T-T-T
    if (length(unique(dat$time)) > 1) {
      glmest.time <- FALSE
      dat$re.time <- dat$time  # add RE
      dat <- dat[, -2]  # remove linear effect
      warning(paste("Arguments 're.time' and 'time.linear' cannot be used", 
          "together. Omitting the linear time effect."))
      warning(paste("Arguments 're.time' and 'time.quadratic' cannot be used", 
          "together. Omitting the quadratic time effect."))
    } else {
      glmest.time <- TRUE
      dat <- dat[, -2]  # remove linear effect
      message("Time effects are ignored because only one time step is present.")
    }
  } else if (time.linear == TRUE && time.quadratic == TRUE && 
      re.time == FALSE) {
    # T-T-F
    glmest.time <- TRUE
    if (length(unique(dat$time)) > 1) {
      dat$time.squared <- dat$time^2
    } else {
      dat <- dat[, -2]  # remove linear effect
      message("Time effects are ignored because only one time step is present.")
    }
  } else if (time.linear == TRUE && time.quadratic == FALSE && 
      re.time == FALSE) {
    # T-F-F
    glmest.time <- TRUE
    if (length(unique(dat$time)) > 1) {
      # OK; do not modify anything
    } else {
      dat <- dat[, -2]  # remove linear effect
      message("Time effects are ignored because only one time step is present.")
    }
  } else if (time.linear == FALSE && time.quadratic == FALSE && 
      re.time == FALSE) {
    # F-F-F
    dat <- dat[, -2]  # remove linear effect
    glmest.time <- TRUE
  } else if (time.linear == FALSE && time.quadratic == TRUE && 
      re.time == FALSE) {
    # F-T-F
    glmest.time <- TRUE
    if (length(unique(dat$time)) > 1) {
      dat$time.squared <- dat$time^2  # create quadratic effect
    } else {
      message("Time effects are ignored because only one time step is present.")
    }
    dat <- dat[, -2]  # remove linear effect
  } else if (time.linear == FALSE && time.quadratic == TRUE && 
      re.time == TRUE) {
    # F-T-T
    if (length(unique(dat$time)) > 1) {
      glmest.time <- FALSE
      dat$re.time <- dat$time  # add RE
      message(paste("Arguments 're.time' and 'time.quadratic' cannot be used", 
          "together. Omitting the quadratic time effect."))
    } else {
      glmest.time <- TRUE
      message("Time effects are ignored because only one time step is present.")
    }
    dat <- dat[, -2]  # remove linear effect
  } else if (time.linear == FALSE && time.quadratic == FALSE && 
      re.time == TRUE) {
    # F-F-T
    if (length(unique(dat$time)) > 1) {
      glmest.time <- FALSE
      dat$re.time <- dat$time  # add RE
    } else {
      glmest.time <- TRUE
      message("Time effects are ignored because only one time step is present.")
    }
    dat <- dat[, -2]  # remove linear effect
  } else if (time.linear == TRUE && time.quadratic == FALSE && 
      re.time == TRUE) {
    # T-F-T
    if (length(unique(dat$time)) > 1) {
      glmest.time <- FALSE
      dat$re.time <- dat$time  # add RE
      dat <- dat[, -2]  # remove linear effect
      warning(paste("Arguments 're.time' and 'time.linear' cannot be used", 
          "together. Omitting the linear time effect."))
    } else {
      glmest.time <- TRUE
      dat <- dat[, -2]  # remove linear effect
      message("Time effects are ignored because only one time step is present.")
    }
  }
  
  if (glmest.node == FALSE || glmest.time == FALSE) {
    glmest <- FALSE
  } else {
    glmest <- TRUE
  }
  
  # estimate!
  if (glmest == TRUE) {  # GLM is necessary; no random effects
    model <- glm(dat, family = family, na.action = na.action, ...)
  } else {  # mixed-effects model (lme4) is required; random effects present
    if (is.character(family)) {
      family <- get(family, mode = "function", envir = parent.frame(2))
    } else if (is.function(family)) {
      family <- family()
    }
    if (isTRUE(all.equal(family, gaussian()))) {  # gaussian link: use lmer
      if (re.node == TRUE && re.time == TRUE) {
        model <- lme4::lmer(response ~ . - re.time - node + (1|re.time) + 
            (1|node), data = dat, na.action = na.action, ...)
      } else if (re.node == TRUE && re.time == FALSE) {
        model <- lme4::lmer(response ~ . - node + (1|node), data = dat, 
            na.action = na.action, ...)
      } else if (re.node == FALSE && re.time == TRUE) {
        model <- lme4::lmer(response ~ . -re.time + (1|re.time), data = dat, 
            na.action = na.action, ...)
      }
    } else {
      if (re.node == TRUE && re.time == TRUE) {  # other link function: glmer
        model <- lme4::glmer(response ~ . - re.time - node + (1|re.time) + 
            (1|node), data = dat, family = family, na.action = na.action, ...)
      } else if (re.node == TRUE && re.time == FALSE) {
        model <- lme4::glmer(response ~ . - node + (1|node), data = dat, 
            family = family, na.action = na.action, ...)
      } else if (re.node == FALSE && re.time == TRUE) {
        model <- lme4::glmer(response ~ . - re.time + (1|re.time), data = dat, 
            family = family, na.action = na.action, ...)
      }
    }
  }
  
  return(model)
}

