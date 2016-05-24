#' Regularise set of profiles
#'
#' This function takes a set of profiles and regularise them by means of a SVM
#'
#' @param profiles Matrix of y values, one column per profile
#' @param x Vector of predictive variable values, common to all profiles
#' @param svm.c SVM parameter (cost)
#' @param svm.eps SVM parameter (epsilon)
#' @param svm.gamma SVM parameter (gamma)
#' @param parsvm.unique Same parameters for all profiles? (logical [TRUE])
#'
#' @return Regularized profiles
#' 
#' @note The package \code{e1071} is needed in order to be able to use this function. SVM Parameters can be vectors of the same lenght as number of profiles, or a single value for all of them
#'
#' @author Javier M. Moguerza and Emilio L. Cano
#'
#' @references Cano, E.L. and Moguerza, J.M. and Prieto Corcoba, M. (2015)
#' \emph{Quality Control with R. An ISO Standards Approach}. Springer.
#' 
#' @export
#' @examples
#' wby.smooth <- smoothProfiles(profiles = ss.data.wby,
#'     x = ss.data.wbx)
#' plotProfiles(profiles = wby.smooth,
#'     x = ss.data.wbx)     
smoothProfiles <- function(profiles, x = 1:nrow(profiles), svm.c = NULL, svm.eps = NULL, svm.gamma = NULL, parsvm.unique = TRUE){
  profiles <- as.matrix(profiles)
  ncolprofiles <- ncol(profiles)
  nrowprofiles <- nrow(profiles)
  if (!is.null(svm.c) & (length(svm.c) != 1 | length(svm.c) != ncol(profiles))){
    stop("Incorrect number of svm.c parameters: there should be 1 or as many as profiles")
  }
  if (!is.null(svm.c) & (length(svm.eps) != 1 | length(svm.c) != ncol(profiles))){
    stop("Incorrect number of svm.eps parameters: there should be 1 or as many as profiles")
  }
  if (!is.null(svm.c) & (length(svm.gamma) != 1 | length(svm.c) != ncol(profiles))){
    stop("Incorrect number of svm.gamma parameters: there should be 1 or as many as profiles")
  }
  paramatrix <- matrix(0, nrow = 3, ncol = ncolprofiles)
  if (!is.null(svm.c)) {
    if (is.numeric(svm.c)){
      paramatrix[1, ] <- ifelse(parsvm.unique, median(svm.c), svm.c)
    } else{
      stop("svm.c should be numeric")
    }
  }
  if (!is.null(svm.eps)) {
    if (!is.null(svm.eps) & is.numeric(svm.eps)){
      paramatrix[2, ] <- ifelse(parsvm.unique, median(svm.eps), svm.eps)
    } else{
      stop("svm.eps should be numeric")
    }
  }
  if (!is.null(svm.gamma)) {
    if (!is.null(svm.gamma) & is.numeric(svm.gamma)){
      paramatrix[3, ] <- ifelse(parsvm.unique, median(svm.gamma), svm.gamma)
    } else{
      stop("svm.gamma should be numeric")
    }
  }
  for (i in 1:ncolprofiles){
    y <- profiles[, i]
    
    ## c SVM parameter
    if (is.null(svm.c)){
      paramatrix[1, i] <- max(c(abs(mean(y) + 3*sd(y)), abs(mean(y) - 3*sd(y))))
    }
    ## eps SVM parameter
    if (is.null(svm.eps)){
      mloess <- loess(y ~ x)
      yhat <- predict(mloess, newdata = x)
      deltas <- y - yhat
      par.sigma <- sd(deltas)
      paramatrix[2, i] <- 3*par.sigma*sqrt(log(nrowprofiles)/nrowprofiles)
    }
    if (is.null(svm.gamma)){
      ## gamma SVM parameter
      par.p <- 0.3*diff(range(x))
      paramatrix[3, i] <- 1/(2*par.p^2)
    }
  }
  if (parsvm.unique){
    parmedian <- apply(X = paramatrix, 1, median)
    paramatrix[1, ] <- parmedian[1]
    paramatrix[2, ] <- parmedian[2]
    paramatrix[3, ] <- parmedian[3]
  }
  
  reg.profiles <-  matrix(0, nrowprofiles, ncolprofiles)
  colnames(reg.profiles) <- colnames(profiles)
  
  # All regularised curves are created, along with a vector of residuals
  
  for(i in 1:ncolprofiles){
    y <- profiles[, i]
    # if (!require(e1071)){
    #   stop("Package e1071 is not installed!")
    # } else{
      
      x.svm <- e1071::svm(x, y, type = "eps-regression", 
                   cost = paramatrix[1, i], 
                   epsilon = paramatrix[2, i], 
                   gamma = paramatrix[3, i], 
                   scale = FALSE)
      reg.profiles[, i] = predict(x.svm, x)
      
    # }
  }
  return(reg.profiles)  
}


#' Compute profiles limits
#' 
#' Function to compute prototype profile and confidence bands for a set of profiles (Phase I)
#' 
#' @param profiles Matrix with profiles in columns
#' @param x Vector for the independent variable
#' @param smoothprof regularize profiles? [FALSE]
#' @param smoothlim regularize confidence bands? [FALSE]
#' @param alpha limit for control limits [0.01]
#'
#' @return a matrix with three profiles: prototype and confidence bands
#' 
#' @author Javier M. Moguerza and Emilio L. Cano
#' 
#' @references Cano, E.L. and Moguerza, J.M. and Prieto Corcoba, M. (2015)
#' \emph{Quality Control with R. An ISO Standards Approach}. Springer.
#' 
#' @export
#' @examples 
#' wby.phase1 <- ss.data.wby[, 1:35]
#' wb.limits <- climProfiles(profiles = wby.phase1,
#'     x = ss.data.wbx,
#'     smoothprof = FALSE,
#'     smoothlim = FALSE)
#'     plotProfiles(profiles = wby.phase1,
#'                  x = ss.data.wbx, 
#'                  cLimits = wb.limits)
climProfiles <- function(profiles, 
                         x = 1:nrow(profiles), 
                         smoothprof = FALSE, 
                         smoothlim = FALSE, 
                         alpha = 0.01){
  nrowProfiles <- nrow(profiles)
  # if (smoothprof | smoothlim){
  #   library(e1071)
  # }
  if (smoothprof){
    profiles <- smoothProfiles(profiles)
  }
  ucLim <- rep(0, nrowProfiles)
  lcLim <- rep(0, nrowProfiles)
  cLine <- rep(0, nrowProfiles)
  
  for(i in 1:nrowProfiles)
  {
    lcLim[i] <- quantile(profiles[i, ], alpha/2)
    ucLim[i] <- quantile(profiles[i, ], 1 - (alpha/2))
    cLine[i] <- median(profiles[i, ])
  }
  cLimits <- cbind(lcLim, ucLim, cLine)
  colnames(cLimits) <- c("LCL", "UCL", "CL")
  if (smoothlim){
    cLimits <- smoothProfiles(cLimits, x)
  }
  return(cLimits)
}


#' Plot Profiles
#' 
#' Plot profiles and optionally control limits
#' 
#' @param profiles matrix with profiles in columns
#' @param x vector with the independent variable
#' @param cLimits matrix with three profiles: prototype and confidence bands (limits)
#' @param outControl identifiers of out-of-control profiles
#' @param onlyout plot only out-of-control profiles? [FALSE]
#' 
#' @return Only graphical output with the profiles
#' 
#' @author Javier M. Moguerza and Emilio L. Cano
#' 
#' 
#' @references Cano, E.L. and Moguerza, J.M. and Prieto Corcoba, M. (2015)
#' \emph{Quality Control with R. An ISO Standards Approach}. Springer.
#' 
#' @export
#' @examples 
#' plotProfiles(profiles = ss.data.wby,
#'     x = ss.data.wbx)     
plotProfiles <- function(profiles, 
                         x = 1:nrow(profiles), 
                         cLimits = NULL, 
                         outControl = NULL, onlyout = FALSE){
  # library(scales)
  ncolProfiles <- ncol(profiles)
  plot(x, profiles[, 1], ylim = range(profiles), type = "n",
       main = "Profiles", xlab = "", ylab="", las = 1)
  if (!onlyout){
    for(i in 1:ncolProfiles){
      lines(x, profiles[, i], col = scales::alpha("black", 0.5))
    }
  }
  if (!is.null(cLimits)){
    points(x, cLimits[, 1], col = "blue", type="l", lwd = 2)
    points(x, cLimits[, 2], col = "blue", type="l", lwd = 2)
    points(x, cLimits[, 3], col = "green3", type="l", lwd = 2)
  }
  if (!is.null(outControl)){
    for (i in 1:length(outControl)){
      points(x, profiles[, outControl[i]], type = "l", col = "red4", lwd = 2)
    }
  }
}

#' Get out-of-control profiles
#' 
#' Returns a list with information about the out-of-control 
#' profiles given a set of profiles and some control limits
#' 
#' @param profiles Matrix of profiles
#' @param x Vector with the independent variable
#' @param cLimits Matrix with the prototype and confidence bands profiles
#' @param tol Tolerance (\%)
#' 
#' @return a list with the following elements:
#' \item{labOut}{labels of the out-of-control profiles}
#' \item{idOut}{ids of the out-of-control profiles}
#' \item{pOut}{proportion of times the profile values are out of the limits}
#' 
#' @references Cano, E.L. and Moguerza, J.M. and Prieto Corcoba, M. (2015)
#' \emph{Quality Control with R. An ISO Standards Approach}. Springer.
#' 
#' @export
#' @examples 
#' wby.phase1 <- ss.data.wby[, 1:35]
#' wb.limits <- climProfiles(profiles = wby.phase1,
#'     x = ss.data.wbx,
#'     smoothprof = TRUE,
#'     smoothlim = TRUE)
#' wby.phase2 <- ss.data.wby[, 36:50]
#' wb.out.phase2 <- outProfiles(profiles = wby.phase2,
#'     x = ss.data.wbx,
#'     cLimits = wb.limits,
#'     tol = 0.8)
#' wb.out.phase2
#' plotProfiles(wby.phase2,
#'     x = ss.data.wbx,
#'     cLimits = wb.limits,
#'     outControl = wb.out.phase2$idOut,
#'     onlyout = TRUE)
outProfiles <- function(profiles, 
                        x = 1:nrow(profiles), 
                        cLimits, 
                        tol = 0.5){
  ncolProfiles <- ncol(profiles)
  nrowProfiles <- nrow(profiles)
  nOut <- rep(0, ncolProfiles)
  for(i in 1:ncolProfiles)
  {
    nOut[i] <- sum((profiles[, i] >= cLimits[, 2]) | (profiles[, i] <= cLimits[, 1]))
  }
  pOut <- nOut/nrowProfiles
  idOut <- which(pOut >= tol)
  if (length(idOut) == 0){
    idOut <- NULL
  }
  labOut <- colnames(profiles)[idOut]
  if (length(labOut) == 0){
    labOut <- NULL
  }
  return(list(labOut = labOut, idOut = idOut, pOut = round(pOut, 2)))
}


#' Profiles control plot
#'
#' Plots the proportion of times that each profile remains 
#' out of the confidence bands
#'
#' @param pOut identifiers of profiles out of control
#' @param tol tolerance for the proportion of times the value of the profile is out of control
#' 
#' @return There is only graphical output
#'  
#' @references Cano, E.L. and Moguerza, J.M. and Prieto Corcoba, M. (2015)
#' \emph{Quality Control with R. An ISO Standards Approach}. Springer.
#' 
#' @author Javier M. Moguerza and Emilio L. Cano
#' 
#' 
#' @export
#' @examples 
#' wby.phase1 <- ss.data.wby[, 1:35]
#' wb.limits <- climProfiles(profiles = wby.phase1,
#'     x = ss.data.wbx,
#'     smoothprof = TRUE,
#'     smoothlim = TRUE)
#' wby.phase2 <- ss.data.wby[, 36:50]
#' wb.out.phase2 <- outProfiles(profiles = wby.phase2,
#'     x = ss.data.wbx,
#'     cLimits = wb.limits,
#'     tol = 0.8)
#' plotControlProfiles(wb.out.phase2$pOut, tol = 0.8)
plotControlProfiles <- function(pOut, tol = 0.5){
  plot(pOut, type = "b", pch = 16, xlab = "Profile", ylab = "Out-of-control rate",
       main = "Profiles control chart", las = 1)
  abline(h = tol, col = "red3", lwd = 2)
}
