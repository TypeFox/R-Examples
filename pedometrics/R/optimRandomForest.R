#' Optimum number of iterations to de-bias a random forest regression
#' 
#' Compute the optimum number of iterations needed to de-bias a random forest 
#' regression.
#' 
#' @param x Data frame or matrix of covariates (predictor variables).
#' 
#' @param y Numeric vector with the response variable.
#' 
#' @param niter Number of iterations. Defaults to \code{niter = 10}.
#' 
#' @param nruns Number of simulations to be used in each iteration. Defaults to
#' \code{nruns = 100}.
#' 
#' @param ntree Number of trees to grow. Defaults to \code{ntree = 500}.
#' 
#' @param ntrain Number (or proportion) of observation to be used as training 
#' cases. Defaults to 2/3 of the total number of observations.
#' 
#' @param nodesize Minimum size of terminal nodes. Defaults to 
#' \code{nodesize = 5}.
#' 
#' @param mtry Number of variables randomly sampled as candidates at each 
#' split. Defaults to 1/3 of the total number of covariates.
#' 
#' @param profile Should the profile of the standardized mean squared prediction
#' error be plotted at the end of the optimization? Defaults to 
#' \code{profile = TRUE}.
#' 
#' @param progress Should a progress bar be displayed. Defaults to 
#' \code{progress = TRUE}.
#' 
#' @details 
#' A fixed proportion of the total number of observations is used to calibrate
#' (train) the random forest regression. The set of calibration observations is
#' randomly selected from the full set of observations in each simulation. The
#' remaining observations are used as test cases (validation). In general, the 
#' smaller the calibration dataset, the more simulation runs are needed to 
#' obtain stable estimates of the mean squared prediction error (MSPE).
#' 
#' The optimum number of iterations needed to de-bias the random forest 
#' regression is obtained observing the evolution of the MSPE as the number of
#' iterations increases. The MSPE is defined as the mean of the squared 
#' differences between predicted and observed values.
#' 
#' @seealso \code{\link[randomForest]{randomForest}}
#' 
#' @references 
#' Breiman, L. Random forests. \emph{Machine Learning}. v. 45, p. 5-32, 2001.
#' 
#' Breiman, L. \emph{Using adaptive bagging to debias regressions}. Berkeley: 
#' University of California, p. 16, 1999.
#' 
#' Liaw, A. & Wiener, M. Classification and regression by randomForest. 
#' \emph{R News}. v. 2/3, p. 18-22, 2002.
#' 
#' Xu, R. \emph{Improvements to random forest methodology}. Ames, Iowa: Iowa 
#' State University, p. 87, 2013.
#' 
#' @author Ruo Xu \email{xuruo.isu@@gmail.com}, with improvements by
#' Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' 
#' @note The original function was published as part of the dissertation of 
#' Ruo Xu, which was developed under the supervision of Daniel S Nettleton 
#' \email{dnett@@iastate.edu} and Daniel J Nordman 
#' \email{dnordman@@iastate.edu}.
#' 
#' @importFrom stats predict
#' @export
# FUNCTION - OPTIMIZATION ######################################################
optimRandomForest <-
  function (x, y, niter = 10, nruns = 100, ntree = 500, ntrain = 2/3, 
            nodesize = 5, mtry = max(floor(ncol(x) / 3), 1), profile = TRUE,
            progress = TRUE) {
    
    # Check if suggested packages are installed
    pkg <- c("randomForest", "utils")
    id <- !sapply(pkg, requireNamespace, quietly = TRUE)
    if (any(id)) {
      pkg <- paste(pkg[which(id)], collapse = " ")
      stop(paste("Package(s) needed for this function to work but not",
                 "installed: ", pkg, sep = ""), call. = FALSE)
    }
    
    # Settings
    nsample <- length(y)
    if (ntrain < 1) ntrain <- round(nsample * ntrain)
    ntest <- nsample - ntrain
    
    # Start simulation
    if (progress) {
      pb <- utils::txtProgressBar(min = 1, max = niter, style = 3) 
    }
    mse <- NULL
    for (k in 1:nruns) {
      id.train <- sample(1:nsample, ntrain, replace = FALSE)
      id.test <- c(1:nsample)[-id.train]
      
      res <- .iterRandomForest(
        xtrain = x[id.train, ], ytrain = y[id.train], xtest = x[id.test, ], 
        ntree = ntree, mtry = mtry, nodesize = nodesize, niter = niter)
      
      # prediction for test cases of this holding-out
      res <- res[, -c(1:ntrain)]
      
      mse <- rbind(mse, apply(((t(res) - y[id.test]) ^ 2), 2, mean))
      
      if (progress) utils::setTxtProgressBar(pb, k)
    }
    if (progress) close(pb)
    
    # Prepare output
    colnames(mse) <- paste("iter-", 1:niter, sep = "")
    res <- list(
      mse = data.frame(
        mean = apply(mse, 2, mean),  sd = apply(mse, 2, stats::sd)),
      call = data.frame(nruns = nruns, ntree = ntree,  ntrain = ntrain, 
                        ntest = ntest, nodesize = nodesize, mtry = mtry))
    
    # Plot mse profile
    if (profile) {
      .profRandomForest(mse = mse, nruns = nruns, niter = niter)
    }
    
    # Output
    return (res)
  }
# INTERNAL FUNCTION - ITERATIONS ###############################################
.iterRandomForest <- 
  function (xtrain, ytrain, xtest, ntree, nodesize, mtry, niter) {
    # xtrain is the design matrix for training cases, n*p; 
    # ytrain is the response for training cases, a vector of length n;
    # xtest is the design matrix for test cases, m*p;
    # ntree is the number of trees per forest;
    # nodesize is the maximal node size per tree;
    # niter is the number of iterations for the bias-correction RFs.
    
    # Initial settings
    ni <- 0
    ntrain <- nrow(xtrain)
    ntest <- nrow(xtest)
    pred.iter <- NULL # This records the predicted value for all X's
    b <- ytrain # Reponse of training data of each iteration
    
    repeat {
      ni <- ni + 1
      RF.temp <- randomForest::randomForest(
        x = xtrain, y = b, ntree = ntree, mtry = mtry, nodesize = nodesize)
      bpred.oob <- RF.temp$predicted
      
      # Predict for test cases
      bpred.test <- as.numeric(predict(object = RF.temp, newdata = xtest))
      pred.iter <- rbind(pred.iter, c(bpred.oob, bpred.test))
      b <- b - bpred.oob
      if (ni >= niter) break
    } # repeat ends
    
    rownames(pred.iter) <- paste("iter-", 1:niter, sep = "")
    colnames(pred.iter) <- c(
      paste("Tr-", 1:ntrain, sep = ""),  paste("Test-", 1:ntest, sep = ""))
    
    # This saves the final prediction results of different iterations of BC
    pred <- pred.iter[1,]
    for (k in 2:niter) {
      pred <- rbind(pred, apply(pred.iter[1:k,], 2, sum))
    } # for k ends
    
    rownames(pred) <- paste("iter-", 1:niter, sep = "")
    colnames(pred) <- c(
      paste("Tr-", 1:ntrain, sep = ""), paste("Test-", 1:ntest, sep = ""))
    
    return(pred)
    
  }
# INTERNAL FUNCTION - PLOT MSE PROFILE #########################################
.profRandomForest <-
  function (mse, nruns, niter) {
    
    # Prepare data for plotting
    eb <- apply(mse, 2, stats::sd) / sqrt(nruns)
    mse <- apply(mse, 2, mean)
    upper <- (mse + eb) / mse[1]
    lower <- (mse - eb) / mse[1]
    mse <- mse / mse[1]
    
    # Plotting
    graphics::plot(
      1:niter, mse, type = "b", ylim = c(min(lower), max(upper)),
      ylab = "Standardized mean squared error", xlab = "Iteration",
      main = "Mean squared error profile")
    graphics::abline(h = mse[1], col = "red")
    graphics::arrows(
      1:niter, lower, 1:niter, upper, length = 0.05, angle = 90, code = 3)
    graphics::points(x = which.min(mse), y = min(mse), col = "blue", lwd = 2)
  }
