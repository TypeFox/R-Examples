plotf <- function(dataf, main = "Functional data", xlab = "args", ylab = "vals", colors = c("red", "blue", "green", "black", "orange", "pink")) {
  
  if(main == "Functional data" && !is.null(dataf$name))
    main = dataf$name
  if(xlab == "args" && !is.null(dataf$args))
    xlab = dataf$args
  if(ylab == "vals" && !is.null(dataf$vals))
    ylab = dataf$vals
  
  ylims = matrix(unlist(lapply(dataf$dataf, function(e) (range(e$vals)))), ncol = 2, byrow = TRUE)
  plot(0, type="n", xlim=range(dataf$dataf[[1]]$args), ylim=c(min(ylims[,1]), max(ylims[,2])), 
       xlab=xlab, ylab=ylab, 
       main = main)
  grid()
  
  labs = sort(unlist(unique(dataf$labels)))
  
  for (i in 1:length(dataf$dataf)){
    ind = match(dataf$labels[[i]],labs)
    lineColor <- colors[ind]
    
    lines(dataf$dataf[[i]]$args, dataf$dataf[[i]]$vals, col=lineColor)
  }
}


dataf.sim.1.CFF07 <- function(numTrain = 100, numTest = 50, numDiscrets = 51, plot = FALSE){
  # Processes:
  # X(t) = m_0(t) + e(t), m_0(t) = 30*(1-t)*t^1.2
  # Y(t) = m_1(t) + e(t), m_1(t) = 30*(1-t)^1.2*t
  # e(t): Gaussian with mean = 0, cov(X(s), X(t)) = 0.2*exp(-abs(s - t)/0.3)
  
  t <- 0:(numDiscrets - 1)/(numDiscrets - 1)
  mean0 <- 30*(1-t)*t^1.2
  mean1 <- 30*(1-t)^1.2*t
  cov <- matrix(nrow=numDiscrets, ncol=numDiscrets)
  for (i in 1:numDiscrets){
    for (j in 1:numDiscrets){
      cov[i,j] <- 0.2*exp(-abs(t[i] - t[j])/0.3)
    }
  }
  
  X <- mvrnorm(n=numTrain+numTest, mu=mean0, Sigma=cov)
  Y <- mvrnorm(n=numTrain+numTest, mu=mean1, Sigma=cov)
  
  datafX <- list()
  datafY <- list()
  labelsX <- as.list(rep(0,numTrain+numTest))
  labelsY <- as.list(rep(1,numTrain+numTest))
  for (i in 1:(numTrain + numTest)){
    datafX[[i]] <- list(args = t, vals = X[i,])
    datafY[[i]] <- list(args = t, vals = Y[i,])
  }
  
  learn <- list(dataf = c(head(datafX, numTrain), head(datafY, numTrain)), 
                labels = c(head(labelsX, numTrain), head(labelsY, numTrain)))
  test <- list(dataf = c(tail(datafX, numTest), tail(datafY, numTest)), 
               labels = c(tail(labelsX, numTest), tail(labelsY, numTest)))
  if (plot){
    plot(0, type="n", xlim=c(0,1), ylim=c(0, 9), 
         main=paste("Model 1 from CuevasFF07: ", 
                    "0 red (", sum(unlist(learn$labels) == 0), "), ", 
                    "1 blue (", sum(unlist(learn$labels) == 1), "), ", sep=""))
    grid()
    for (i in 1:length(learn$dataf)){
      if (learn$labels[[i]] == 0){
        lineColor <- "red"
        lineType <- 1
      }
      if (learn$labels[[i]] == 1){
        lineColor <- "blue"
        lineType <- 2
      }
      lines(learn$dataf[[i]]$args, learn$dataf[[i]]$vals, col=lineColor, lty=lineType)
    }
  }
  
  return (list(learn = learn, test = test))
}

dataf.sim.2.CFF07 <- function(numTrain = 100, numTest = 50, numDiscrets = 51, plot = FALSE){
  # Processes
  # X(t) = m_0(t) + e(t), m_0(t) = 30*(1-t)*t^2 + 0.5*abs(sin(20*pi*t))
  # Y(t) = smooth.spline with 8 knots
  # e(t): Gaussian with mean = 0, cov(X(s), X(t)) = 0.2*exp(-abs(s - t)/0.3)
  
  t <- 0:(numDiscrets - 1)/(numDiscrets - 1)
  mean0 <- 30*(1 - t)*t^2 + 0.5*abs(sin(20*pi*t))
  cov <- matrix(nrow=numDiscrets, ncol=numDiscrets)
  for (i in 1:numDiscrets){
    for (j in 1:numDiscrets){
      cov[i,j] <- 0.2*exp(-abs(t[i] - t[j])/0.3)
    }
  }
  
  X <- mvrnorm(n=numTrain+numTest, mu=mean0, Sigma=cov)
  Y <- NULL
  for (i in 1:nrow(X)){
    Y <- rbind(Y, smooth.spline(t, X[i,], nknots = 8)$y)
  }
  
  datafX <- list()
  datafY <- list()
  labelsX <- as.list(rep(0,numTrain+numTest))
  labelsY <- as.list(rep(1,numTrain+numTest))
  for (i in 1:(numTrain + numTest)){
    datafX[[i]] <- list(args = t, vals = X[i,])
    datafY[[i]] <- list(args = t, vals = Y[i,])
  }
  
  learn <- list(dataf = c(head(datafX, numTrain), head(datafY, numTrain)), 
                labels = c(head(labelsX, numTrain), head(labelsY, numTrain)))
  test <- list(dataf = c(tail(datafX, numTest), tail(datafY, numTest)), 
               labels = c(tail(labelsX, numTest), tail(labelsY, numTest)))
  
  if (plot){
    plot(0, type="n", xlim=c(0,1), ylim=c(0, 7), 
         main=paste("Model 2 from CuevasFF07: ", 
                    "0 red (", sum(unlist(learn$labels) == 0), "), ", 
                    "1 blue (", sum(unlist(learn$labels) == 1), "), ", sep=""))
    grid()
    for (i in 1:length(learn$dataf)){
      if (learn$labels[[i]] == 0){
        lineColor <- "red"
        lineType <- 1
      }
      if (learn$labels[[i]] == 1){
        lineColor <- "blue"
        lineType <- 2
      }
      lines(learn$dataf[[i]]$args, learn$dataf[[i]]$vals, col=lineColor, lty=lineType)
    }
  }
  
  return (list(learn = learn, test = test))
}