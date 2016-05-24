thresholdIdentification <- function(g, x, y, n.cand = 3, covtype = "matern5_2", KM = NULL){
  
  if(is.null(KM)){
    KM <- km(~1, design = data.frame(x), response = y, covtype = covtype)
  }
  
  if(is.data.frame(y)){
    if (ncol(y) != 1) {stop("y must be one-dimensional")}
    y <- y[,1]
  }

  looAdditive <- function (x, y, parameter, covtype , eps.R = 1e-08, 
                           cl, iso = FALSE, se.compute = FALSE) 
  {
    n.cl <- length(cl)
    if (identical(iso, FALSE)) {
      iso <- rep(FALSE, n.cl)
    }
    if (n.cl == 1 & class(parameter) == "km") { 
      message("full clique: used DiceKriging::leaveOneOut.km")
      result <- leaveOneOut.km(parameter, type="UK")
    }
    else {
      if (n.cl == 1 & class(parameter) != "km") 
        stop("for only one clique a kriging model of class km is required for 'parameter'")
      parameter <- paramList2Vect(parameter, cl, iso)
      alpha <- parameter[(1:n.cl)]
      theta <- parameter[-(1:n.cl)]
      n <- length(y)
      x <- as.matrix(x)
      DM <- matrix(1, ncol = 1, nrow = n)
      p <- 1
      R <- Rfunc(as.matrix(x), theta, alpha, covtype, n, n.cl, 
                               cl, iso) + diag(eps.R, ncol = n, nrow = n)
      y.cv <- numeric(n)
      r.newdata <- rfunc(x, x, theta, alpha, covtype, n.cl, cl, iso)
      for (i in 1:n){
        newdata <- matrix(x[i,], ncol=ncol(x))
        x.but.i <- x[-i,]
        R.but.i <- R[-i,-i]
        DM.but.i <- DM[-i,]
        y.but.i <- y[-i]
        Rinvs <- solve(R.but.i)
        beta <- solve(t(DM.but.i) %*% Rinvs %*% DM.but.i) %*% t(DM.but.i) %*% Rinvs %*% y.but.i
        factor2 <- Rinvs %*% (y.but.i - DM.but.i %*% beta)
        r.i <- r.newdata[-i,i]
        y.cv[i] <- as.numeric(beta) + r.i %*% factor2
      }
      result <- data.frame(mean = y.cv)
    }
    return(result)
  }
  
  chooseDelta <- function(g, number = 3){
    tii.sort <- sort(g$tii.scaled[,1])
    jumps <- diff(tii.sort)
    vote <- which(jumps %in%  sort(jumps, decreasing=TRUE)[1:number])
    l <- tii.sort[vote]
    u <- tii.sort[vote+1]
    m <- as.numeric(1/2*(l+u)) # no names
    digit <- rep(0,length(l))
    while (any(a <- (l > round(m, digit) | u < round(m, digit)))){
      digit[which(a)] <- digit[which(a)] + 1
    }
    return(round(m, digit))
  }
  
  delta <- c(chooseDelta(g, number = n.cand), 1)
  
  graphs <- list(NULL)
  models <- list(NULL)
  y.cv <- list(NULL)
  RMSE <- numeric(length(delta)+1)

  op <- par("mfrow")
  nc <- round(sqrt(length(delta)+1))
  nl <- ceiling((length(delta)+1)/nc)
  mfrow <- c(nc, nl)
  par(mfrow=mfrow)
  
  models[[1]] <- KM
  y.cv[[1]] <- leaveOneOut.km(KM, type="UK")$mean
  RMSE[1] <- sqrt(sum((y - y.cv[[1]])^2))
  plot(y, y.cv[[1]], ylab= "crossvalidated y", main="threshold = 0", asp=TRUE)
  abline(0,1,col="red")
  
  for (i in 1:length(delta)){
    graphs[[i]] <- threshold(graphlist=g, delta=delta[i])
    models[[i+1]] <- kmAdditive(x, y, cl = graphs[[i]]$cliques, n.initial.tries=20, covtype = covtype, eps.R=1e-4,max.it=100)
    y.cv[[i+1]] <- looAdditive(x,y,models[[i+1]],cl=graphs[[i]]$cliques, covtype = covtype)$mean
    RMSE[i+1] <- sqrt(sum((y - y.cv[[i+1]])^2))
    plot(y, y.cv[[i+1]], ylab= "crossvalidated y", main=paste0("threshold = ", delta[i]), asp=TRUE)
    abline(0,1,col="red")
  }
  
  par(mfrow=op)
  print(cbind("threshold"=c(0,delta), "RMSE"=RMSE))
  return(list(threshold = c(0,delta), models = models, y.cv = y.cv, RMSE = RMSE))
} 
