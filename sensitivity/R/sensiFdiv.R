# Sensitivity indices based on Csiszar f-divergence (Da Veiga 2014)
#
# Sebastien Da Veiga 2014

Compute_fdiv <- function(t,fdiv){
  switch(fdiv,
         TV={
           res <- abs(t-1)
         },
         KL={
           res <- -log(t)
         },
         Hellinger={
           res <- (sqrt(t)-1)^2
         },
         Chi2={
           res <- (1-t)^2/t
         },
         stop("invalid argument 'fdiv', not a valid f-divergence")
  )
  return(res)
}

Compute_ratio <- function(x,y){
  if (requireNamespace("ks", quietly = TRUE)){
    fx <- ks::kde(x,eval.points=x)$estimate
    fy <- ks::kde(y,eval.points=y)$estimate
    fxy <- ks::kde(cbind(x,y),eval.points=cbind(x,y))$estimate
  }
  return(fxy/(fx*fy))
}

sensiFdiv <- function(model = NULL, X, fdiv="TV", nboot = 0, conf = 0.95, ...) {
  
  
  if (is.data.frame(X)){
    X <- as.matrix(unname(X))
  }
  
  else if(!is.matrix(X)){
    stop("The sample X must be a matrix or a data frame")
  }
  
  x <- list(model = model, X = X, fdiv = fdiv, nboot = nboot,
            conf = conf, call = match.call()) 
  class(x) <- "sensiFdiv"
  
  #calcul of the response for explicit model
  if (! is.null(x$model)) {
    response(x, ...)
    x=tell(x, ...)
  }  
  return(x)
}


estim.sensiFdiv <- function(data, i=1:nrow(data), fdiv) {
  
  ptot <- ncol(data)
  p <- ptot - 1
  X <- data[i,1:p]
  Y <- data[i,ptot]
  nf <- length(fdiv)
  S = matrix(0,nrow=p,ncol=nf)
  
  # Csiszar f-divergence indices
  for (i in 1:p){
    for (k in 1:nf){
      Xtemp <- as.matrix(X[,i])
      Ytemp <- as.matrix(Y)
      ratio <- Compute_ratio(Xtemp,Ytemp)
      fdiv_sample <- Compute_fdiv(1/ratio,fdiv[k])
      S[i,k] <- mean(fdiv_sample,na.rm = TRUE)
    }
  }
  dim(S) <- c(p*nf,1)
  return(S)
}


tell.sensiFdiv <- function(x, y = NULL, ...) {
  
  id <- deparse(substitute(x))
  if (! is.null(y)) {
    x$y <- y
  } 
  else if (is.null(x$y)) {
    stop("y not found")
  }
  
  n <- nrow(x$X)
  p <- ncol(x$X)
  nf <- length(x$fdiv)
  data <-cbind(x$X,x$y)
  
  
  if (x$nboot == 0) {
    res <- estim.sensiFdiv(data, 1:n, x$fdiv)
    dim(res) <- c(p,nf)
    if (nf == 1) {
      x$S <- data.frame(res)
      colnames(x$S) <- "original"
    }
    else{
      x$S <- list()
      for (i in 1:nf){
        x$S[[i]] <-  data.frame(res[,i])
        colnames(x$S[[i]]) <- "original"
      }
      names(x$S) <- x$fdiv
    }
  } 
  else {
    S.boot <- boot(data, estim.sensiFdiv, fdiv=x$fdiv, R = x$nboot)
    Stemp <- bootstats(S.boot, x$conf, "basic")
    if (nf == 1){
      x$S <- Stemp
    }
    else{
      x$S <- list()
      for (i in 1:nf){
        x$S[[i]] <-  Stemp[((i-1)*p+1):(i*p),]
      }
      names(x$S) <- x$fdiv
    }
  }
  
  
  rownames <- c()
  for (i in 1:p) {
    rownames <- c(rownames,paste("S",i,sep=""))
  }
  if (nf ==1) {
    rownames(x$S) <- rownames
  }
  else{
    for (i in 1:nf){
      rownames(x$S[[i]]) <- rownames
    }
  }
  
  
  assign(id, x, parent.frame())
  return(x)
}


print.sensiFdiv <- function(x, ...) {
  
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    if (! is.null(x$S)) {
      nf <- length(x$fdiv)
      if (nf==1){
        cat(paste("\n\n\nCsiszar f-divergence indices with ",x$fdiv," \n",sep=""))
        print(x$S)
      }
      else {
        for (i in 1:nf){
          cat(paste("\n\n\nCsiszar f-divergence indices with ",x$fdiv[i]," \n",sep=""))
          print(x$S[[i]])
        }
      }
    }
    else{
      cat("(empty)\n")
    }
  }
}


plot.sensiFdiv <- function(x, ylim = c(0, 1), ...) {
  
  if (! is.null(x$y)) {
    nf <- length(x$fdiv)
    if (nf==1){
      nodeplot(x$S, ylim = ylim)
      legend(x = "topright", legend = paste("Csiszar f-divergence indices with", x$fdiv,sep=" "))
    }
    else {
      for (i in 1:nf){
        nodeplot(x$S[[i]], ylim = ylim)
        legend(x = "topright", legend = paste("Csiszar f-divergence indices with", x$fdiv[i],sep=" "))
      }
    }
  }
}
