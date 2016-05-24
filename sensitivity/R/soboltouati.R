
# Khalid Boumhaout and Taieb Touati (2016)

soboltouati <- function(model = NULL, X1, X2, conf = 0.95, ...) {
  if ((ncol(X1) != ncol(X2)) | (nrow(X1) != nrow(X2)))
    stop("The samples X1 and X2 must have the same dimensions")
  p <- ncol(X1)
  
  X <- rbind(X1,X2)
  for (i in 1:p) {
    Xb <- X1
    Xb[,i] <- X2[,i]
    X <- rbind(X, Xb) 
  }
  
  x <- list(model = model, X1 = X1, X2 = X2, conf = conf, X = X,
            call = match.call())
  class(x) <- "soboltouati"
  
  if (!is.null(x$model)) {
    response(x, ...)
    tell(x)
  }
  
  return(x)
}


Pmv = function(x,Y,Z){
  cov11 = cov(Y,Z)
  cov21 = cov(Y^2,Z)
  cov12 = cov(Y,Z^2)
  cov22 = cov(Y^2,Z^2)
  cov31 = cov(Y^3,Z)
  cov13 = cov(Y,Z^3)
  
  EY = mean(Y)
  EY2 = mean(Y^2)
  EZ = mean(Z)
  EZ2 = mean(Z^2)
  
  sdY = sd(Y)
  sdZ = sd(Z)
  corYZ = cor(Y,Z)
  
  K11 = var((Y-EY)*(Z-EZ))
  K12 = (cov31 - 3*EY*cov21 - cov11*(EY2-4*EY^2))/(2*sdY)
  K13 = (cov13 - 3*EZ*cov12 - cov11*(EZ2-4*EZ^2))/(2*sdZ)
  K22 = (mean((Y-EY)^4) - sdY^4)/(4*sdY^2)
  K23 = (cov((Y-EY)^2,(Z-EZ)^2) + 2*EZ*cov((Y-EY)^2,(Z-EZ)) - 2*EZ*cov((Y-EY)^2, Y-EY))/(4*sdY*sdZ)
  K33 = (mean((Z-EZ)^4) - sdZ^4)/(4*sdZ^2)
  
  return(abs((K22/(sdY^2) + K33/(sdZ^2) + 2*K23/(sdY*sdZ))*x^2
             - 2*x*(K12/sdY + K13/sdZ)/(sdY*sdZ)
             + K11/((sdY*sdZ)^2)))
}


estim.soboltouati <- function(data, i = 1 : nrow(data), conf=0) {
  d <- as.matrix(data[i, ]) # as.matrix for colSums
  n <- nrow(d)
  p <- ncol(d)-2
  
  V <- var(d[, 1])
  ecor <- rep(0,p) ; ecorcompl <- rep(0,p)
  if(conf != 0) {
    VV <- matrix(V,nrow=1,ncol=3,dimnames=list(1,c("estim","CIinf","CIsup")))
    estcor <- matrix(0,nrow=p,ncol=3,dimnames=list(2:(p+1),c("estim","CIinf","CIsup")))
    estcorcompl <- matrix(0,nrow=p,ncol=3,dimnames=list((p+2):(2*p+1),c("estim","CIinf","CIsup")))
  }
  for(ii in 1:p) {
    ecor[ii] <- cor(d[,2],d[,ii+2],use="pairwise.complete.obs")
    ecorcompl[ii] <- cor(d[,1],d[,ii+2],use="pairwise.complete.obs")
    
    if(conf != 0) { ### Here goes PMV stuff
      estcor[ii,1] <- ecor[ii]
      tau2 = Pmv(estcor[ii,1],d[,2],d[,ii+2])
      estcor[ii,2] <- estcor[ii,1] - sqrt(tau2)*qnorm((1+conf)/2)/sqrt(n)
      estcor[ii,3] <- estcor[ii,1] + sqrt(tau2)*qnorm((1+conf)/2)/sqrt(n)
      
      estcorcompl[ii,1] <- ecorcompl[ii]
      tau2 = Pmv(estcorcompl[ii,1],d[,1],d[,ii+2])
      estcorcompl[ii,2] <- estcorcompl[ii,1] + sqrt(tau2)*qnorm((1+conf)/2)/sqrt(n)
      estcorcompl[ii,3] <- estcorcompl[ii,1] - sqrt(tau2)*qnorm((1+conf)/2)/sqrt(n)
    }
  }
  if(conf != 0) 
  { return(rbind(VV, estcor, estcorcompl))}
  else 
  { return(c(V, ecor, ecorcompl))}
}


tell.soboltouati <- function(x, y = NULL, return.var = NULL, ...) {
  id <- deparse(substitute(x))
  
  if (! is.null(y)) {
    x$y <- y
  } else if (is.null(x$y)) {
    stop("y not found")
  }
  
  p <- ncol(x$X1)
  n <- nrow(x$X1)
  
  data <- matrix(x$y, nrow = n)
  
  # estimation of the partial variances (V, D1 and Dt)
  
  V <- data.frame(original = estim.soboltouati(data, 1:n, x$conf))
  if(x$conf){colnames(V) <- c("original","min. c.i.","max. c.i.")}
  else{colnames(V) <- c("original")}
  
  # estimation of the Sobol' indices (S1 and St)
  
  S <- V[2:(p + 1),, drop = FALSE]
  T <- 1 - V[(p + 2):(2 * p + 1),, drop = FALSE]
    
  rownames(S) <- colnames(x$X1)
  rownames(T) <- colnames(x$X1)
  
  # return
  x$V <- V
  x$S <- S
  x$T <- T
  
  assign(id, x, parent.frame())
}


print.soboltouati <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (!is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    cat("\nFirst order indices:\n")
    if(!is.null(x$S$original.CIinf)) {
      df=data.frame(pointEstimate=x$S[,1],  minCI=x$S[,2], maxCI=x$S[,3])
      colnames(df)=c("estimate","min. c.i.","max. c.i.")
      print(df)
    } else {
      print(x$S)
    }
    cat("\nTotal indices:\n")
    if(!is.null(x$T$original.CIinf)) {
      df=data.frame(pointEstimate=x$T[,1], minCI=x$T[,2], maxCI=x$T[,3])
      colnames(df)=c("estimate","min. c.i.","max. c.i.")
      print(df)
    } else {
      print(x$T)
    }
  }
}


plot.soboltouati <- function(x, ylim = c(0, 1), ...) {
  if (!is.null(x$y)) {
    p <- ncol(x$X1)
    pch = c(21, 24)
    nodeplot(x$S, xlim = c(1, p + 1), ylim = ylim, pch = pch[1])
    nodeplot(x$T, xlim = c(1, p + 1), ylim = ylim, labels = FALSE,
             pch = pch[2], at = (1:p)+.3, add = TRUE)
    legend(x = "topright", legend = c("main effect", "total effect"), pch = pch)
  }
}