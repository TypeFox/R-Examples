#
cthreshold <- function(alpha,nbtest) return(1-(1-alpha)^(nbtest^-1))
#
minpermute <- function(alpha,nbtest,margin=1,ru=3) {
  return(round(floor(margin*(1-(1-alpha)^(nbtest^-1))^-1),-ru)+(10^ru)-1)
}
#
dtau <- function(x, nu, tol = .Machine$double.eps ^ 0.5) {
  res <- numeric(length(x))
  nu <- rep(nu, length.out = length(x))
  # Integrand function to be integrated over positive real numbers.
  f <- function(z, x, nu)
    dt(z, nu) * dt(x / z, nu) / z   # abs(z) always == z for real positive.
  # The integrand being symmetric, only the positive is integrated and the result multiplied by two.
  for(i in 1L:length(x))
    res[i] <- 2 * integrate(f, lower = 0, upper = Inf, x = x[i], nu = nu[i], rel.tol = tol)$value
  res
}
#
ptau <- function(q, nu, lower.tail = TRUE, tol = .Machine$double.eps^0.5) {
  res <- rep(0.5, length(q))
  nu <- rep(nu, length.out = length(q))
  # Only the positive 
  signq <- sign(q)
  q <- abs(q)
  # Code avoid calculate PDF(0) because the function has a singularity at that point.
  if(lower.tail) {
    for(i in which(!!signq)) {
      if(signq[i]==1)
        res[i] <- 1 - integrate(dtau, lower = q[i], upper = Inf, nu = nu[i], tol = tol, rel.tol = tol)$value
      else
        res[i] <- integrate(dtau, lower = q[i], upper = Inf, nu = nu[i], tol = tol, rel.tol = tol)$value
    }
  } else {
    for(i in which(!!signq)) {
      if(signq[i]==1)
        res[i] <- integrate(dtau, lower = q[i], upper = Inf, nu = nu[i], tol = tol, rel.tol = tol)$value
      else
        res[i] <- 1 - integrate(dtau, lower = q[i], upper = Inf, nu = nu[i], tol = tol, rel.tol = tol)$value
    }
  }
  res
}
#
MCA <- function(y, x, emobj) {
  if(!inherits(emobj,"eigenmap")) stop("Parameter 'emobj' must be a 'eigenmap' object!")
  if(length(y) != length(x)) stop("Number of observations in y and x do not match!")
  if(nrow(emobj$U) != length(y)) stop("Number of observations in y does not match the number of lines in U.")
  #
  mssd <- c(my=NA,mx=NA,ssdy=NA,ssdx=NA)
  mssd[1:2] <- c(mean(y),mean(x))
  yxc <- cbind(yc=y-mssd[1], xc=x-mssd[2])
  mssd[3:4] <- c(t(yxc[,1]) %*% yxc[,1], t(yxc[,2]) %*% yxc[,2])  
  Upyxcb <- matrix(NA, dim(emobj$U)[2], 4)  
  colnames(Upyxcb) <- c("Upy", "Upx", "C", "B")
  rownames(Upyxcb) <- colnames(emobj$U)
  Upyxcb[,1:2] <- t(emobj$U) %*% yxc
  Upyxcb[,3] <- (Upyxcb[,1] * Upyxcb[,2]) / sqrt(mssd[3] * mssd[4])
  Upyxcb[,4] <- (Upyxcb[,1] / Upyxcb[,2])
  #
  ## Output block.
  return(structure(list(
    data = cbind(y=y,x=x),
    emobj = emobj,
    Upyxcb = Upyxcb,
    test = NULL),
    class = "mca"))
}
#
test.mca <- function (mcaobj, alpha = 0.05, max.step = Inf) {
  if (class(mcaobj) != "mca")
    stop("Parameter 'mcaobj' must be of class 'mca'.")
  if (!is.finite(max.step[1L]))
    max.step <- ncol(mcaobj$emobj$U)
  us <- matrix(NA, nrow(mcaobj$emobj$U), 0L)
  uspyx <- matrix(NA, 0L, 2L)
  yxc <- cbind(yc = mcaobj$data[,1L] - mean(mcaobj$data[,1L]), 
               xc = mcaobj$data[,2L] - mean(mcaobj$data[,2L]))
  ord <- order(abs(mcaobj$Upyxcb[,3L]), decreasing = TRUE)
  ttable <- matrix(NA, 0L, 4L)
  colnames(ttable) <- c("tau", "ddf", "Testwise p", "Familywise p")
  step <- 1L
  while (step != 0L) {
    us <- cbind(us, mcaobj$emobj$U[,ord[step]])
    uspyx <- rbind(uspyx, mcaobj$Upyxcb[ord[step],1L:2])
    ddfr <- nrow(mcaobj$data) - step - 1L
    ryx <- yxc - (us %*% uspyx)
    tau <- ddfr * (uspyx[step,1L] * uspyx[step,2L]) / sqrt(sum(ryx[,1L]^2) * sum(ryx[,2L]^2))
    ttable <- rbind(ttable, c(tau, ddfr, NA, NA))
    ttable[step,3L] <- 2 * ptau(abs(tau), ddfr, lower.tail = FALSE)
    ttable[step,4L] <- 1 - (1 - ttable[step,3L])^(ncol(mcaobj$emobj$U) - step + 1L)
    if (ttable[step,4L] > alpha || step >= max.step) {
      rownames(ttable) <- colnames(mcaobj$emobj$U)[ord[1L:step]]
      step <- 0L
    } else step <- step + 1L
  }
  signif <- ord[which(ttable[,4L] <= alpha)]
  return(structure(list(data = mcaobj$data,
                        emobj = mcaobj$emobj,
                        Upyxcb = mcaobj$Upyxcb,
                        test = list(permute = FALSE,
                                    significant = signif,
                                    test.table = ttable,
                                    details = NULL)),
                   class = "mca"))
}
#
permute.mca <- function (mcaobj, permute = NA, alpha = 0.05, max.step = Inf) {
  if (!is.finite(max.step[1L]))
    max.step <- ncol(mcaobj$emobj$U)
  if (is.na(permute[1L]))
    permute <- minpermute(alpha, ncol(mcaobj$emobj$U), 10L, 3L)
  us <- matrix(NA, nrow(mcaobj$emobj$U), 0L)
  uspyx <- matrix(NA, 0L, 2L)
  yxc <- cbind(yc = mcaobj$data[,1L] - mean(mcaobj$data[,1L]), 
               xc = mcaobj$data[,2L] - mean(mcaobj$data[,2L]))
  ord <- order(abs(mcaobj$Upyxcb[, 3]), decreasing = TRUE)
  ttable <- matrix(NA, 0L, 4L)
  colnames(ttable) <- c("tau", "ddf", "Testwise p", "Familywise p")
  details <- matrix(NA, 0L, 3L)
  colnames(details) <- c("tau* <= -|tau|", "-|tau| < tau* < |tau|", 
                         "tau* >= |tau|")
  step <- 1L
  while (step != 0L) {
    us <- cbind(us, mcaobj$emobj$U[,ord[step]])
    uspyx <- rbind(uspyx, mcaobj$Upyxcb[ord[step],1L:2])
    ddfr <- nrow(mcaobj$data) - step - 1L
    ryx <- yxc - (us %*% uspyx)
    tau0 <- (uspyx[step,1L] * uspyx[step,2L]) / sqrt(sum(ryx[,1L]^2) * sum(ryx[,2L]^2))
    ttable <- rbind(ttable, c(ddfr * tau0, ddfr, NA, NA))
    details <- rbind(details, c(0L, 0L, 1L))
    details[step,] <- .C("mcapermute",
                         as.double(tau0),
                         as.double(ryx[,1L]),
                         as.double(ryx[,2L]),
                         as.double(us[,step]),
                         as.integer(nrow(us)),
                         details = as.integer(details[step,]),
                         as.integer(permute))$details
    ttable[step,3L] <- (details[step,1L] + details[step,3L])/sum(details[step,])
    ttable[step,4L] <- 1 - (1 - ttable[step, 3])^(ncol(mcaobj$emobj$U) - step + 1L)
    if (ttable[step,4L] > alpha || step >= max.step) {
      rownames(ttable) <- rownames(details) <- colnames(mcaobj$emobj$U)[ord[1L:step]]
      step <- 0L
    } else step <- step + 1L
  }
  signif <- ord[which(ttable[, 4] <= alpha)]
  return(structure(list(data = mcaobj$data,
                        emobj = mcaobj$emobj,
                        Upyxcb = mcaobj$Upyxcb,
                        test = list(permute = permute,
                                    significant = signif,
                                    test.table = ttable,
                                    details = details)),
                   class = "mca"))
}
#
print.mca <- function(x, ...) {
  cat("\nMulti-scale Codependence Analysis\n---------------------------------\n\n")
  cat("Coefficients:\n")
  print(cbind(round(x$Upyxcb[,3:4],5), Lambda=round(x$emobj$lambda,5)))
  cat("\n")
  return(invisible(NULL))
}
#
summary.mca <- function(object, ...) {
  if(is.null(object$test)) {
    cat("\nNo testing informations available\n\n")
  } else {
    cat("\nTest table:\n")
    print(object$test$test.table)
    cat("\n")
  }
  return(invisible(NULL))
}
#
plot.mca <- function(x, ...) {
  cc <- rep(grey(0.5),nrow(x$Upyxcb))
  if(!is.null(x$test)) cc[x$test$significant] <- grey(0)
  barplot(x$Upyxcb[,3],names.arg=rownames(x$Upyxcb),ylab="C",
          ylim=c(-1,1)*max(abs(x$Upyxcb[,3])), las=2, space=0, col = cc)
  return(invisible(NULL))
}
#
fitted.mca <- function(object, which=NA, components=FALSE, ...) {
  if(!is.null(object$test)) which <- object$test$significant
  else if(is.na(which[1])) stop("No testing informations available: user must identify relevant coefficients.")
  fit <- matrix(0,nrow(object$data),1)
  if(components) {
    cpns <- matrix(NA, nrow(object$data), length(which))
    rownames(cpns) <- rownames(object$emobj$U)
  }
  if(length(which) != 0) {
    by <- object$Upyxcb[which,4] * object$Upyxcb[which,2]
    fit <- object$emobj$U[,which] %*% cbind(by) + mean(object$data[,1])
    if (components) {
      for (i in 1:length(which)) cpns[,i] <- object$emobj$U[,which[i]] * by[i]
      colnames(cpns) <- paste("Component", which)
    }
  }
  colnames(fit) <- "fitted" ; rownames(fit) <- rownames(object$emobj$U)
  if(components) {
    return(list(fitted=fit, components=cpns))
  } else return(fit)
}
#
residuals.mca <- function(object, which=NA, ...) {
  if(!is.null(object$test)) which <- object$test$significant
  else if(is.na(which[1])) stop("No testing informations available: user must identify relevant coefficients.")
  res <- object$data[,1]
  if(length(which) != 0) {
    by <- object$Upyxcb[which,4] * object$Upyxcb[which,2]
    res <- res - object$emobj$U[,which] %*% cbind(by) - mean(object$data[,1])
  }
  colnames(res) <- "residuals" ; rownames(res) <- rownames(object$emobj$U)
  return(res)
}
#
predict.mca <- function(object, which=NA, newdata=NA, components=FALSE, ...) {
  if(is.na(newdata[1])) return (fitted.mca(object, which=which))
  if(!is.null(object$test)) which <- object$test$significant
  else if(is.na(which[1])) stop("No testing informations available: user must identify relevant coefficients.")
  if(is.matrix(newdata)) {
    newdata <- newdata[,1] ; warning("Only the first row of the matrix provided as 'newdata' is used.")
  }
  if(length(newdata) != nrow(object$emobj$U)) {
    stop("Number of observations in 'newdata' does not match the number of lines in U.")
  }
  mnew <- mean(newdata) ; newc <- cbind(newc=newdata-mnew)
  Upnew <- t(object$emobj$U) %*% newc
  pred <- matrix(0,nrow(object$data),1)
  if(components) {
    cpns <- matrix(NA, nrow(object$data), length(which))
    rownames(cpns) <- rownames(object$emobj$U)
  }
  if(length(which) != 0) {
    by <- object$Upyxcb[which,4] * Upnew[which]
    pred <- object$emobj$U[,which] %*% cbind(by) + mnew * mean(object$data[,1]) / mean(object$data[,2])
    if(components) {
      for (i in 1:length(which)) cpns[,i] <- object$emobj$U[,which[i]] * by[i]
      colnames(cpns) <- paste("Component", which)
      }
  }
  colnames(pred) <- "predicted" ; rownames(pred) <- rownames(object$emobj$U)
  if(components) {
    return(list(predicted=pred, components=cpns))
  } else return(pred)
}
#
