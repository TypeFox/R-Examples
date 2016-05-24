chi <- 
  # Much of the code in here was written by Alec Stephenson
  # and comes from his 'chiplot' function in his 'evd' package.
  # Minor differences between the evd implementation and this are: 
  # use of ties.method="first" here as oppsed to the evd method 
  # which used ties.method="average"; lower bound on chibar wrong 
  # in evd package.
function (data, nq = 100, qlim = NULL, alpha = 0.05, trunc = TRUE) {
    
    theCall <- match.call()

    eps <- .Machine$double.eps^0.5

    data <- na.omit(data)
    n <- nrow(data)

    # Get the EDFs
    t.method <- "first"

	if (is.R()){
	    data <- cbind(rank(data[, 1],ties.method = t.method)/(n + 1), 
    	              rank(data[, 2],ties.method = t.method)/(n + 1))
	}
	else {
		data <- cbind(rank(data[, 1])/(n + 1), 
    	              rank(data[, 2])/(n + 1))
	}

    rowmax <- apply(data, 1, max)
    rowmin <- apply(data, 1, min)

    qlim2 <- c(min(rowmax) + eps, max(rowmin) - eps)

    if (!is.null(qlim)) {
        if (qlim[1] < qlim2[1]){ stop("lower quantile limit is too low") }
        if (qlim[2] > qlim2[2]){ stop("upper quantile limit is too high") }
        if (qlim[1] > qlim[2]){ stop("lower quantile limit is less than upper quantile limit") }
    }
    else{
        qlim <- qlim2
    }

    u <- seq(qlim[1], qlim[2], length = nq)

    # HS. Replaced 2 for loops with calls to sapply

    cu <- sapply(1:nq, function(i, x, y){ mean(y < x[i]) }, y=rowmax, x=u )
    cbaru <- sapply(1:nq, function(i, x, y){ mean(y > x[i]) }, y=rowmin, x=u )

    # Get \chi and \bar\chi
    chiu <- 2 - log(cu)/log(u) # 3.2 of Coles, Heffernan, Tawn
    chibaru <- (2 * log(1 - u))/log(cbaru) - 1 # Page 348 of Coles, Heffernan, Tawn

    # Get confidence limits
    varchi <- ((1/log(u)^2 * 1)/cu^2 * cu * (1 - cu))/n
    varchi <- qnorm(1 - alpha/2) * sqrt(varchi)

    varchibar <- (((4 * log(1 - u)^2)/(log(cbaru)^4 * cbaru^2)) * cbaru * (1 - cbaru))/n
    varchibar <- qnorm(1 - alpha/2) * sqrt(varchibar)

    chiu <- cbind(chilow = chiu - varchi,           # Lower
                  chi = chiu,                       # Point est.
                  chiupp = chiu + varchi)           # Upper

    chibaru <- cbind(chiblow = chibaru - varchibar, # Lower
                     chib = chibaru,                # Point est.
                     chibupp = chibaru + varchibar) # Upper

    chiulb <- 2 - log(pmax(2 * u - 1, 0))/log(u)

    if (trunc) {
        chiu[chiu > 1] <- 1
        chiu <- apply(chiu, 2, function(x, z){ pmax(x, z) }, z = chiulb)
        chibaru[chibaru > 1] <- 1
        chibaru[chibaru < -1] <- -1
    }

    res <- list(chi=chiu, chibar = chibaru, quantile=u, call=theCall, 
                qlim=qlim, chiulb = chiulb)
    oldClass(res) <- "chi"
    res
} # Close chi <- function


print.chi <- function(x, ...){
    print(x$call,...)
    cat("Values of chi and chi-bar obtained and",
         length(x$quantile), "quantiles.\n")
    invisible()
}

summary.chi <- function(object, digits=3, ...){
    print(object$call)
    cat("Values of chi and chi-bar obtained and",
         length(object$quantile), "quantiles.\n")

    wh <- quantile(object$quantile, prob=c(.05, .5, .95))
    wh <- sapply(wh, function(i, u){
                        d <- abs(u - i)
                        min(u[d == min(d)])
                     }, u=object$quantile)

    chiQ <- object$chi[object$quantile %in% wh, 2]
    chibarQ <- object$chibar[object$quantile %in% wh, 2]

    res <- rbind(wh, chiQ, chibarQ)
    dimnames(res) <- list(c("quantile", "chi", "chi-bar"), rep("", 3))
    print(res, digits=digits, ...)
    invisible(res)
}


plot.chi <- function(x, show=c("Chi"=TRUE,"ChiBar"=TRUE), lty = 1, cilty = 2, col = 1, spcases = TRUE, cicol = 1,
                     xlim = c(0, 1), ylimChi = c(-1, 1), ylimChiBar = c(-1, 1),
                     mainChi = "Chi", mainChiBar = "Chi Bar",
                     xlab = "Quantile", 
                     ylabChi = expression(chi(u)),#"Chi",
                     ylabChiBar = expression(bar(chi)(u)), #"Chi Bar",
                     ask, ...){

  lty <- c(cilty, lty, cilty)
  col <- c(cicol, col, cicol)
  nb.fig <- prod(par("mfcol"))

	if (is.R() & missing(ask)){
	  ask <- nb.fig < sum(show) && dev.interactive()
	}	else {
		ask <- FALSE
	}
  if (ask) {
     op <- par(ask = TRUE)
     on.exit(par(op))
  }
  ChiBarAsympIndep <- prod(tail(x$chibar[,3]) < 1)
  if (show["ChiBar"]) {
    matplot(x$quantile, x$chibar, type = "l", lty = lty, col = col, 
            xlim = xlim, ylim = ylimChiBar, main = mainChiBar, xlab = xlab, 
            ylab = ylabChiBar, ...)
    if (spcases) {
      segments(x$qlim[1], 0, x$qlim[2], 0, lty = 5, col = "grey")
      segments(x$qlim[1], 1, x$qlim[2], 1, lty = 5, col = "grey")
      segments(x$qlim[1],-1, x$qlim[2],-1, lty = 5, col = "grey")
    }
  }
  if (show["Chi"]) {
    if (ChiBarAsympIndep) {
      col <- "grey"
      cols <- "grey"
    } else {
      cols <- "black"
    }
  
    matplot(x$quantile, x$chi, type = "l", lty = lty, col = col, xlim = xlim, 
            ylim = ylimChi, main = mainChi, xlab = xlab, ylab = ylabChi, 
            col.lab=cols,col.main=cols,col.sub=cols,axes=FALSE,
             ...)
     box(col=cols)
     axis(1,col=cols,col.axis = cols)
     axis(2,col=cols,col.axis = cols)
             
     if (spcases) {
       segments(x$qlim[1], 0, x$qlim[2], 0, lty = 5, col = "grey")
       segments(x$qlim[1], 1, x$qlim[2], 1, lty = 5, col = "grey")
       lines(x$quantile, x$chiulb, lty = 5, col = "grey")
     }
  }

  invisible()
}

