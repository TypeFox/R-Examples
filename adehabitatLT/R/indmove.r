indmove <- function(ltr, nrep = 200, conflim = seq(0.95, 0.5, length=5),
                    sep = ltr[[1]]$dt[1], units=c("seconds", "minutes",
                                            "hours", "days"),
                    plotit = TRUE)
  {
    if (!inherits(ltr, "ltraj"))
      stop("ltr should be of class 'ltraj'")
    units <- match.arg(units)

    if ((!is.regular(ltr))&attr(ltr,"typeII"))
        stop("ltr should be regular or of type I")

    ## warning donnees manquantes
    na <- unlist(lapply(ltr, function(i) length(i$x[is.na(i$x)])))
    if (any(na>0)) {
      messa <- "Missing values have been interpolated in the data:\n"
      messa2 <- paste(unlist(lapply(1:length(ltr), function(i) paste("Burst ", attr(i, "burst"), ": ", na[i], "\n"))), collapse="")
      warning(paste(messa, messa2, "\n", sep=""))
    }

    ## Interpolation donnees manquantes
    ltr <- do.call("c.ltraj", lapply(ltr, function(x) {
      x1 <- approx(unclass(x$date)[!is.na(x$x)], x$x[!is.na(x$x)], xout=unclass(x$date))$y
      y1 <- approx(unclass(x$date)[!is.na(x$x)], x$y[!is.na(x$y)], xout=unclass(x$date))$y
      return(as.ltraj(data.frame(x=x1, y=y1), x$date, id=attr(x,"id"), burst=attr(x,"burst")))
    }))

    foo <- function(xy, xy2, dt, nrep = 500, units, main, plotit)
      {
        dt <- dt[-nrow(xy)]  # n-1
        if (units == "hours") dt <- dt/3600
        if (units == "days") dy <- dt/(3600*24)
        if (units == "minutes") dy <- dt/60

        xy <- xy[-nrow(xy),] # n-1
        toto <- .C("permutR2n", as.double(t(xy)), as.integer(nrow(xy)),
                   as.integer(nrep), double(nrep*nrow(xy)),
                   as.double(dt), double(nrep*nrow(xy)), PACKAGE = "adehabitatLT")
        dtr <- toto[[6]]
        toto <- toto[[4]]
        res <- matrix(toto, ncol=nrep, byrow = TRUE)
        dtr <- matrix(dtr, ncol=nrep, byrow = TRUE)
        di2 <- sapply(1:(nrow(xy2)),
                      function(i) return(((xy2[i,1] - xy2[1,1])^2) +
                                         ((xy2[1,2] - xy2[i,2])^2)))
        res <- cbind(di2[-1], res)
        dtr <- cbind(cumsum(dt), dtr)
        dtr <- rbind(rep(0,nrep+1), dtr)
        res <- rbind(rep(0,nrep+1), res)
        dtu <- unique(c(dtr))
        dtu <- dtu[order(dtu)]

        res2 <- .C("prepquart", as.double(dtu), as.integer(length(dtu)),
                   as.double(t(dtr)), as.double(t(res)),
                   as.integer(nrow(res)), as.integer(ncol(res)),
                   double(length(dtu)*ncol(res)), PACKAGE = "adehabitatLT")[[7]]
        res2 <- matrix(res2, nrow=length(dtu), byrow = TRUE)
        res2[res2<0] <- NA
        u <- res2[,-1]

        qua <- lapply(conflim, function(x) {
          lim3 <- c((1-x)/2, 1-(1-x)/2)
          return(t(apply(res2[,-1], 1, function(y)
                         quantile(y, probs=lim3, na.rm=TRUE))))
               }
          )

        ma <- max(unlist(lapply(qua,
                                function(x) max(c(c(x),
                                                  res2[!is.na(res2[,1]),1])))))

        if (plotit) {
          opar <- par(mar=c(5,5,4,2))
          on.exit(par(opar))
          plot(dtr[,1], res[,1], ty="n", ylim=c(0, max(ma)), xlab="time",
               ylab=expression(R(t)^2), main=main)
          lapply(1:length(qua), function(i)
                 polygon(c(dtu, dtu[length(dtu):1]),
                         c(qua[[i]][,1], qua[[i]][length(dtu):1,2]),
                         col=gray((length(qua)-i+1)/(length(qua)+length(qua)/20)), border=NA))
          lines(dtr[,1],res[,1], lwd=4, col="white")
          lines(dtr[,1],res[,1], lwd=2)
        }
        colnames(res) <- c("Obs",paste("Rep", 1:nrep, sep="."))
        colnames(dtr) <- c("Obs",paste("Rep", 1:nrep, sep="."))
        res <- list(time=dtr, R2n=res)
        return(res)
      }
    if (plotit) {
      opar <- par(mfrow=n2mfrow(length(ltr)))
      on.exit(par(opar))
    }
    oo <- lapply(1:length(ltr), function(x) {
      i <- ltr[[x]]
      uu <- foo(i[,c("dx","dy")],i[,c("x","y")],i$dt,
                nrep=nrep, units = units, main=attr(ltr[[x]],"burst"),
                plotit)
      return(uu)
    })
    names(oo) <- unlist(lapply(ltr, function(x) attr(x, "burst")))
    invisible(oo)
  }

