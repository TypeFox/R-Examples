"EMPIRsim" <-
function(n=100, empgrid=NULL, kumaraswamy=FALSE, na.rm=TRUE, keept=FALSE,
                graphics=TRUE, ploton=TRUE, points=TRUE, snv=FALSE,
                infsnv.rm=TRUE, trapinfsnv=.Machine$double.eps, ...) {
  if(! graphics) {
     ploton <- FALSE
     points <- FALSE
  }
  empinv <- EMPIRgridderinv(empgrid=empgrid, kumaraswamy=kumaraswamy)
  rows <- as.numeric(attributes(empinv)$rownames)
  ix <- 1:length(rows)
  cols <- attributes(empinv)$colnames
  U <- vector(mode="numeric", length=n); V <- U
  Ts <- runif(n)
  for(i in 1:n) {
    u <- runif(1); t <- Ts[i]
    ix.needed1 <- max(ix[rows <= u])
    ix.needed2 <- min(ix[rows >= u])
    if(ix.needed1 == 1) ix.needed1 <- 2
    if(ix.needed1 == ix.needed2) {
      v.available <- empinv[ix.needed1,]
      v <- approx(cols, y=v.available, xout=t, rule=2)$y
    } else {
      v.available1 <- empinv[ix.needed1,]
      v1 <- approx(cols, y=v.available1, xout=t, rule=2)$y
      v.available2 <- empinv[ix.needed2,]
      v2 <- approx(cols, y=v.available2, xout=t, rule=2)$y
      w1 <- u - rows[ix.needed1]
      w2 <- rows[ix.needed2] - u
      tw <- 1/w1 + 1/w2
      v <- (v1/w1 + v2/w2)/tw
    }
    U[i] <- u; V[i] <- v
  }

  # Because z is a data.frame, it must be assigned within the ifelse()
  ifelse(keept, z <- data.frame(U=U, V=V, T=Ts), z <- data.frame(U=U, V=V))
  if(na.rm) {
     z <- z[complete.cases(z), ]
     m <- length(z[,1])
     if(m != n) {
        warning("user requested n=",n," simulations but only m=",m,
                " could be made without NA from EMPIRgridderinv")
        row.names(z) <- NULL # reset the rows to "1:m"
     }
  }
  if(snv) {
     if(infsnv.rm) {
        z <- z[z$U != 0, ]; z <- z[z$U != 0, ]
        z <- z[z$V != 1, ]; z <- z[z$V != 1, ]
        row.names(z) <- NULL
     } else if(trapinfsnv) {
        z$U[z$U == 0] <-   trapinfsnv; z$V[z$V == 0] <-   trapinfsnv
        z$U[z$U == 1] <- 1-trapinfsnv; z$V[z$V == 1] <- 1-trapinfsnv
     }
     z$U <- qnorm(z$U)
     z$V <- qnorm(z$V)
  }
  if(ploton) {
     if(snv) {
        plot(z$U, z$V, type="n",
             xlab="STANDARD NORMAL SCORE FOR U", ylab="STANDARD NORMAL SCORE FOR V")
     } else {
        plot(NA, NA, type="n", xlim=c(0,1), ylim=c(0,1),
             xlab="U, NONEXCEEDANCE PROBABILITY", ylab="V, NONEXCEEDANCE PROBABILITY")
     }
  }
  if(points & ! is.null(dev.list())) points(z$U, z$V,...)

  return(z)
}

