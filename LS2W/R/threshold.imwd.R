`threshold.imwd` <-
function (x, levels = 3:(x$nlevels - 1), type = "hard", policy = "universal",
    by.level = FALSE, value = 0, dev = var, verbose = FALSE,
    return.threshold = FALSE, compression = TRUE, Q = 0.05, ...)

{
#
#
#Check class of imwd
#
if(verbose == TRUE) cat("Argument checking\n")
   ctmp <- class(x)
if(is.null(ctmp))
   stop("data supplied has no class")
else if(ctmp != "imwd")
   stop("data supplied is not of class imwd")
if(policy != "universal" && policy != "manual" && policy != "probability" && policy != "fdr" && policy != "LSuniversal")
   stop("Only policys are universal, manual, fdr and probability at present")
if(type != "hard" && type != "soft")
   stop("Only hard or soft thresholding at  present")
r <- range(levels)
if(r[1.] < 0.)
   stop("levels out of range, level too small")
if(r[2.] > x$nlevels - 1.)
   stop("levels out of range, level too big")
if(r[1.] > x$nlevels - 1.) {
   warning("no thresholding done")
   return(x)
}
if(r[2.] < 0.) {
   warning("no thresholding done")
   return(x)
}
nthresh <- length(levels)
d <- NULL
#       Decide which policy to adopt
#               The next if-else construction should define a vector called
#               "thresh" that contains the threshold value for each level
#               in "levels". This may be the same threshold value
#               a global threshold.
#
n <- 2.^(2. * x$nlevels)
if(policy == "universal") {
   if(verbose == TRUE)
   cat("Universal policy...")
   if(by.level == FALSE) {
      if(verbose == TRUE)
         cat("All levels at once\n")
      for(i in 1.:nthresh) {
         d <- c(d, x[[lt.to.name(levels[i], "CD")]],
         x[[lt.to.name(levels[i], "DC")]],
         x[[lt.to.name(levels[i], "DD")]])
      }
      noise.level <- sqrt(dev(d))
      thresh <- sqrt(2. * log(n)) * noise.level
      if(verbose == TRUE)
         cat("Global threshold is: ", thresh, "\n")
         thresh <- rep(thresh, length = nthresh)
      }
      else {
      if(verbose == TRUE)
         cat("Level by level\n")
      thresh <- rep(0., length = nthresh)
      for(i in 1.:nthresh) {
         d <- c(x[[lt.to.name(levels[i], "CD")]],x[[lt.to.name(levels[i], "DC")]],x[[lt.to.name(levels[i], "DD")]])
         noise.level <- sqrt(dev(d))
         thresh[i] <- sqrt(2. * log(n)) * noise.level
         if(verbose == TRUE)
            cat("Threshold for level: ", levels[i], " is ", thresh[i], "\n")
      }
  }
}
if(policy == "LSuniversal") {
   if(verbose == TRUE)
      cat("Universal policy...")
   if(by.level == FALSE) {
      if(verbose == TRUE)
         cat("All levels at once\n")
      for(i in 1.:nthresh) {
         d <- c(d, x[[lt.to.name(levels[i], "CD")]],x[[lt.to.name(levels[i], "DC")]],x[[lt.to.name(levels[i], "DD")]])
      }
   noise.level <- sqrt(dev(d))
   thresh <- log(n) * noise.level
   if(verbose == TRUE)
   cat("Global threshold is: ", thresh, "\n")
   thresh <- rep(thresh, length = nthresh)
}
else {
   if(verbose == TRUE)
      cat("Level by level\n")
   thresh <- rep(0., length = nthresh)
   for(i in 1.:nthresh) {
      d <- c(x[[lt.to.name(levels[i], "CD")]],x[[lt.to.name(levels[i], "DC")]],x[[lt.to.name(levels[i], "DD")]])
      noise.level <- sqrt(dev(d))
      thresh[i] <- log(n) * noise.level
      if(verbose == TRUE)
         cat("Threshold for level: ", levels[i], " is ", thresh[i], "\n")
      }
   }
}
else if(policy == "manual") {
   if(verbose == TRUE)
      cat("Manual policy...\n")
   thresh <- rep(value, length = nthresh)
   if(length(value) != 1. && length(value) != nthresh)
      warning("your threshold is not the same length as number of levels")
}
else if(policy == "fdr") {
#
#
#               Threshold chosen by FDR-procedure
#
   if(verbose == TRUE) 
      cat("FDR policy...")
   if(by.level == FALSE) {
      if(verbose == TRUE)
         cat("All levels at once\n")
      for(i in 1.:nthresh) {
         d <- c(d, x[[lt.to.name(levels[i], "CD")]],x[[lt.to.name(levels[i], "DC")]],x[[lt.to.name(levels[i], "DD")]])
      }
      if(length(value) != 1.)
         stop("Length of value should be 1")
      noise.level <- sqrt(dev(c(x[[lt.to.name(levels[nthresh], "CD")]], x[[lt.to.name(levels[nthresh], 
"DC")]], x[[lt.to.name(levels[nthresh], "DD")]])))
      minit <- n
      dinit <- d
      thinit <- qnorm(1. - Q/2.) * noise.level
      if(log(n, 2.) > 15.)
      ninit <- 4.
      else {
         if(log(n, 2.) > 12.)
            ninit <- 3.
         else {
            if(log(n, 2.) > 10.)
               ninit <- 2.
            else ninit <- 1.
          }
      }
   for(k in seq(1., ninit)) {
      dinit1 <- dinit[abs(dinit) >= thinit]
      minit <- length(dinit1)
      if(minit == 0.)
         thresh <- max(abs(d)) * 1.0001
      else {
      thinit <- qnorm(1. - (Q * minit)/(2. * n)) * noise.level
      minit1 <- length(dinit1[abs(dinit1) >=thinit])
      if(minit1 == minit || minit1 == 0.)
         break
         dinit <- dinit1
      }
   }
   if(noise.level > 0.) {
      m <- length(d)
      minit <- length(dinit)
      p <- (2. - 2. * pnorm(abs(dinit)/noise.level))
      index <- order(p)
      j <- seq(1., minit)
      m0 <- max(j[p[index] <= (Q * j)/m])
      if(m0 != "NA" && m0 < minit)
         thresh <- abs(dinit[index[m0]])
      else {
        if(m0 == "NA")
           thresh <- max(abs(dinit)) *1.0001
      else thresh <- 0.
     }
   }
else thresh <- 0.
thresh <- rep(thresh, length = nthresh)
if(verbose == TRUE)
   cat("Global threshold is: ", thresh[1.], "\n","sigma is: ", noise.level, "\n")
}
else {
   if(verbose == TRUE)
     cat("Level by level\n")
     thresh <- rep(0., length = nthresh)
     for(i in 1.:nthresh) {
        d <- c(x[[lt.to.name(levels[i], "CD")]],x[[lt.to.name(levels[i], "DC")]],x[[lt.to.name(levels[i], "DD")]])
        m <- length(d)
        noise.level <- sqrt(dev(d))
        thinit <- qnorm(1. - Q/2.) * noise.level
        dinit <- d[abs(d) >= thinit]
        minit <- length(dinit)
        if(minit == 0.)
           thresh[i] <- max(abs(d)) * 1.0001
        else {
           if(noise.level > 0.) {
              p <- (2. - 2. * pnorm(abs(dinit)/noise.level))
             index <- order(p)
             j <- seq(1., minit)
             m0 <- max(j[p[index] <= (Q *j)/m])
             if(m0 != "NA" && m0 < minit)
                thresh[i] <- abs(dinit[index[m0]])
             else {
               if(m0 == "NA")
                  thresh[i] <-max(abs(dinit)) *1.0001
               else thresh[i] <- 0.
             }
        }
       else thresh[i] <- 0.
   }
   if(verbose == TRUE)
      cat("Threshold for level: ", levels[i], "is", thresh[i], "\n")
    }
  }
}
else if(policy == "probability") {
   if(verbose == TRUE)
      cat("Probability policy...")
   if(by.level == FALSE) {
      if(verbose == TRUE)
      cat("All levels at once\n")
         for(i in 1.:nthresh) {
            d <- c(d, x[[lt.to.name(levels[i], "CD")]],x[[lt.to.name(levels[i], "DC")]],x[[lt.to.name(levels[i], "DD")]])
         }
      if(length(value) != 1.)
         stop("Length of value should be 1")
      thresh <- rep(quantile(abs(d), prob = value), length = nthresh)
      if(verbose == TRUE)
         cat("Global threshold is: ", thresh[1.], "\n")
   }
   else { 
      if(verbose == TRUE)
         cat("Level by level\n")
      thresh <- rep(0., length = nthresh)
      if(length(value) == 1.)
         value <- rep(value, nthresh)
      if(length(value) != nthresh)
         stop("Wrong number of probability values")
      for(i in 1.:nthresh) {
         d <- c(x[[lt.to.name(levels[i], "CD")]],x[[lt.to.name(levels[i], "DC")]],x[[lt.to.name(levels[i], "DD")]])
         thresh[i] <- quantile(abs(d), prob = value[i])
         if(verbose == TRUE)
             cat("Threshold for level: ", levels[i], " is ", thresh[i], "\n")
      }
   }
}
if(return.threshold == TRUE)
   return(thresh)
for(i in 1.:nthresh) {
   dCD <- x[[lt.to.name(levels[i], "CD")]]
   dDC <- x[[lt.to.name(levels[i], "DC")]]
   dDD <- x[[lt.to.name(levels[i], "DD")]]
   if(type == "hard") {
      dCD[abs(dCD) <= thresh[i]] <- 0.
      dDC[abs(dDC) <= thresh[i]] <- 0.
      dDD[abs(dDD) <= thresh[i]] <- 0.
   if(verbose == TRUE) {
      cat("Level: ", levels[i], " there are ", sum(dCD == 0.), ":", sum(dDC == 0.), ":",sum(dDD == 0.), " zeroes and: ")
      cat(sum(dCD != 0.), ":", sum(dDC != 0.), ":",sum(dDD != 0.), " nonzeroes\n")
   }
}
else if(type == "soft") {
   dCD <- sign(dCD) * (abs(dCD) - thresh[i]) * (abs(dCD) >thresh[i])
   dDC <- sign(dDC) * (abs(dDC) - thresh[i]) * (abs(dDC) >thresh[i])
   dDD <- sign(dDD) * (abs(dDD) - thresh[i]) * (abs(dDD) >thresh[i])
   if(verbose == TRUE) {
      cat("Level: ", levels[i], " there are ", sum(dCD == 0.), ":", sum(dDC == 0.), ":",sum(dDD == 0.), " zeroes and: ")
      cat(sum(dCD != 0.), ":", sum(dDC != 0.), ":",sum(dDD != 0.), " nonzeroes\n")
     }
   }
   x[[lt.to.name(levels[i], "CD")]] <- dCD
   x[[lt.to.name(levels[i], "DC")]] <- dDC
   x[[lt.to.name(levels[i], "DD")]] <- dDD
}
if(compression == TRUE)
      return(compress(x, verbose = verbose))
else return(x)
}

