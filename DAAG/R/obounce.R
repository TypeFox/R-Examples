"obounce" <-
     function (x, d)
{
     ord <- order(x)
     xsort <- x[ord]
     n <- length(x)
     xnew <- xsort
     if (n > 1) {
         i1 <- 1
         while (i1 < n) {
             x1 <- xsort[i1]
             i2 <- i1 + 1
             for (j in i2:n) {
                 nobounce <- TRUE
                 jn <- n - j + i2
                 dj <- xsort[jn] - x1
                 dsought <- (jn - i1) * d
                 if (dj < dsought) {
                     jot <- (dsought - dj)/2
                     for (k in i1:jn) xnew[k] <- x1 - jot + (k -
                                                             i1) * d
                     i1 <- jn + 1
                     nobounce <- FALSE
                     break
                 }
             }
             if (nobounce)
                 i1 <- i1 + 1
         }
         if (min(diff(xnew)) < d * 0.999) {
             n1 <- (1:(n - 1))[diff(xnew) < d]
             cat("Error in bounce().  Improperly separated points are:",
                 fill = TRUE)
             cat(paste(n1, ":", n1 + 1, sep = ""), fill = TRUE)
             cat(paste(xnew[n1], ":", xnew[n1 + 1], sep = ""),
                 fill = TRUE)
         }
     }
     x[ord] <- xnew
     x
}
