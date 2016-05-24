DiscrFact <-
function (x, threshold = 1/10)
  {
    p <- x$int$dim[2]

    if (!any (class (x) == "tclust"))
      stop ("parameter x: expected object of type \x22tclust\x22")

    alpha <- x$par$alpha

    n <- nrow (x$par$x)
    ll <-matrix(ncol=x$k, nrow = n)
    disc <- disc2 <- ind <- array (NA, n)

    no.trim = floor(n*(1-alpha))

    for (k in 1:x$k)
      ll[,k] <- (x$size[k] / no.trim)* .dmnorm (x$par$x,x$centers[,k],
        as.matrix (x$cov[,,k]))

    llo <- apply (-ll, 1, order)  # the row-wise order of matrix ll 
    if (!is.matrix (llo))
      llo <- matrix (llo, ncol = nrow (ll)) ## -> transposed..


    disc <- ll[cbind (1:n, llo[1, ])]
    ind <- llo[1,]

    if (nrow (llo) >= 2)
    {
      disc2 <- ll[cbind (1:n, llo[2, ])]
      ind2<- llo[2,]
    }
    else
    {
      ind2 <- rep (0, n)
      disc2 <- disc * threshold
     }

    mropt <- sort (disc)[n - floor (n * (1 - alpha)) + 1] 
    idx.out <- disc < mropt

    ind2[idx.out] <- ind[idx.out]
    disc2[idx.out] <- disc[idx.out]
    disc[idx.out] <- mropt
    ind[idx.out] <- 0

    assignfact <- log (disc2 / disc)
    idx.fin <- is.finite (assignfact)
    ylimmin <- min (assignfact[idx.fin]) * 1.5
    assignfact[!idx.fin] <- ylimmin * 2

    mean.DiscrFact <- array (,x$k)
    for (i in 0:x$k)
      mean.DiscrFact[i + 1] = mean (assignfact[ind == i])
      
    names (mean.DiscrFact) <- c ("O", 1:x$k)

    ret <- list (x = x, ylimmin = ylimmin, ind = ind, ind2 = ind2, assignfact = 
                assignfact, disc = disc, disc2 = disc2, threshold = 
                log (threshold), mean.DiscrFact = mean.DiscrFact)
    class (ret) <- "DiscrFact"
    ret
  }

