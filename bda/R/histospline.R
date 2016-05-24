## last updated 2015/04/16:  

.histospline <- function(x, from, to, gridsize=512L)
    {
        stopifnot(gridsize > 10)
        nbin <- length(x$mids)

        F <- x$counts; X <- x$mids; n <- sum(F)
        A <- x$breaks[1:nbin]; B <- x$breaks[-1]
        x0 <- runif(n,rep(A,F),rep(B,F))
        h <-  (1/(4*pi))^(1/10)*(243/(35*n))^(1/5)*sqrt(var(x0))*15^(1/5)
        
        a <- ifelse(missing(from), x$breaks[1] - 3*h, from)
        b <- ifelse(missing(to), x$breaks[nbin+1] + 3*h, to)
        
        out <- spline(x$mids, x$density, n=gridsize,
                      xmin = a, xmax = b)
        f0 <- out$y
        f0[f0<0] <- 0
        gpoints <- out$x
        
        structure(list(y = f0, x = gpoints,
                       conf.level = NULL,
                       type="histospline",
                       pars = NULL,
                       ucb=NULL, lcb=NULL,
                       call = match.call()
                       ),
                  class = 'histosmooth')

    }
