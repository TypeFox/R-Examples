##This file contains all varieties of the combine functions
##Two numeric vectors by Bayes
combine.bs <- function(x,y)
    {
        lx <- log(x)
        ly <- log(y)
        xy <- exp(lx + ly)
        z <- xy / sum(xy)
        return(z)
    }

##Two numeric vectors by Dempster-Shafer
combine.ds <- function(x, y)
    {
        N <- length(x)
        z <- tmp1 <- tmp2 <- numeric(N)
        for(i in 1:N)
            {
                tmp1 <- x[i] * y
                tmp2 <- x * y[i]
                z[i] <- z[i] + tmp1[i] + tmp1[N] + tmp2[N]
            }
        z[N] <- 1 - sum(z[1:(N-1)])
        names(z) <- c(1:(N-1), Inf)
        return(z)
    }

##Two bpa objects by Bayes
combine.bpa.bs <- function(b1, b2)
    {
        ##Check for validity: Different Sample spaces!!!
        stopifnot(identical(b1$get.setlist() , b2$get.setlist()))

        lx <- log(b1$get.full.m())
        ly <- log(b2$get.full.m())
        xy <- exp(lx + ly)
        z <- xy / sum(xy)
        b <- bpa(n = b1$get.N(), setlist = b1$get.setlist(), mlist = z)
        return(b)
    }

##Two bpa objects by Dempster-Shafer
combine.bpa.ds <- function(b1, b2)
    {
        ##Check for validity: Different Sample spaces!!!
        stopifnot(identical(b1$get.setlist() , b2$get.setlist()))
        N <- b1$get.N()
        x <- b1$get.full.m()
        y <- b2$get.full.m()

        z <- tmp1 <- tmp2 <- numeric(N + 1)
        names(z) <- names(tmp1) <- names(tmp2) <- c(1:N, Inf)

        for(i in c(1:N,Inf))
            {
                tmp1 <- x[i] * y
                tmp2 <- x * y[i]
                z[i] <- z[i] + tmp1[i] + tmp1['Inf'] + tmp2['Inf']
            }
        z['Inf'] <- 1 - sum(z)
        
        b <- bpa(n = N, setlist = b1$get.setlist(), mlist = z)
        return(b)
    }

##List of bpa objects by Bayes
combine.bpa.list.bs <- function(blist)
    {
        k <- length(blist)
        y <- blist[[1]]
        for(i in 2:k)
            y <- combine.bpa.bs(y, blist[[i]])
        return(y)
    }

##List of bpa objects by Dempster-Shafer
combine.bpa.list.ds <- function(blist)
    {
        k <- length(blist)
        y <- blist[[1]]
        for(i in 2:k)
            y <- combine.bpa.ds(y, blist[[i]])
        return(y)
    }

##Two bpamat objects by Bayes
combine.bpamat.bs <- function(bmat1, bmat2)
    {
        ##Check for validity: Different Sample spaces!!!
        stopifnot(identical(bmat1$get.setlist() , bmat2$get.setlist()))
        ##Check for validity: Different Points
        stopifnot(identical(bmat1$get.pointlist() , bmat2$get.pointlist()))
        
        m <- bmat1$get.mat() * bmat2$get.mat()
        m1 <- m / rep(colSums(m), each = nrow(m))

        b <- bpamat(info = 'C', mat = m1)
        return(b)        
    }

##Two bpamat objects by Dempster-Shafer
combine.bpamat.ds <- function(bmat1, bmat2)
    {
        pts  <- bmat1$get.pointlist()
        sets <- bmat1$get.setlist()
        m <- matrix(0, nrow = length(sets), ncol = length(pts))
        row.names(m) <- sets
        colnames(m) <- pts
        ##Check for validity: Different Sample spaces!!!
        stopifnot(identical(sets , bmat2$get.setlist()))
        ##Check for validity: Different Points
        stopifnot(identical(pts , bmat2$get.pointlist()))

        for(i in pts)
            {
                p1 <- bmat1$get.point(i)
                p2 <- bmat2$get.point(i)
                m[,i] <- combine.ds(p1$get.full.m(), p2$get.full.m())
            }

        b <- bpamat(info = 'C', mat = m)
        return(b)        
    }

##List of bpamat object by Bayes
combine.bpamat.list.bs <- function(bmatlist)
    {
        k <- length(bmatlist)
        y <- bmatlist[[1]]
        for(i in 2:k)
            y <- combine.bpamat.bs(y, bmatlist[[i]])
        return(y)
    }

##List of bpamat objects by Dempster-Shafer
combine.bpamat.list.ds <- function(bmatlist)
    {
        k <- length(bmatlist)
        y <- bmatlist[[1]]
        for(i in 2:k)
            y <- combine.bpamat.ds(y, bmatlist[[i]])
        return(y)
    }
