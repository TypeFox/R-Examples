# plot.gibbs.A0.R

# Plots densities of the A(0) parameters for a B-SVAR model object.
# Also includes plotting of ther HDRs for each parameter

# 20120120 : Renamed function to conform to classing.


# This function from hdrcde package adds an hdr to a density figure.
# It is hidden in the pkg, so we include it here so it can be use to
# draw the hdrs.

# The next three functions are from Rob Hyndman's hdrcde package for
# computing HPDs for a density.  They are from version 2.07

"add.hdr" <- function(hdr, pos, width, col, horiz=FALSE, border=TRUE)
{
    ## Routine to add a single HDR on a graph
    ## Called by plot.hdr.conf().

    nint <- length(hdr[!is.na(hdr)])/2
    if(nint == 0)
        return(invisible())

    for(i in 1:nint)
    {
        l <- i*2-1  # lower
        tempx <- pos + c(-0.5,-0.5,0.5,0.5) * width
        tempy <- c(hdr[l],hdr[l+1],hdr[l+1],hdr[l])
        if(horiz==TRUE)
            polygon(tempy, tempx, col=col, border=border)
        else
            polygon(tempx, tempy, col = col, border=border)
    }
}

hdr <- function(x=NULL,prob=c(50,95,99),den=NULL,h=NULL,nn=5000,all.modes=FALSE)
{
    if(!is.null(x))
    {
        r <- diff(range(x))
        if(r==0)
            stop("Insufficient data")
    }
##     if(is.null(den))
##         den <- den.1d(x,h)
    alpha <- sort(1-prob/100)
    falpha <- calc.falpha(x,den,alpha,nn=nn)
    hdr.store <- matrix(NA,length(alpha),100)
    for(i in 1:length(alpha))
    {
        junk <- hdr.ends(den,falpha$falpha[i])$hdr
        if(length(junk) > 100)
        {
            junk <- junk[1:100]
            warning("Too many sub-intervals. Only the first 50 returned.")
        }
        hdr.store[i,] <- c(junk,rep(NA,100-length(junk)))
    }
    cj <- colSums(is.na(hdr.store))
    hdr.store <- matrix(hdr.store[,cj < nrow(hdr.store)],nrow=length(prob))
    rownames(hdr.store) <- paste(100*(1-alpha),"%",sep="")
    if(all.modes)
    {
        y <- c(0,den$y,0)
        n <- length(y)
        idx <- ((y[2:(n-1)] > y[1:(n-2)]) & (y[2:(n-1)] > y[3:n])) | (den$y == max(den$y))
        mode <- den$x[idx]
    }
    else
        mode <- falpha$mode
    return(list(hdr=hdr.store,mode=mode,falpha=falpha$falpha))
}

"calc.falpha" <-
function(x=NULL, den, alpha, nn=5000)
{
    # Calculates falpha needed to compute HDR of density den.
    # Also finds approximate mode.
    # Input: den = density on grid.
    #          x = independent observations on den
    #      alpha = level of HDR
    # Called by hdr.box and hdr.conf

    if(is.null(x))
        calc.falpha(x=sample(den$x, nn, replace=TRUE, prob=den$y), den, alpha)
    else
    {
        fx <- approx(den$x,den$y,xout=x,rule=2)$y
        falpha <- quantile(sort(fx), alpha)
        mode <- den$x[den$y==max(den$y)]
        return(list(falpha=falpha,mode=mode,fx=fx))
    }
}

"hdr.ends" <-
function(den,falpha)
{
    miss <- is.na(den$x) # | is.na(den$y)
    den$x <- den$x[!miss]
    den$y <- den$y[!miss]
    n <- length(den$x)
    if(falpha > max(den$y))
        return(list(falpha=falpha,hdr=NA) )
    dd <- den$y - falpha
    dd <- dd[2:n]*dd[1:(n-1)]
    index <- (1:(n-1))[dd<=0]
    index <- index[!is.na(index)]
    ni <- length(index)
    intercept <- numeric(ni)
    if(ni>0)
    {
        for(j in 1:ni)
        {
            idx <- c(index[j],index[j]+1)
            intercept[j] <- approx(den$y[idx],den$x[idx],xout=falpha)$y
        }
    }
    intercept <- sort(unique(intercept))
    ni <- length(intercept)
    if(ni == 0)
        intercept <- c(den$x[1],den$x[n])
    x1 <- 0.5*(intercept[1] + den$x[1])
    x2 <- 0.5*(intercept[ni] + den$x[n])
    if(approx(den$x,den$y,xout=x1)$y > falpha)
        intercept <- c(NA,intercept)
    if(approx(den$x,den$y,xout=x2)$y > falpha)
        intercept <- c(intercept,NA)
    return(list(falpha=falpha,hdr=intercept))
}


"plot.gibbs.A0" <- function(x, hpd=0.68, varnames=attr(x, "eqnames"), ...)
{
    # Get constants
    m <- ncol(x$ident)
    ident <- t(x$ident)
    prob <- hpd*100
    # Convert to a matrix
    x <- A02mcmc(x)

    # Set location in ident counter
    k <- 1

    # Plot setup
    par(mar=c(2,2,1,1))
    split.screen(c(m,m))

    # Loop over the elements of the ident matrix of free parameters,
    # grab the correct part of the A(0) posterior and plot the density
    # in the matrix of densities.

    for(i in 1:m)
    { for (j in 1:m)
      {
      # Both non-zero
      if(ident[i,j]==1)
      {
          # Set the screen
          screen((i-1)*m + j)

          # Find the densities
          den1 <- density(x[,k])

          # Compute hdr's for densities
          hdr1 <- hdr(x[,k], prob=prob, den=den1)

          # Put larger margins on the left and top of plots in figure
          if(i==1 | j==1) par(omi=c(0.15, 0.5, 0.5, 0.15))

          # Actual plotting
          plot(den1, lty=1, main="", ylab="", cex.axis=0.75)
          abline(v=0)

          # Plot the first HDR
          nregions <- nrow(hdr1$hdr)
          maxden1 <- max(den1$y)
          for(l in 1:nregions)
          {
              lines(range(den1$x), rep(hdr1$falpha[l],2), lty=2)
              for(n in 1:length(hdr1$hdr[l,]))
                  lines(rep(hdr1$hdr[l,n],2),c((0.01+(l-1)*0.02)*maxden1,hdr1$falpha[l]), lty=2)
          }
          for(l in 1:nrow(hdr1$hdr))
              add.hdr(hdr1$hdr[l,], (0.01+(l-1)*0.02)*maxden1,
                      0.1*maxden1, col="black", horiz=TRUE, border=FALSE)

          # Label figure as necessary
          if(i==1) mtext(varnames[j], side=3, line=1, outer=FALSE)
          if(j==1) mtext(varnames[i], side=2, line=3, outer=FALSE)

          k <- k+1
      }

      # Handle cases where there is no free parameter to plot, but we
      # need to label the row or column
      if(ident[i,j]==0 & (i==1 || j==1))
      {
          screen((i-1)*m + j)
          plot.new()
          if(i==1) mtext(varnames[j], side=3, line=1, outer=FALSE)
          if(j==1) mtext(varnames[i], side=2, line=3, outer=FALSE)


      }
  }
  }
    close.screen(all.screens = TRUE)
}
