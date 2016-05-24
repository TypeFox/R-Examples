#---------------------------------------------------------------------------
#
#
#Author...									Date: 23-Feb-2012
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#


#================================================================================
#  1. method for class montePop...
#
setMethod('hist',
          signature(x = 'montePop'),
function(x,
         col = 'gray90',
         ...
        )
{
    hg = hist(x@popVals, col=col, ...)

    return(invisible(hg))

}    #hist for ''montePop'
) #setMethod




#================================================================================
#  2. method for class monteSample...
#
setMethod('hist',
          signature(x = 'monteSample'),
function(x,
         n = NA,
         main = NULL,
         col = 'gray90',
         ...
        )
{
#------------------------------------------------------------------------------
#   histogram(s) for sample size(s) n...
#------------------------------------------------------------------------------
#
    if(is.null(n) || all(is.na(n)))
      n = x@n
    else {
      n = sort(round(n))  #integer sample sizes only
      n = intersect(n, x@n)
      if(length(n) == 0)
        stop('Please choose sample size(s), n, that are in the monteSample object!')
    }
    n.names = paste('n',n,sep='.')
    

    nn = length(n)
    if(nn > 1) {
      snn = sqrt(nn)
      nrow = ceiling(snn)
      ncol = round(snn)
      opar = par(mfrow=c(nrow, ncol))
    }

#
#   determine the overall break points and maximum in y so each histo is scaled the same...
#
    means = x@means
    breaks = hist(as.vector(as.matrix(means)), plot=FALSE)$breaks #for full range histo
    maxCounts = rep(0, nn)
    for(i in seq_len(nn)) {
      n.i = n.names[i]
      maxCounts[i] = max(hist(means[,n.i], breaks = breaks, plot=FALSE)$counts)
    }
    y.max = max(maxCounts)

#
#   now we can plot them...
#
    hg = vector('list',nn)
    for(i in seq_len(nn)) {
      n.i = n.names[i]
      main = n.names[i]
      hg[[i]] = hist(means[,n.i], breaks = breaks, ylim=c(0,y.max*1.1), col=col,
                     main = main, ...)
    }

    if(nn > 1)
      par(opar)
         
    return(invisible(hg))

}    #hist for ''monteSample'
) #setMethod





#================================================================================
#  3. method for class monte...
#
setMethod('hist',
          signature(x = 'monte'),
function(x,
         n = NA,
         xlab = '',
         col = 'gray90',
         type = c('normal', 'bootstrap', 'population'),
         ...
        )
{
#------------------------------------------------------------------------------
#   histogram(s) for sample size(s) n; doesn't plot anything, just passes along...
#------------------------------------------------------------------------------
#
    type = match.arg(type)

    hg = NULL
    if(nchar(xlab) == 0)
      xlab = x@estimate
    if(type == 'normal' && !is.null(x@NTsamples)) 
      hg = hist(x@NTsamples, n=n, col=col, xlab=xlab, ...)
    else if(type == 'bootstrap' && !is.null(x@BSsamples))
      hg = hist(x@BSsamples, n=n, col=col, xlab=xlab, ...)
    else if(type == 'population')
      hg = hist(x@pop, col=col, xlab=xlab, ...)
    else
      cat('\nDesired histogram information (',type,') is not availble in this object\n', sep='')

    return(invisible(hg))

}    #hist for ''monte'
) #setMethod
 
