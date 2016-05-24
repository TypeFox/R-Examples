setClass('longCat',
         representation(   data        = "matrix",
                           data.sorted = "matrix",
                           dim         = "integer",
                           times       = "matrix",
                           endt        = "matrix",
                           times.sorted= "matrix",
                           endt.sorted = "matrix",
                           labels      = "character",
                           factors     = "numeric",
                           IndTime     = "logical",
                           nfactors    = "integer",
                           sorted      = "logical",
                           ascending   = "logical",			   
                           group       = "matrix",
                           groupLabels = "character",
                           order.data  = "matrix") 
)
setMethod("summary",
          signature(object = "longCat"),
          definition = function (object, ...) 
          {
            if( object$sorted) temp <- object[c(3,8:15)]
            if(!object$sorted) temp <- object[c(3,8:12 )]
            temp$group <- table(temp$group)
            print(temp)
          }
)


longContPlot <- function(y, times=NULL, jog=FALSE, ylim=NULL, xlim=NULL, ...)
{
  # check the inputs and set graphing parameters
  if( is.null(ylim) ){ ylim=range(y, na.rm=T) }
  if( is.null(times) )
  { 
    txx <- 1:ncol(y)
    xlim=range(txx) 
    times <- data.frame( matrix(txx, 1, ncol(y)) )
  }
  if( !is.null(times) & is.null(dim(times)) )
  { 
    times <- data.frame( matrix(times, 1, ncol(y)) ) 
  }
  if( all( dim(y)==dim(times) ) & !is.null(times) ){ txx <- times[1,] }
  if( any( dim(y)!=dim(times) ) & !is.null(times) ){ txx <- times[1,] }
  if( !is.null(times) & is.null(xlim) ){ xlim=range(times) }
  if(  is.null(times) & is.null(xlim) ){ xlim=c(1,ncol(y)) }
  
  if(jog)
  {
    j <- matrix( runif(nrow(y), -.25, .25), nrow(y), ncol(y) )
    y <- y + j 
  }
  
  # initiate a blank plot
  plot( unlist(txx), unlist(y[1,]), col='transparent', ylim=ylim, xlim=xlim, type='n', ...)
  
  # loop through subjects adding them to the plot
  if( any( dim(y)!=dim(times) ) )
  {
    for(r in 1:nrow(y)){ lines(txx, y[r,])  }
  }
  if( all( dim(y)==dim(times) ) ) 
  {
    for(r in 1:nrow(y))
    { 
      txx <- times[r,]
      lines(txx, y[r,])  
    }
  }
}

lunique <- function(y)
{
  u <- unique(y)
  sum(!is.na(u))
}
levelCheck <- function(y)
{
  lu <- apply(y, 2, lunique)
  maxlu <- max(lu)
  if( maxlu > 9 )
  {
    warning('One or more variables in y has 10 or more categories.\nConsider using continuous methods\nsuch as longContPlot().')
  }
  # determine unique values
  factors <- as.numeric(levels(factor(unlist(y))))
  return( factors )
}
dimCheck <- function(y, times)
{
  dimErr1 <- "The dimension of times does not equal the dimension of y"
  dimErr2 <- "The length of times does not equal ncol(y)"
  ry <- nrow(y)
  cy <- ncol(y)
  rt <- nrow(times)
  ct <- ncol(times)
  
  if( rt>1 & rt!=ry ) stop(dimErr1) 
  if( cy!=ct ) stop(dimErr2)
  if( ry==rt & cy==ct )
  {
    IndTime <- T
  }
  if( rt==1 & cy==ct )
  {
    IndTime <- F
  }
  return(IndTime)
}


longCat <- function(y, times=NULL, Labels=NULL, tLabels=NULL, id=NULL, endt=NULL)
{
  # convert data frame to matrix and character to numeric
  y <- apply(y, 2, as.numeric)
  # check id
  if( !is.null(id) & nrow(y) != length(id) )
  {
    stop('The number of IDs length(id) does not match the number of rows in y nrow(y)')
  }
  # create the order.data matrix
  order.data <- matrix(NA, nrow(y), 3)
  if( !is.null(id) ) order.data[,1] <- id
  # rescale inputs to positive sequential integers from 1 to the maximum 
  # number of categories
  u <- unique(c(y))
  u <- u[ !is.na(u) ]
  u <- u[order(u)]
  if( !all(u == 1:length(u)) )
  {
    temp <- y
    for(i in 1:length(u))
    {
      temp[y==u[i]] <- i
    }
    y <- temp; rm(temp)
  }
  rm(u)
  
  # check the levels using levelCheck()
  factors <- levelCheck(y)
  
  # check the times input and force to data.frame
  if( is.null(times) )
  { 
    times <- data.frame( matrix(1:ncol(y), 1, ncol(y)) ) 
  }
  if(!is.null(times))
  { 
    # if times is a vector and endt is NULL
    if( is.null(dim(times)) & is.null(endt) )
    {
      times <- data.frame( matrix(times, 1, ncol(y)) )
    }
    # if times is a matrix
    if(!is.null(dim(times))) 
    {
      if( ncol(times) == ncol(y) & is.null(endt) )
      {
        times <- data.frame( times )
      }
      if( ncol(times) == (ncol(y)+1) & is.null(endt) )
      {
        endt  <- times[,ncol(times)]
        times <- data.frame( times[,1:ncol(y)] )
      }
      if( ncol(times) == (ncol(y)+1) & !is.null(endt) & length(endt)==nrow(y) )
      {
        if( all( times[,ncol(times)]==endt ) )
        {
          endt  <- times[,ncol(times)]
          times <- data.frame( times[,1:ncol(y)] )
        }
        if( all( times[,ncol(times)]==endt ) )
        {
          warning( paste("times has an extra column and endt is provided.\n",
                         "Make sure this is intended.\n",
                         "The last column of times is ignored.") )
          times <- data.frame( times[,1:ncol(y)] )
        }      
      }
    }
  }
  if( !is.null(endt) & length(endt)==1 )
  {
    endt <- as.matrix( rep(endt, nrow(y) ) )
  }
  if( is.null(endt) )
  {
    avgt <- ( max(times, na.rm=T) - min(times, na.rm=T) ) / ncol(y)
    endt <- as.matrix( rep(avgt, nrow(y) ) )
  }
  endt <- as.matrix( endt )
  
  # check the dimensions of the data (y) and the times input
  IndTime <- dimCheck(y, times)
  
  # create Labels if not provided
  if(is.null(Labels)){ Labels <- factors }
  
  # count the number of factors returned by leveCheck()
  nfactors <- length(factors)
  
  # check that Labels, nfactors, tLabels, and times conform
  if( length(Labels) != nfactors )
  {
    warning(paste('The number of labels in Labels does not equal the\n',
                  'number of unique values in the data set.'))
  }
  if( !is.null(tLabels) & length(unique(times)) != length(tLabels)  )
  {
    warning(paste('The number of labels in tLabels does not equal the\n',
                  'number of unique values in times.'))
  }
  
  # if individually varying times of observation, make sure times is a matrix
  if(!is.null(times) & IndTime) times <- as.matrix(times)
  
  # apply class and return output
  lc =  list(data=y,
             data.sorted=NULL,
             dim=dim(y),
             times=times,
             endt=endt,
             times.sorted=NULL,
             endt.sorted=NULL,
             labels=Labels,
             tLabels=tLabels,
             factors=factors,
             IndTime=IndTime,
             nfactors=nfactors,
             sorted=FALSE,
             ascending = NULL,					   
             group = NULL,
             groupLabels = NULL,
             order.data = order.data)
  class(lc) = 'longCat'
  return(lc) 
}


norpt <- function( alist = c(1,2,2,3,3,3,4,4,4,4,5) )
{
  outlist <- alist[1]
  for(i in 2:length(alist))
  {
    if( !is.na(alist[i-1]) & !is.na(alist[i]) )
    {
      if(alist[i-1] != alist[i]){ outlist <- c(outlist, alist[i]) }
    }
  }
  outlist <- c(outlist, rep(NA, (length(alist)-length(outlist))))
  outlist
}
makePatterns <- function(dat, times=NULL, num=TRUE, mindur=NULL, igrpt=FALSE)
{
  # set times if null
  if( is.null(times) ) times <- 1:ncol(dat)
  # reduce the effect of short durations
  if(!is.null(mindur) & length(dim(times))==2)
  {
    times <- times - times[,1]
    mintime <- times <= mindur & !is.na(times)
    dat[mintime] <- NA
  }
  # if desired, (ig)nore (r)e(p)ea(t) observations, 
  #   i.e., c(1, 2, 2, 3) becomes c(1, 2, 3)
  if(igrpt){ dat <- t( apply(dat, 1, norpt) )  }
  # concatenate rows into a string
  out <- apply(dat, 1, paste, collapse="")  
  # if desired, turn to numeric and reduce the scale
  if(num)
  {
    g    <- gsub('NA', '', out)
    nc   <- nchar(g)
    tens <- 10^nc
    out  <- ( as.numeric(g)/tens )*10
  }
  as.matrix(out, nrow(dat), 1)
}
sorter <- function(lc, ascending=TRUE, whichColumns=NULL, num=TRUE, mindur=NULL, 
                   igrpt=FALSE, customSort=NULL, initFirst=FALSE, group=NULL, groupLabels=NULL, ggap=10)
{
  # if object is already sorted, don't attempt to resort
  if(lc$sorted)
  {
    sortWarn <- paste('lc object already sorted, use the following approach\n',
                      'lc <- (y)\n',
                      'lc_sort1 <- sorter(lc, group=g1)\n',
                      'lc_sort2 <- sorter(lc, group=g2)\n',
                      'lc_customSort <- sorter(lc, customSort=cs)\n'
    )
    stop(sortWarn)    
  }
  
  # if customSort contains missing data, stop
  if( !is.null(customSort) )
  {
    if( any(is.na(customSort)) )
    {
      stp <- paste('\ncustomSort cannot contain missing values \n',
                   'consider trying one of:\n',
                   '1. Subsetting the data first\n',
                   '2. Setting missing to an extreme value and adding a group variable, such as\n',
                   '     lc <- longCat(y)\n',
                   '     customSort[is.na(customSort)] <- -99\n',
                   '     group <- customSort==-99\n',
                   "     lc_customSort <- sorter(lc, group=group, \n",
                   "          groupLabels=c('nonMissing customSort', 'Missing customSort'),\n",
                   "          customSort=customSort\n")
      stop(stp)
    }  
  }
  # check for missing values on the grouping variable
  if( !is.null(group) )
  {
    if( any(is.na(group)) ) 
    {
      ndeleted <- sum(is.na(group))
      nag <- !is.na(group)
      lc$data <- lc$data[nag,]
      lc$order.data <- lc$order.data[nag,]
      if(lc$IndTime) lc$times <-  lc$times[nag,]
      customSort <- customSort[nag]
      group <- group[nag]
      w <- paste(ndeleted, 
                 ' row(s) with missing membership on the group variable\n',
                 'have been deleted.\n\n', 
                 'If a large number of cases have missing data on the\n',
                 'group variable, consider recoding the missings into\n',
                 'their own group, e.g.:\n\n',
                 '     group[is.na(group)] <- -999 \n\n',
                 'and add a missing label to groupLabels, e.g.:\n\n',
                 "     groupLabels=c('Missing', 'Group1', 'Group2', 'Etc.')\n",
                 sep='')
      warning(w)
    }
  }
  
  # check inputs and set additional sorting parameters
  if(is.null(whichColumns)) whichColumns = 1:ncol(lc$data)
  if( lc$IndTime) pats <- makePatterns(lc$data[,whichColumns], lc$times[,whichColumns], num, mindur, igrpt)
  if(!lc$IndTime) pats <- makePatterns(lc$data[,whichColumns], NULL, num, mindur, igrpt)
  if(lc$IndTime) tpat <- do.call(order, data.frame(lc$times[,whichColumns]) )
  if( is.null(group) & !lc$IndTime) o <- order(pats, decreasing = !ascending)
  if(!is.null(group) & !lc$IndTime) o <- order(group, pats, decreasing = !ascending)
  if( is.null(group) &  lc$IndTime) o <- order(pats, tpat, decreasing = !ascending)
  if(!is.null(group) &  lc$IndTime) o <- order(group, pats, tpat, decreasing = !ascending)
  # if there is a custom sorting variable
  if(!is.null(customSort) & !initFirst)
  {
    if( is.null(group) & !lc$IndTime) o <- order(customSort, pats, decreasing = !ascending)
    if(!is.null(group) & !lc$IndTime) o <- order(group, customSort, pats, decreasing = !ascending)
    if( is.null(group) &  lc$IndTime) o <- order(customSort, pats, tpat, decreasing = !ascending)
    if(!is.null(group) &  lc$IndTime) o <- order(group, customSort, pats, tpat, decreasing = !ascending) 
  }
  if(!is.null(customSort) &  initFirst)
  {
    if( is.null(group) & !lc$IndTime) o <- order(lc$data[,1], customSort, pats, decreasing = !ascending)
    if(!is.null(group) & !lc$IndTime) o <- order(lc$data[,1], group, customSort, pats, decreasing = !ascending)
    if( is.null(group) &  lc$IndTime) o <- order(lc$data[,1], customSort, pats, tpat, decreasing = !ascending)
    if(!is.null(group) &  lc$IndTime) o <- order(lc$data[,1], group, customSort, pats, tpat, decreasing = !ascending) 
  }
  data.sorted <- lc$data[o,]
  group <- group[o]
  if(lc$IndTime) times.sorted <- lc$times[o,]
  endt.sorted <- lc$endt[o]
  lc$order.data[,2] <- o
  lc$order.data[,3] <- pats
  
  # check grouping parameters
  if( !is.null(group) )
  {
    if( nrow(as.matrix(group)) != nrow(lc$data) )
    {
      stop('group has a length that does not equal the number of rows in y')
    }
    group <- as.numeric(group)
    u <- unique(group)
    if( is.null(groupLabels) ) groupLabels <- paste('Group',u,sep='')
    if( length( u ) != length(groupLabels) )
    {
      stop(paste('The number of labels in groupLabels does not equal the\n',
                 'number of unique values in the variable group.'))
    }
  }
  # if grouping/stratification is present, augment the data with empty rows
  # which will visually dilineate groups
  if( !is.null(group) )
  {
    temp <- vector('list', length(u) )
    gtemp <- vector('list', length(u) )
    ttemp <- vector('list', length(u) )
    etemp <- vector('list', length(u) )
    blines <- matrix(NA, ggap, ncol(data.sorted))
    glines <- blines[,1]
    for(i in 1:length(u))
    {
      gdat <- data.sorted[group==u[i],]
      pats <- makePatterns(gdat[,whichColumns], NULL, num, mindur)
      o    <- order(pats, decreasing = !ascending)
      if(lc$IndTime)
      {
        tdat <- times.sorted[group==u[i],]
        edat <- endt.sorted[group==u[i]]
        pats <- makePatterns(gdat[,whichColumns], tdat[,whichColumns], num, mindur)
        tpat <- do.call(order, data.frame(tdat) )
        o    <- order(pats, tpat, decreasing = !ascending)
      }     
      temp[[i]]  <- rbind(blines, gdat[o,]) 
      gtemp[[i]] <- as.matrix(c(glines, group[group==u[i]]))
      if(lc$IndTime) 
      {
        ttemp[[i]] <- rbind(blines    , tdat[o,] )
        etemp[[i]] <-     c(blines[,1], edat[o ] )
      }
    }
    data.sorted <- do.call(rbind, temp); rm(temp)
    group <- do.call(rbind, gtemp); rm(gtemp)
    if(lc$IndTime){ times.sorted <- do.call(rbind, ttemp); rm(ttemp) }
    endt.sorted <- unlist(etemp); rm(etemp)
  }                                          
  
  # assign NULL if not previously assigned
  if(!lc$IndTime) times.sorted = NULL
  
  # return modified lc object
  lc = list( data = lc$data,
             data.sorted = data.sorted,
             dim = lc$dim,
             times = lc$times,
             endt = lc$endt,
             times.sorted = times.sorted,
             endt.sorted = endt.sorted,
             labels = lc$labels,
             tLabels = lc$tLabels,
             factors = lc$factors,
             IndTime = lc$IndTime,
             nfactors = lc$nfactors,
             sorted = TRUE,
             ascending = ascending,			   
             group = group,
             groupLabels = groupLabels, 
             order.data = lc$order.data)
  class(lc) = 'longCat'
  return(lc)        
}


colChoose <- function(colScheme, nfactors, reverse=FALSE)
{
  if(colScheme=='gray'){cols <- gray((nfactors-1):0/nfactors)}
  if(colScheme=='rainbow'){cols <- rainbow(nfactors)}
  if(colScheme=='heat'){cols <- heat.colors(nfactors, alpha = 1)[nfactors:1]}
  if(colScheme=='terrain'){cols <- terrain.colors(nfactors, alpha = 1)}
  if(colScheme=='topo'){cols <- topo.colors(nfactors, alpha = 1)}
  if(colScheme=='cm'){cols <- cm.colors(nfactors, alpha = 1)}  
  # some finessing to make sure contrast is sufficient when nFactors < 9
  if(reverse) cols <- cols[nfactors:1]
  return(cols)
}

longCatPlot <- function(lc, xlab="Days",
                        ylab=NULL, cols=NULL,
                        colScheme='heat', reverse=FALSE, lwd=.5, lcex=1, llwd=3, 
                        legendBuffer=.12, groupBuffer=0, groupRotation=90, gcex=1, 
                        seg.len=1, xlas=0, xcex=1, ...)
{
  if( is(lc) != 'longCat' ){ stop('longCatPlot requires an object of class longCat.')  }
  if(is.null(cols)){ cols <- colChoose(colScheme, lc$nfactors, reverse) }
  if(legendBuffer < 0 | legendBuffer > 1){stop('legendBuffer must be in [0,1]')}
  if(groupBuffer < 0 | groupBuffer > 1){stop('groupBuffer must be in [0,1]')}
  if(is.null(ylab)) ylab = paste("Each Line Represents a Participant, N =", lc$dim[1])
  
  # on the fly sorting
  if( !lc$sorted & ( any(is.na(lc$data)) |  lc$IndTime) ) lc <- sorter(lc, num=TRUE)
  if( !lc$sorted &  !any(is.na(lc$data)) & !lc$IndTime  ) lc <- sorter(lc, num=FALSE)
  
  # set up plot, pre-allocating a region for labels using ymax
  lo = min(lc$times, na.rm=T)
  up = max(lc$times, na.rm=T)
  xrange = lo:up
  
  # reps is used to automatically scale the x-axis
  reps = nrow(lc$data.sorted)-length(xrange)    
  # fix reps if negative, will occur when the number of cases is fewer than
  # the number of time points. This happens when plotting an individual
  if( reps < 0 ) reps <- 0
  
  # set additional plotting parameters
  xbuffer <- .25*mean(xrange[2:length(xrange)]-xrange[1:(length(xrange)-1)])
  groupBuffer <- ceiling( groupBuffer*up )
  tx <- c(lo-.5, xrange, rep( lo, reps), up+.5 )
  if( !is.null(lc$group) ){ tx[1] <- tx[1] - groupBuffer*xbuffer }
  ymax <- nrow(lc$data.sorted) + ceiling( legendBuffer*nrow(lc$data.sorted) )
  
  # initiate the empty plot
  plot(tx,y=rep(NA,length(tx)),col='transparent',ylim=c(0,ymax), 
       xlab=xlab, ylab='', axes=FALSE, ...)
  if(groupBuffer >0 | !is.null(lc$group)) title(ylab=ylab, mgp=c(1   ,1,0))
  if(groupBuffer==0 &  is.null(lc$group)) title(ylab=ylab, mgp=c(0.25,1,0))  
  
  # add axes
  if(!is.null(lc$tLabels)) axis( 1, at = unique(lc$times), 
                                 labels=lc$tLabels, las=xlas, cex.axis=xcex )
  if( is.null(lc$tLabels)) axis( 1, at = NULL )
  
  # plot loops
  for(r in 1:nrow(lc$data.sorted)) # loop over cases
  {
    # select plotting data for the r^th case
    pdat <- lc$data.sorted[r,]
    if( lc$IndTime & !lc$sorted) txx <- lc$times[r,]
    if( lc$IndTime &  lc$sorted) txx <- lc$times.sorted[r,]
    if(!lc$IndTime) txx <- lc$times
    tempy <- rep(r,2)
    for(j in 1:length(pdat)) # loop over observations
    {
      if( !is.na( pdat[j] ) ) # if observation is NA, skip
      {
        # define the x-values for any but the last segment
        if( j <length(pdat) )
        { 
          tempx <- c(txx[j], txx[j+1])
          # correct for missing endpoint
          if(is.na(txx[j+1])) tempx[2] <- txx[j]+xbuffer 
        }
        # define the x-values for the last segment
        if( j==length(pdat) ){ tempx <- c(txx[j], txx[j]+lc$endt[j]) }
        # horizontal line plot
        lines(tempx, tempy, lwd=lwd, col=cols[ unlist(pdat[j]) ] )
      }
    }
  }
  
  # add legend at top 
  if(legendBuffer > 0)
  {
    legend(min(lc$times, na.rm=T), ymax + .1, legend=lc$labels, lty=1,
           cex=lcex, col=cols, bty='n', lwd=llwd, seg.len=seg.len, horiz=T)  
  }
  
  # if there is grouping add group labels
  if( !is.null(lc$group) )
  {
    g <- cbind(lc$group, 1:length(lc$group))
    gag <- aggregate(g[,2] ~ g[,1], g, mean, na.rm=T)
    u <- unique(lc$group)
    for(i in 1:length(u))
    {
      text(tx[1], gag[i,2], labels = lc$groupLabels[i], srt=groupRotation, cex=gcex)
    }
  }
  # return colors
  invisible(cols)
}
