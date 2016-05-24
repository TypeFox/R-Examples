#' 
#' @method plot errorLocation
#' @param x errorLocation object
#' @param topn Number of variables to show in 'errors per variable plot'. Only the top-n are are shown. 
#'      By default the top-20 variables with the most errors are shown.
#' @param ... other arguments that will be transferred to \code{barplot}
#' @rdname errorLocation
#' @export
plot.errorLocation <- function(x, topn=min(10,ncol(x$adapt)), ...){
  lgcrit <- 50
  oldpar <- par(mfcol=c(2,2))
  N <- nrow(x$adapt)
  if ( N <= 1 ) stop('Nothing to plot (single datapoint)')
  #########################################################################
  # PLOT ERRORS PER VARIABLE (1,1)
  #########################################################################
  vf <- colSums(x$adapt)
  # break variable names at 5th char.
  names(vf) <-  gsub('(.{5})','\\1\n',names(vf))

  varfreq <- sort(vf[vf>0]/N,decreasing=TRUE)
  barplot( varfreq[1:min(topn,length(varfreq))]
         , main = paste("Errors per variable",paste('(top ',topn,')',sep=''))
         , xlab = "Frequency"
         , horiz = TRUE
         , las = 1
         , ...
         )
  #########################################################################
  # PLOT ERRORS PER OBSERVATION (2,1)
  #########################################################################
   cnt <- table(rowSums(x$adapt))
   ner <- as.numeric(names(cnt))

   cnt <- as.numeric(cnt)
   noerr <- ner==0
   nnoer <-sum(cnt[noerr],0)
   cnt <- cnt[!noerr]
   ner <- ner[!noerr]
    
    lg <- ''
    if ( max(ner) > lgcrit ) lg <- paste(lg,'x',sep='')
    if ( max(cnt) > lgcrit ) lg <- paste(lg,'y',sep='')
  plot( ner, as.numeric(cnt),
         , main = "Errors per record"
         , ylab = "Count"
         , xlab = "Number of errors"
         , log=lg,
         , pch='.'
         , cex=5
         , ...
         )
    grid()
    mtext(paste(nnoer,'records with no errors'),side=3,line=0)
  #########################################################################
  # PLOT DURATION DENSITY AND CUMULATIVE DENSITY (1,2)
  #########################################################################
    # compute densisty in log-space
    eps <- sqrt(.Machine$double.eps)
    du <- density(log(x$status$user+eps))
    de <- density(log(x$status$elapsed+eps))
    du$F <- cumsum(du$y)*(du$x[2]-du$x[1])
    de$F <- cumsum(de$y)*(de$x[2]-de$x[1])
    du$x  <- exp(du$x)
    de$x  <- exp(de$x)
    par(mar=c(5,4,4,5)+.1)
    
    plot(de$x,de$y,
        log='x',
        type='l',
        lwd=2,
        col='black',
        xlab='Duration [seconds]',
        ylab='Density',
    )
    lines(du,lwd=2,col='blue')
    par(new=TRUE)
    
    plot(de$x,de$F,
        col='black',
        lwd=2,
        type='l',
        lty=2,
        log='x',
        xaxt="n",
        yaxt="n",
        xlab="",
        ylab="")
    axis(4)
    mtext("Cumulative",side=4,line=3)
    lines(du$x,du$F,lwd=2,lty=2,col='blue')
    grid()
    # default horiz abline used by plot.density crosses image border.
    lines(c(min(de$x),max(de$x)),c(0,0),col='gray',lwd=0.1)

    # 'legend': text colour corresponds to line colors
    mtext(paste('Elapsed time (',secToHuman(sum(x$status$elapsed)),')', sep=''),side=3,line=2,adj=0,font=2)
    mtext(paste('User time (',secToHuman(sum(x$status$user)),')', sep=''),side=3,line=1,col='blue',adj=0,font=2)
    mtext(paste(sum(x$status$maxDurationExceeded),'exceeded max. duration'),side=3,line=0)

  if(!("mip" %in% x$method)){
      #########################################################################
      # PLOT DEGENERACY
      #########################################################################
        tb <- table(x$status$degeneracy)
    #    br <- integer(max(x$status$degeneracy))
    #    br[as.numeric(names(tb))] <- tb
        lg <- ''
        dg <- as.integer(names(tb))
        ct <- as.integer(tb)
        if ( max(dg) > lgcrit ) lg <- paste(lg,'x',sep='')
        if ( max(ct) > lgcrit ) lg <- paste(lg,'y',sep='')
    
        plot(dg, ct,
            log=lg,
            main='Number of degenerate solutions',
            xlab='Degeneracy',
            ylab='Count',
            pch='.',
            cex=5,
        )
      grid()
  }
    # reset plot parameters
  par(oldpar)
}


# translate time in seconds to human-readable format
secToHuman <- function(x){
    if ( x < 1 )    return(paste(round(1000*x),'ms'))
    if ( x < 60 )   return(paste(round(x,2),'s'))
    if ( x < 3600 ) return(paste(round(x/60,2),'min'))
    if ( x < 24*3600 ) return(paste(round(x/3600,2),'hours'))
    if ( x < 24*3600*365 ) return(paste(round(x/3600/34,2),'days'))
    paste(round(x/24/3600/365,2),'years')
}



