##
## S3 method to plot 'TukeyC' object
##

plot.TukeyC <- function(x,
                        result=TRUE,
                        replicates=TRUE, 
                        pch=19,
                        col=NULL,
                        xlab=NULL,
                        ylab=NULL,
                        xlim=NULL,
                        ylim=NULL,
                        id.lab=NULL,
                        id.las=1,
                        rl=TRUE,
                        rl.lty=3,
                        rl.col='gray',
                        mm=TRUE,
                        mm.lty=1,
                        title='', ...)
{
  fun <- function(m) {
    a <- rep('\n', length(m))
    a[which(m != '')[1]] <- ''
    return(paste(a, m, sep=''))
  }

  if(!inherits(x,
               'TukeyC'))
    stop("Use only with 'TukeyC' objects!")

  if(is.null(xlab)) xlab <- 'Levels'

  if(is.null(ylab)) ylab <- 'Means'

  means <- x$Means[, 1]

  if(replicates)
    r <- x$Replicates

  groups <- 1:length(means)

  m.res <- t(x$Result[, 2:ncol(x$Result)])

  if(dim(m.res)[1] != 1) {
    m.res <- apply(m.res, 2, fun)
    id.groups <- c(apply(m.res,
                         2,
                         paste,
                         collapse=''))
  }
  else
    id.groups <- m.res 

  minmax <- x$Means[, 2:3]

  if(is.null(col))
    col <- 'black'

  if(is.null(xlim))
    xlim <- c(1,
              max(groups))

  if(is.null(ylim))
    if(replicates)
      ylim <- c(.90 * min(minmax[, 1]),
                max(minmax[, 2]))
    else
      ylim <- c(min(minmax[, 1]),
                max(minmax[, 2]))

  old.par <- par(mai=c(1, 1, 1.25, 1))

  plot(groups,
       means,
       pch=pch,
       col=col,
       xlab=xlab,
       ylab=ylab,
       xlim=xlim,
       ylim=ylim,
       axes=FALSE, ...)

  if(rl == TRUE)       
    segments(rep(-0.5,
                 length(means)),
             means,
             groups,
             means,
             lty=rl.lty,
             col=rl.col, ...) 

  if(mm == TRUE)
    segments(groups,
             minmax[, 2],
             groups,
             minmax[, 1],
             lty=mm.lty,
             col=col, ...)

  axis(2,
       at=round(seq(ylim[1],
                    ylim[2],
                    length.out=5),
                1))

  if(is.null(id.lab))
    id.lab <- rownames(x$Means)

  axis(1,
       at=1:length(means),
       labels=id.lab,
       las=id.las,
       col.axis=FALSE, ...)    

  if(result) 
    axis(3,
         at=1:length(means),
         labels=id.groups, ...)

  if(replicates)
    text(x=1:length(means),
         y=min(ylim),
         labels=r,
         pos=3, ...)

  mtext(text=id.lab,
        side=1,
        line=1,
        at=1:length(means),
        las=id.las, ...)

  title(title, ...)

  par(old.par)
}

