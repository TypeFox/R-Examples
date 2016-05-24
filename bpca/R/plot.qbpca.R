plot.qbpca <- function(x,
                       xlab='Index',
                       ylab='r',
                       pch=c(1,8),
                       col=c(4,2), ...)
{
  if (!inherits(x, 'qbpca'))
    stop("Use this function only with 'qbpca' class!")

  plot(x$obs,
       ylim=c(-1,1),
       xlab=xlab,
       ylab=ylab,
       pch=pch[1],
       col=col[1], ...)

  points(x$var.rb,
         col=col[2],
         pch=pch[2])

  legend('bottomleft',
         c('r', 'r.rb'),
         pch=pch,
         col=col,
         title='legend')
}
