.plt.legend <-
function(colnms, horiz, stroke, fill, shape, usr) {

  par(xpd=NA) 
  
  ll <- legend(0,0, legend=colnms, cex=.7, pt.cex=0.9,
               horiz=TRUE, plot=FALSE)  # get coordinates

  if (horiz) {
    ytop <- usr[3] - (1.25 * ll$rect$h) 

    axis.horiz <- usr[2] - usr[1]
    lgnd.hhalf <- (ll$rect$w) / 2
    xleft <- usr[1] + axis.horiz/2 - lgnd.hhalf

    legend(xleft, ytop, legend=colnms, horiz=TRUE, box.lwd=.5, 
           box.col="gray30", cex=.7, pt.cex=1, pt.bg=fill,
           col=stroke, pch=shape)  # display legend
  }

  else {

    # just evaluate largest col name of y, as they are additive
    max.width <- which(nchar(colnms) == max(nchar(colnms)))
    ll <- legend(0,0, legend=max.width, cex=.7, pt.cex=0.9,
                 horiz, plot=FALSE)  # get coordinates

    size <- (par("cxy")/par("cin"))  # 1 inch in user coordinates 
    epsilon <- (size[1] - ll$rect$w) / 3

    axis.vert <- usr[4] - usr[3]
    xleft <- usr[2] + epsilon   # usr[2] user coordinate of right axis
    lgnd.vhalf <- (ll$rect$h) / 2
    axis.cntr <- axis.vert / 2  + usr[3]
    ytop <- axis.cntr + lgnd.vhalf  # user coordinate of legend top

    legend(xleft, ytop, legend=colnms, horiz=FALSE, box.lwd=.5, 
           box.col="gray30", cex=.7, pt.cex=1, pt.bg=fill,
           col=stroke, pch=shape)  # display legend
  }

}
