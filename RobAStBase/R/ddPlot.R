setMethod("ddPlot", signature = signature(data = "matrix"),
          function(data, dist.x = NormType(), dist.y  = NormType(),
       cutoff.x, cutoff.y, ...,
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, lab.pts, adj, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff, text.abline = TRUE,
       text.abline.x = NULL, text.abline.y = NULL,
       cex.abline = par("cex"), col.abline = col.cutoff,
       font.abline = par("font"), adj.abline = c(0,0),
       text.abline.x.x = NULL, text.abline.x.y = NULL,
       text.abline.y.x = NULL, text.abline.y.y = NULL,
       text.abline.x.fmt.cx = "%7.2f", text.abline.x.fmt.qx = "%4.2f%%",
       text.abline.y.fmt.cy = "%7.2f", text.abline.y.fmt.qy = "%4.2f%%",
       jitt.fac){
       mc <- as.list(match.call(expand.dots = TRUE, 
                                call = sys.call(sys.parent(1)))[-1])
       mc$data <- data
       do.call(RobAStBase:::.ddPlot.MatNtNtCoCo, args = mc)
})

setMethod("ddPlot", signature = signature(data = "data.frame"),
          function(data, dist.x = NormType(), dist.y  = NormType(),
       cutoff.x, cutoff.y, ...,
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, lab.pts, adj, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff, text.abline = TRUE,
       text.abline.x = NULL, text.abline.y = NULL,
       cex.abline = par("cex"), col.abline = col.cutoff,
       font.abline = par("font"), adj.abline = c(0,0),
       text.abline.x.x = NULL, text.abline.x.y = NULL,
       text.abline.y.x = NULL, text.abline.y.y = NULL,
       text.abline.x.fmt.cx = "%7.2f", text.abline.x.fmt.qx = "%4.2f%%",
       text.abline.y.fmt.cy = "%7.2f", text.abline.y.fmt.qy = "%4.2f%%",
       jitt.fac){

         mc <- match.call(call = sys.call(sys.parent(1)))
         mc$data <- t(as.matrix(data))
         do.call(ddPlot,args = as.list(mc[-1]))
       })

setMethod("ddPlot", signature = signature(data = "numeric"),
          function(data, dist.x = NormType(), dist.y  = NormType(),
       cutoff.x, cutoff.y, ...,
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, lab.pts, adj, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff,
       text.abline = TRUE,
       text.abline.x = NULL, text.abline.y = NULL,
       cex.abline = par("cex"), col.abline = col.cutoff,
       font.abline = par("font"), adj.abline = c(0,0),
       text.abline.x.x = NULL, text.abline.x.y = NULL,
       text.abline.y.x = NULL, text.abline.y.y = NULL,
       text.abline.x.fmt.cx = "%7.2f", text.abline.x.fmt.qx = "%4.2f%%",
       text.abline.y.fmt.cy = "%7.2f", text.abline.y.fmt.qy = "%4.2f%%",
       jitt.fac){

         mc <- match.call(call = sys.call(sys.parent(1)))
         mc$data <- matrix(data,nrow=1)
         do.call(ddPlot,args = as.list(mc[-1]))
       })
