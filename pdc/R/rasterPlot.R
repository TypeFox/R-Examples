
rasterPlot <- function(cl, raw, monochrome=FALSE, aspect=.2, ...)
{
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  
  totalw <- aspect
	
  x <- cl
  X <- cl$data
  type <- "rectangle"
        par(fig = c(0, 0.3*totalw, 0, 1), mar = c(0, 0, 0, 0), new = F)
        class(x) <- "hclust"
        plot(as.dendrogram(x), horiz = TRUE, leaflab = "none", 
            frame.plot = FALSE, edgePar = list(col = "black", 
                lw = 2), type = type)
        scale <- 0.9
        offset <- 0.05
        for (i in 1:x$N) {
            par(fig = c(0.35*totalw, 1*totalw, offset + (i - 1)/x$N * scale, 
                offset + (i/x$N * scale)), mar = c(0, 0, 0, 0), 
                new = TRUE)
            if (x$multichannel) {
                dat <- X[, x$order[i], 1]
            }
            else {
                dat <- X[, x$order[i]]
            }
          
           myimg <- raw[[x$order[i]]]

        plot(c(0,3),c(0,1),type="n",xaxt='n',yaxt='n',ann=F, frame.plot=F)
		pimg <- raw[[x$order[i]]]
		if (monochrome) pimg <- sw(pimg)
        rasterImage(pimg,0,0,1,1)
		
  }

	invisible()
}