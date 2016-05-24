xyplotPCArchetypes <- function(x, y, data.col = 1, data.pch = 19, data.bg = NULL, atypes.col = 2, atypes.pch = 19, 
                      ahull.show = FALSE, ahull.col = atypes.col,chull = NULL, chull.col = gray(0.7), 
                      chull.pch = 19, adata.show = FALSE, adata.col = 3, adata.pch = 13, link.col = data.col, 
                      link.lty = 1, ...){
 
 zs <- x
 data <- y
 
 plot(data, col = data.col, pch = data.pch, bg = data.bg, ...)
 points(zs, col = atypes.col, pch = atypes.pch, ...) #change respect to xyplot.

 if(!is.null(chull)){
  points(data[chull, ], col = chull.col, pch = chull.pch, ...)
  lines(data[c(chull, chull[1]), ], col = chull.col, ...)
  }

 if(ahull.show)
  lines(ahull(zs), col = ahull.col)
 if(adata.show){
  adata <- fitted(zs)
  link.col <- rep(link.col, length = nrow(adata))
  link.lty <- rep(link.lty, length = nrow(adata))
  points(adata, col = adata.col, pch = adata.pch, ...)
  for (i in seq_len(nrow(data)))
   lines(rbind(data[i, ], adata[i, ]), col = link.col[i], lty = link.lty[i], ...)
  }
 invisible(NULL)
}
