### This files contains functions for ploting density map for clusters
### Written: Wei-Chen Chen on 2008/11/04.


plotlikeprob <- function(DA, L, file = NULL,
    file.type = c("jpeg", "pdf"), ...){
  xlim <- range(DA[, 1])
  ylim <- range(DA[, 2])

  if(!is.null(file)){
#    bitmap(file = file, type = "jpeg", width = 10, height = 11)
    if(file.type[1] == "pdf"){
      pdf(file = file, width = 6, height = 8)
    } else{
      jpeg(filename = file, width = 600, height = 800, quality = 100)
    }
  }

  layout(matrix(c(1,2), nrow = 2, ncol = 1, byrow = TRUE), , c(3, 1))

  L.index <- order(L, decreasing = TRUE)
  L <- L[L.index]
  DA <- DA[L.index,]

  L.col <- create.colors(L, n = 10)
  L.cex <- create.cex(sqrt(L), n = 10)
  L.col.cl <- match.colcex(L, L.col)
  L.cex.cl <- match.colcex(L, L.cex)

  plot(NULL, NULL, xlim = xlim, ylim = ylim, type = "n",
       xlab = "X", ylab = "Y", ...)
  points(DA, col = L.col.cl, pch = 16, cex = L.cex.cl, ...)

  plotbar(L.col, L)

  if(!is.null(file)){
    dev.off()
  }
}


plotbar <- function(col, cluster, main = NULL){
  z <- matrix(seq(along = col), ncol = 1)

  ticket <- NULL
  ticket$at <- seq(0, 1, length = min(6, length(unique(cluster))))
  ticket$labels <- c(min(cluster),
                     min(cluster) + ticket$at[c(-1, -length(ticket$at))] * diff(range(cluster)),
                     max(cluster))
  ticket$labels <- signif(ticket$labels, 4)

  image(z, col = col, axe = FALSE, main = main)
  box()
  axis(1, at = ticket$at, labels = ticket$labels)
}

create.colors <- function(x, alpha = 0.3, n = NULL){
  if(is.null(n)){
    tl.x <- length(x)
  } else{
    tl.x <- n
  }
  col <- rainbow(tl.x, start = 0, end = 2 / 6, alpha = alpha)
  col
}
create.cex <- function(x, cex = c(0.3, 1), n = NULL){
  if(is.null(n)){
    tl.x <- length(x)
  } else{
    tl.x <- n
  }
  cex <- seq(cex[1], cex[2], length = tl.x)
  cex
}
match.colcex <- function(x, colcex){
  tl.colcex <- length(colcex)
  tl.x <- length(x)
  if(tl.colcex < tl.x){
    tl.x.colcex <- tl.x %/% tl.colcex
    colcex <- c(rep(colcex, rep(tl.x.colcex, tl.colcex)),
                rep(colcex[tl.colcex], tl.x %% tl.colcex))
  }
  colcex <- colcex[order(x)]
  colcex
}


