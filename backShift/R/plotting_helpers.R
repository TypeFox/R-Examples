autocurve.edges2 <-function (graph, start = 0.5)
{
  cm <- count.multiple(graph)
  mut <-is.mutual(graph)  #are connections mutual?
  el <- apply(get.edgelist(graph, names = FALSE), 1, paste,
              collapse = ":")
  ord <- order(el)
  res <- numeric(length(ord))
  p <- 1
  while (p <= length(res)) {
    m <- cm[ord[p]]
    mut.obs <-mut[ord[p]] #are the connections mutual for this point?
    idx <- p:(p + m - 1)
    if (m == 1 & mut.obs==FALSE) { #no mutual conn = no curve
      r <- 0
    }
    else {
      r <- seq(-start, start, length = m)
    }
    res[ord[idx]] <- r
    p <- p + m
  }
  res
}



## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

## Add an alpha value to a colour
vary.alpha <- function(alpha, col){  
  col <- col2rgb(col)/255
  sapply(alpha, function(x) rgb(col[1], col[2], col[3], alpha=x))  
}

convert <- function(valvec, range.min = 0.2, range.max = 0.9){
  min.val <- min(valvec)
  max.val <- max(valvec)
  sapply(valvec, function(i) (((i - min.val)*(range.max - range.min))/(max.val - min.val)) + range.min )
}

convert.given.min <- function(valvec, min.val, max.val, range.min = 0.2, range.max = 0.9){
  sapply(valvec, function(i) (((i - min.val)*(range.max - range.min))/(max.val - min.val)) + range.min )
}