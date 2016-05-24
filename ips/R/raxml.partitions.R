## PACKAGE: ips
## CALLED BY: USER
## AUTHOR: Christoph Heibl (at gmx.net)
## LAST UPDATE: 2014-08-07

raxml.partitions <- function(...){
  
  x <- list(...)
  p <- cbind(rep(1, length(x)), sapply(x, ncol))
  for ( i in 2:nrow(p) ){
    p[i, 1] <- p[i - 1, 2] + 1
    p[i, 2] <- p[i, 1] + p[i, 2] -1
  }
  p <- data.frame(type = "DNA", 
             locus = rownames(p),
             begin = p[, 1],
             end = p[, 2])
  rownames(p) <- NULL
  p
}