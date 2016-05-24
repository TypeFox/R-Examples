# Reorders the latent variables to build a chain.
# The input is the adjacency matrix D from the inner model.
reorder <- function(D, n){
  if(!missing(n)) n <- min(n, ncol(D), nrow(D))
  else n <- min( ncol(D), nrow(D))
  chain <- NULL
  level <- NULL
  Dn <- D
  Dn[Dn!=0] <- 0
  for(i in n:1){
    tmp <- eval(parse(text=paste(rep("D", i), collapse=" %*% ")))
    if(all(tmp==0)) next
    chain <- unique(append(chain, sort(rownames(which(tmp!=0, TRUE)))))
    Dn <- Dn + tmp
  }
  if(length(chain) < ncol(D)) {
    indx <- which(! colnames(D) %in% chain)
    chain <- append(chain, sort(colnames(D)[indx]))
  }
  D <- D[chain, chain]
  source <- dimnames(D)[[1]][which(D!=0, TRUE)[,1]]
  target <- dimnames(D)[[2]][which(D!=0, TRUE)[,2]]
  sm <- cbind(source, target)
  attr(chain, "level") <- level
  # Dn is the sum of the 1, ..., n - step transition matrices
  return(list(chain=chain, Dn=Dn, n=n, strucmod=sm))
}
