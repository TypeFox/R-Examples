
plot.tscgm <- function(x, mat=c("precision","autoregression"),...){
  mat=match.arg(mat)
  if (mat == "precision") {
      prec <- x$theta
      #colnames(prec) <- rownames(prec) <- colnames(data)
      nw_full <- network(prec)
      plot.network(nw_full,label = network.vertex.names(nw_full), usearrows = FALSE,
         displayisolates = FALSE,...)
   }
 else if (mat == "autoregression") {
      autoR <- x$gamma
      #colnames(autoR) <- rownames(autoR) <- colnames(data)
      nw_full <- network(autoR, loops=TRUE)
      plot.network(nw_full,label = network.vertex.names(nw_full), usearrows = TRUE,
        displayisolates = FALSE,...)
   }
}