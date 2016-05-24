
plot.tscgm.ar2 <- function(x, mat=c("precision","autoregression1", 
     "autoregression2"),...){
  mat=match.arg(mat)
  prec <- x$theta
  p=dim(prec)[2]
  if (mat == "precision") {
      prec <- prec
      nw_full <- network(prec)
      plot.network(nw_full,label = network.vertex.names(nw_full), 
         usearrows = FALSE,
         displayisolates = FALSE,...)
   }
 else if (mat == "autoregression1") {
      autoR1 <- x$gamma[1:p,] 
      nw_full <- network(autoR1, loops=TRUE)
      plot.network(nw_full,label = network.vertex.names(nw_full), 
        usearrows = TRUE,
        displayisolates = FALSE,...)
       }
else if (mat == "autoregression2") {
        d1= p + 1
        d2= 2* p
        autoR2 <- x$gamma[d1:d2,] 
        nw_full <- network(autoR2, loops=TRUE)
      plot.network(nw_full,label = network.vertex.names(nw_full), 
        usearrows = TRUE,
        displayisolates = FALSE, ...)
        }
}