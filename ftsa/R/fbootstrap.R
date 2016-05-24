fbootstrap = function(data, estad = func.mean, alpha = 0.05, nb = 200, suav = 0.0, media.dist = FALSE, 
                      graph = FALSE, ...)
{
  plotx <- data$x
  data <- t(data$y)
  nrow <- dim(data)[1]
  ncol <- dim(data)[2]
  estmues <- estad(data, ...)
  if(is.null(nrow) || is.null(ncol)) 
     stop("One of the data is not a matrix")
  distboot <- matrix(NA, nrow = nb)
  estboot <- matrix(NA, nrow = nb, ncol = ncol)
  for(i in 1:nb){
      bdata <- data[sample(1:nrow, size = nrow, replace = T),]
      if(suav > 0){
         bdata <- bdata + mvrnorm(n = nrow, rep(0, ncol), var(data) * suav)
      }
      estboot[i,] <- estad(bdata, ...)
  }
  if(media.dist == TRUE){
     centro = func.mean(estboot)
  } 
  else{
       centro <- estmues
  }
  for(i in 1:nb){
      distboot[i] <- metri.p(centro,estboot[i,])
  }
  dist <- max(distboot[rank(distboot) <= floor((1 - alpha) * nb)])
  if(graph){
     dev.new()
     plot(fts(seq(plotx[1], plotx[length(plotx)], length = ncol), t(estboot)),type = "n")
         if(distboot[i] <= dist){
            lines(fts(seq(plotx[1], plotx[length(plotx)], length = ncol), t(estboot)), lty = 2, col = 3, 
                  pch = c(1:9, 0, letters, LETTERS))
         }
         else{
            lines(fts(seq(plotx[1], plotx[length(plotx)], length = ncol), t(estboot)), lty = 2, col = 1, 
                  pch = c(1:9, 0, letters, LETTERS))
        }
     lines(seq(plotx[1], plotx[length(plotx)], length = ncol), estmues, lwd = 3, lty = 1, col = 2)
     legend("bottomleft", legend = c("Data", "Estimate bootstrap data", "Estimate function"), 
            lty = c(1,2,2), lwd = c(3,1,1), col = c(2,3,1), cex = .8)
  }
  return(list("estimate" = estmues, "max.dist" = dist, "rep.dist" = distboot, "resamples" = estboot,
         "center" = centro))    
}
