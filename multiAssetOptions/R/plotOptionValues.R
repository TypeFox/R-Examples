plotOptionValues <- function(Y, fps = 10) {
  
  interval <- 1 / fps
  
  for (i in 1:dim(Y$value)[2]) {
    if (length(Y$dimS) == 1) {
      x <- Y$S[[1]]
      plot(Y$S[[1]], Y$value[,i], xlim = range(Y$S[[1]]), ylim = 
            range(Y$value[,i]), xlab = "Asset Price", ylab = "Option Value", 
            main = "Option Value", sub = paste("Time = ", round(Y$time[i], 
            digits = 4)), type="l", axes = TRUE, col = "black")
      Sys.sleep(interval)
    } else if (length(Y$dimS) == 2) {
      persp(Y$S[[1]], Y$S[[2]], array(Y$value[,i], dim=Y$dimS), 
            xlim = range(Y$S[[1]]), ylim = range(Y$S[[2]]), zlim = 
            range(Y$value), xlab = "Asset 1 Price", ylab = "Asset 2 Price", 
            zlab = "Option Value", main = "Dual Asset Option Value", sub = paste
            ("Time = ", round(Y$time[i], digits = 4)), theta = 45, phi = 25, 
            axes = TRUE, ticktype = "detailed",col = "lightblue", ltheta = 180, 
            shade = 0.75)
      Sys.sleep(interval)
    } else {
      stop("cannot plot higher than 3-D graphs")
    }
  }
  
}