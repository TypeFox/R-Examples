simplify <-
function(classification, x, y, Hz, D, width_px, width_mm){
  class <- rle(classification)
  simple <- data.frame(class$values, class$lengths, 
                       c(1, cumsum(class$lengths) + 1)[-(length(class$values) + 1)], 
                       cumsum(class$lengths))
    if(length(which(class$values == 'f')) > 0){
      x_start <- y_start <- x_end <- y_end <- mean_x <- mean_y <- POGvar <- numeric()
      for(i in 1:dim(simple)[1]){
        x_start <- c(x_start, x[simple[i, 3]])
        y_start <- c(y_start, y[simple[i, 3]])
        x_end <- c(x_end, x[simple[i, 4]])
        y_end <- c(y_end, y[simple[i, 4]])
        mean_x <- c(mean_x, mean(x[simple[i,3] : simple[i,4]]))
        mean_y <- c(mean_y, mean(y[simple[i,3] : simple[i,4]]))
        POGvar <- c(POGvar, mean(as.matrix(dist(cbind(c(mean_x[length(mean_x)], x[simple[i,3] : simple[i,4]]), c(mean_y[length(mean_y)], y[simple[i,3] : simple[i,4]]))))[-1,1]))
      }
      ## Calculate saccade amplitude and transform POGvar from pixels to degrees of visual angle
      ss <- which(class$values == 's')
      POGvar[ss] <- sqrt((x_start[ss] - x_end[ss]) ^ 2 + (y_start[ss] - y_end[ss]) ^ 2)
      POGvarSacAmp <- atan((POGvar / 2) / mean(D, na.rm = T)) * (180 / pi) * (width_mm / width_px) * 2
      
      simple <- data.frame(class$values, class$lengths * (1000/Hz), 
                           c(1, cumsum(class$lengths * (1000/Hz)) + 1)[-(length(class$values) + 1)], 
                           cumsum(class$lengths * (1000/Hz)), x_start, y_start, x_end, y_end, mean_x, mean_y, POGvarSacAmp)
      names(simple)[1:4] <- c('Value', 'Dur', 'Start', 'End')
    }
  # Remove NA values
  if(length(which(is.na(simple[,1]))) != 0){
    simple <- simple[-which(is.na(simple[,1])),]
  }
  return(simple)
}
