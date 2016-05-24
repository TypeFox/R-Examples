posthocCheck <-
function(classification, x, y){
  class <- rle(classification)
  simple <- data.frame(class$values, class$lengths, 
                       c(1, cumsum(class$lengths) + 1)[-(length(class$values) + 1)], 
                       cumsum(class$lengths))
  lf <- length(which(class$values == 'f'))
  if(lf > 0){
    x_start <- y_start <- x_end <- y_end <- mean_x <- mean_y <- POGvar_px <- numeric()
    for(i in 1:dim(simple)[1]){
      x_start <- c(x_start, x[simple[i, 3]])
      y_start <- c(y_start, y[simple[i, 3]])
      x_end <- c(x_end, x[simple[i, 4]])
      y_end <- c(y_end, y[simple[i, 4]])
      mean_x <- c(mean_x, mean(x[simple[i,3] : simple[i,4]]))
      mean_y <- c(mean_y, mean(y[simple[i,3] : simple[i,4]]))
      POGvar_px <- c(POGvar_px, mean(dist(cbind(x[simple[i,3] : simple[i,4]], y[simple[i,3] : simple[i,4]]))))
    }
    simple <- cbind(simple, x_start, y_start, x_end, y_end, mean_x, mean_y, POGvar_px)
  }
  ## check for miss-classifications
  if(length(which(simple[,1] == 'f')) > 1){
    dis <- dist(simple[which(simple[,1] == 'f'),9:10])
    if(length(dis) == 1){
      dis.between <- dis
    } else {
      dis.between <- diag(as.matrix(dis)[-1,-dim(as.matrix(dis))])
    }
    dis.within <- simple[which(simple[,1] == 'f'),11]
    for(i in 1:length(dis.between)){
      if(dis.between[i] < dis.within[i] + dis.within[i + 1]){
        index <- c(simple[which(simple[,1] == 'f')[i],3], simple[which(simple[,1] == 'f')[i + 1],4])
        classification[index[1] : index[2]] <- 'f'
        x[index[1] : index[2]][is.na(x[index[1]: index[2]])] <- mean(x[index[1] : index[2]], na.rm = T)
        y[index[1] : index[2]][is.na(y[index[1]: index[2]])] <- mean(y[index[1] : index[2]], na.rm = T)
      }
    }
  }
  return(list(classification, x, y))
}
