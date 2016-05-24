summary.gazepath <-
function(object, ..., complete_only = FALSE){
  output <- numeric()
  for(i in 1:length(object[[16]])){
    sim <- object[[16]][[i]]
    l <- length(which(sim[,1] == 'f'))
    if(l != 0){
      if(complete_only == TRUE){
        if(length(which(sim[,1] == 's')) != 0){
          index <- sort(c(complete(sim, 'f'), complete(sim, 's')))
          if(length(index) != 0){
            output <- rbind(output, cbind(sim[index, c(1:4, 9:11)], 1:length(index), i))
          }
        }         
      } else {
        output <- rbind(output, cbind(sim[sim[,1] == 'f', c(1:4, 9:11)], 1:l, i))
      }
    }
  }
  if(length(output) == 0){
    print('There were no fixations or saccades classified, probably data quality of this particpant is very low')
    output <- data.frame(matrix(NA, 1, 7))
    names(output)[3:5] <- c('mean_x', 'mean_y', 'POGvar')
  } 
  names(output)[c(1:4, 8:9)] <- c('Value', 'Duration', 'Start', 'End', 'Order', 'Trial')
  row.names(output) <- 1:dim(output)[1]
  return(output)
}
