plot.gazepath <-
function(x, ..., i = 1){
  object <- x
  if(length(which(!is.na(object[[2]][[i]]))) < object[[10]] | length(which(!is.na(object[[3]][[i]]))) < object[[10]]){
    warning('There is not enough data to identify fixations and saccades in this trial')
  } else {
    layout(matrix(c(1, 1:3), 2, 2))
    plot(object[[2]][[i]], object[[3]][[i]], xlab = "X", ylab = 'Y', type = 'l', las = 1, ylim = c(max(object[[3]][[i]], na.rm = T), min(object[[3]][[i]], na.rm = T)))
    sim <- object[[16]][[i]]
    points(sim[sim[,1] == 'f', 9:10], pch = letters, cex = 3, col = 4)
    
    plot(object[[2]][[i]], ylim = c(-50, object[[14]][[i]]), type = 'l', xlab = 'Time (msec)', ylab = 'position', las = 1, xaxt = 'n')
    lines(object[[3]][[i]])
    axis(1, at = seq(0, length(object[[2]][[i]]), length.out = 6), labels = round(seq(0, length(object[[2]][[i]]) * (1000 / object[[10]]), length.out = 6)))
    fix <- sim[sim[,1] == 'f',3:4] / (1000 / object[[10]])
    rect(fix[,1], -50, fix[,2], 0, col = 3)
    
    if(object[[4]] != 'Tobii' & object[[4]] != 'Eyelink'){
      plot(object[[9]][[i]], typ = 'l', xlab = 'Time (msec)', ylab = 'Speed (deg/s)', las = 1, xaxt = 'n')
      axis(1, at = seq(0, length(object[[9]][[i]]), length.out = 6), labels = round(seq(0, length(object[[9]][[i]]) * (1000 / object[[10]]), length.out = 6)))
      segments(0, object[[7]][[i]], length(object[[9]][[i]]), object[[7]][[i]], col = 2, lwd= 2)
    }
  }
}
