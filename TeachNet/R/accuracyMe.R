accuracy.me <- function(obs, predict, thres=0.5) {# calculate accurancy of modell

  # confusion matrix
  conf <- confusion(predict,obs=obs, threshold=thres)
  
  leng <- length(obs)
  if(leng==1){leng<-nrow(obs)}
  
  zero <- (conf[2,1] + conf[3,1])/leng
  network <- (conf[2,1]+conf[3,2])/leng
  one <- 1-zero

      return(cbind(zero,one,network))
      
}
