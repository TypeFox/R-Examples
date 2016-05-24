Tobii <-
function(data, distance, offset = .9, thres_dur = 100, Hz){
  clas <- 'f'
  for(i in 1:dim(data)[1]){
    if(is.na(data[i,1]) | is.na(data[i,2]) | is.na(distance[i])){
      clas[i] <- 'u'
    } else {
      if(i == 1){
        clas <- 'f'
        index <- 1
        pos.fix <- data[1,]
      } else {
        pos.fix <- apply(data[index : i,], 2, function(x) mean(x, na.rm = T))
      }
      check <- dist(rbind(pos.fix, data[i,]))
      margin <- tan((offset / 2) * pi/180) * distance[i] * 2 * (1024 / 272)
      clas[i] <- ifelse(check < margin, 'f', 'u')
    }
    if(clas[i] == 'u'){
      index <- i
    }
  }
  thres_dur <- thres_dur * Hz / 1000
  
  rle <- rle(clas)
  classify <- numeric()
  for(i in 1:length(rle$values)){
    if(is.na(rle$values[i])){
      classify <- c(classify, rep(NA, rle$lengths[i]))
    } else{
      if(rle$values[i] == 'f' & rle$lengths[i] >= thres_dur){
        classify <- c(classify, rep('f', rle$lengths[i]))
      }
      if(rle$values[i] == 'f' & rle$lengths[i] < thres_dur){
        classify <- c(classify, rep('u', rle$lengths[i]))
      }
      if(rle$values[i] == 'u'){
        classify <- c(classify, rep('s', rle$lengths[i]))
      }
    }
  }
  return(classify)
}
