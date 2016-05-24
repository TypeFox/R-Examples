threshs <-
function(my.plant, min.thre=200, panel, ci=1.9){
  #my.plant <- new.whole.data[[1]]
  sta <- big.peaks.col(my.plant$yy, min.thre)
  sta$wei <- my.plant$xx[sta$pos]
  # subset only the ones participating in the panel
  hoho <- which(sta$wei < panel[2] & sta$wei > panel[1])
  if(length(hoho) > 0){
    sta <- list(pos=sta$pos[hoho], hei=sta$hei[hoho], wei= sta$wei[hoho])
    #plot(x=my.plant$xx, y=my.plant$yy, type="l", xlim=c(panel[1],panel[2]+50))
    #abline(h=mm);  abline(h=mm-se2);  abline(h=mm+se2)
    z <- which(sta$hei == max(sta$hei))
    mm <- mean(sta$hei)
    if(length(sta$hei) > 2){
      se2 <- ci * (sd(sta$hei[-z])/sqrt(length(sta$hei[-z])))
    }else{ se2 <- ci * (sd(sta$hei)/sqrt(length(sta$hei)))}
    
    
    newthre <- mm-se2
  }else{newthre <- 500}
  
  return(newthre)
}
