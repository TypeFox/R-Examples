Ymean_compare <-
function(Y1_mean,Y2_mean,matrixY1,matrixY2,matrixZ_chain,matrixW1,matrixW2,matrixW_chain,matrixX,matrixX_chain,result_file_name,plot_name){


  MEANY1<-apply(matrixY1,1,mean)
  MEANY2<-apply(matrixY2,1,mean)

  SDY1<-apply(matrixY1,1,sd)
  SDY2<-apply(matrixY2,1,sd)

  scaled_Y1_mean <- (Y1_mean-MEANY1)/SDY1
  scaled_Y2_mean <- (Y2_mean-MEANY2)/SDY2

  total_iter<-dim(matrixZ_chain[[1]])[1]
  num_Y1matrix<-length(scaled_Y1_mean)
  num_Y2matrix<-length(scaled_Y2_mean)
  est_Y1_mean_chain <- matrix(0,total_iter,num_Y1matrix)
  est_Y2_mean_chain <- matrix(0,total_iter,num_Y2matrix)
  RMSE_Y_MEAN<-matrix(0,total_iter,2)



  for(i in 1:total_iter){

    EST_W1<-matrix(matrixW_chain[[1]][i,],nrow=dim(matrixW1)[1],ncol=dim(matrixW1)[2])
    EST_W2<-matrix(matrixW_chain[[2]][i,],nrow=dim(matrixW2)[1],ncol=dim(matrixW2)[2])
    EST_X<-matrix(matrixX_chain[i,],nrow=dim(matrixX)[1],ncol=dim(matrixX)[2])

    est_Y1_mean <- EST_W1 %*% EST_X
    est_Y2_mean <- EST_W2 %*% EST_X

    est_Y1_mean_chain[i,] <- est_Y1_mean
    est_Y2_mean_chain[i,] <- est_Y2_mean


    RMSE_Y_MEAN[i,1] <- sqrt(sum((est_Y1_mean-scaled_Y1_mean)^2)/length(est_Y1_mean))
    RMSE_Y_MEAN[i,2] <- sqrt(sum((est_Y2_mean-scaled_Y2_mean)^2)/length(est_Y2_mean))
  }



  save(scaled_Y1_mean,scaled_Y2_mean,est_Y1_mean_chain,est_Y2_mean_chain,RMSE_Y_MEAN,file=result_file_name)

  pdf(plot_name,20,10)
  par(mfrow=c(1,2))
  plot(RMSE_Y_MEAN[,1])
  plot(RMSE_Y_MEAN[,2])
  dev.off()

}

