bag.aucoob <- function(bag_pltr, xdata, Y.name){
  Index = sort(unique(unlist(bag_pltr$IND_OOB)))
  YOOB = xdata[, Y.name][Index]
  Ypred = sapply(1:length(bag_pltr$CUT), function(vv) YOOB*(bag_pltr$LOST[,vv]==0) +(1 - YOOB)*(bag_pltr$LOST[,vv]== 1))
  TPROOB <- c()
  FPROOB <- c()
  j <- 0
  for(cut in bag_pltr$CUT)
  {
    j <- j+1
    predj <- Ypred[,j]
    conf1 <- table(predj, YOOB)
    if(sum(predj) == length(YOOB)){
      TPROOB <- c(TPROOB, 1)
      FPROOB <- c(FPROOB, 1)
    } else if (sum(predj) == 0){
      TPROOB <- c(TPROOB, 0)
      FPROOB <- c(FPROOB, 0)
    } else{
      TPROOB <- c(TPROOB, conf1[4]/(conf1[4] + conf1[3]))
      FPROOB <- c(FPROOB, conf1[2]/(conf1[2] + conf1[1]))
    }
  }
  TPR = sort(c(0,1,TPROOB))
  FPR <- sort(c(0,1,FPROOB))
  AUC <- sum(diff(FPR) * (TPR[-1] + TPR[-length(TPR)]) / 2)
  OOB <- bag_pltr$EOOB
  return(list(AUCOOB = AUC, TPR = TPR, FPR = FPR, OOB = OOB))
}
