CombinationEntropy <-
function(data){
  # compute combination entropy of given variables. There variables can be 
  # one SNP, the class label, SNP-combination, or SNP-combination-class label.
  #
  # input
  #    data: the given variables. for example, data <- cbind(pts[,factor],t(class))
  #
  # output
  #    Combination_Entropy_Value: combination entropy value
  #
  # Junliang Shang
  # 3.10/2014

  data <- t(data)
  data <- data-min(data)
  RowColNum <- dim(data)
  ZZ <- 3^c(RowColNum[1]:1-1)
  tempIndex <- ZZ%*%data+1
  HistNum <- hist(tempIndex,breaks=
                    seq(min(tempIndex),max(tempIndex),
                        (max(tempIndex)-min(tempIndex))/3^RowColNum[1]),
                  plot=FALSE)
  Result <- HistNum$counts
  PResult=Result/sum(Result)
  EntropyFactor <- -PResult%*%log2(PResult+2.2204e-016)
  rm(RowColNum,ZZ,tempIndex,HistNum,Result,PResult,data)
  
  list(Combination_Entropy_Value=EntropyFactor)
}
