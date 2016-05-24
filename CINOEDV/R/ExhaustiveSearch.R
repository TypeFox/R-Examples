ExhaustiveSearch <-
function (pts,class,MaxOrder=3,measure=1,alpha=0){
  # Exhaustively compute main effects of all SNPs and interaction effects of 
  # all SNP-combinations within the maximum order.
  #
  # input
  #    pts: Row -> Sample, Column -> SNP
  #         1 -> AA
  #         2 -> Aa
  #         3 -> aa
  #    class: Row -> 1, Column -> class label
  #         1 -> case
  #         2 -> control
  #    MaxOrder: the specified maximum order, must be setted as 1,2,3,4 or 5.
  #             By default, 3.
  #    measure: the label of current used evaluation measure
  #            1 -> The classic co-information measure
  #            2 -> The Normalized co-information measure
  #            3 -> TingHu's Co-Information Measure
  #    alpha: the lower threshold of effects, either main effects or interaction
  #           effects, must be higher or equal to 0, and by default, 0.
  #
  # output
  #    SingleEffect: There are 2 columns. The first column saves all SNPs, and the 
  #                  second column saves their corresponding effects. Descending 
  #                  save according to their effects.                   
  #    TwoEffect: There are 3 columns. The first 2 columns save all 2-SNP combinations,
  #               and the last column saves their corresponding effects. Descending 
  #               save according to their interaction effects.
  #    ThreeEffect: There are 4 columns. The first 3 columns save all 3-SNP combinations,
  #                 and the last column saves their corresponding effects. Descending 
  #                 save according to their interaction effects.
  #    FourEffect: There are 5 columns. The first 4 columns save all 4-SNP combinations,
  #                and the last column saves their corresponding effects. Descending 
  #                save according to their interaction effects.
  #    FiveEffect: There are 6 columns. The first 5 columns save all 5-SNP combinations,
  #                and the last column saves their corresponding effects. Descending
  #                save according to their interaction effects.
  #
  # Junliang Shang
  # 3.26/2014
  
  SNPNum <- ncol(pts)
  if (alpha<0) alpha <- 0
  
  SingleEffect <- array(0,dim=c(1,1,1))
  TwoEffect <- array(0,dim=c(1,1))
  ThreeEffect <- array(0,dim=c(1,1))
  FourEffect <- array(0,dim=c(1,1))
  FiveEffect <- array(0,dim=c(1,1))
  
  cat("    Exhaustive Search !\n\n")
  
  #############################
  # 1-order effect, that is, main effect
  #############################
  if (MaxOrder>=1){
    
    tic()   
    cat("   1-order effect computing ... \n")
    
    A <- combn(1:SNPNum,1)
    Effect <- sapply(A,EvaluationMeasure,pts,class,measure)
    Value <- sapply(A,function(x,Effect) Effect[[x]],Effect)
    SingleEffect <- t(rbind(A,Value))
    SingleEffect[,2] <- abs(SingleEffect[,2])
    dim(SingleEffect) <- c(length(SingleEffect)/2,2)
    colnames(SingleEffect) <- c("SNP","Value")
    
    SingleEffect <- SingleEffect[order(SingleEffect[,2],decreasing=T),]
    
    toc()
  }
  
  #############################
  # 2-order effect, that is, interaction effect of two SNPs
  #############################
  if (MaxOrder>=2){
    tic()
    cat("   2-order effect computing ... \n")
    
    A <- combn(1:SNPNum,2)
    Effect <- apply(X=A,MARGIN=2,EvaluationMeasure,pts,class,measure)
    Value <- sapply(1:ncol(A),function(x,Effect) Effect[[x]]$Value,Effect)
    TwoEffect <- t(rbind(A,Value))
    TwoEffect <- TwoEffect[TwoEffect[,3]>alpha]
    dim(TwoEffect) <- c(length(TwoEffect)/3,3)
    colnames(TwoEffect) <- c("SNP1","SNP2","Value")
    
    TwoEffect <- TwoEffect[order(TwoEffect[,3],decreasing=T),]
    
    toc()
  }
  
  #############################
  # 3-order effect, that is, interaction effect of three SNPs
  #############################
  if (MaxOrder>=3){
    tic()
    cat("   3-order effect computing ... \n")
    
    A <- combn(1:SNPNum,3)
    Effect <- apply(X=A,MARGIN=2,EvaluationMeasure,pts,class,measure)
    Value <- sapply(1:ncol(A),function(x,Effect) Effect[[x]]$Value,Effect)
    ThreeEffect <- t(rbind(A,Value))
    ThreeEffect <- ThreeEffect[ThreeEffect[,4]>alpha]
    dim(ThreeEffect) <- c(length(ThreeEffect)/4,4)
    colnames(ThreeEffect) <- c("SNP1","SNP2","SNP3","Value")
    
    ThreeEffect <- ThreeEffect[order(ThreeEffect[,4],decreasing=T),]
    
    toc()
  }
  
  #############################
  # 4-order effect, that is, interaction effect of four SNPs
  #############################
  if (MaxOrder>=4){
    tic()
    cat("   4-order effect computing ... \n")
    
    A <- combn(1:SNPNum,4)
    Effect <- apply(X=A,MARGIN=2,EvaluationMeasure,pts,class,measure)
    Value <- sapply(1:ncol(A),function(x,Effect) Effect[[x]]$Value,Effect)
    FourEffect <- t(rbind(A,Value))
    FourEffect <- FourEffect[FourEffect[,5]>alpha]
    dim(FourEffect) <- c(length(FourEffect)/5,5)
    colnames(FourEffect) <- c("SNP1","SNP2","SNP3","SNP4","Value")
    
    FourEffect <- FourEffect[order(FourEffect[,5],decreasing=T),]
    
    toc()
  }
  
  #############################
  # 5-order effect, that is, interaction effect of five SNPs
  #############################
  if (MaxOrder>=5){
    tic()
    cat("   5-order effect computing ... \n")
    
    A <- combn(1:SNPNum,5)
    Effect <- apply(X=A,MARGIN=2,EvaluationMeasure,pts,class,measure)
    Value <- sapply(1:ncol(A),function(x,Effect) Effect[[x]]$Value,Effect)
    FiveEffect <- t(rbind(A,Value))
    FiveEffect <- FiveEffect[FiveEffect[,6]>alpha]
    dim(FiveEffect) <- c(length(FiveEffect)/6,6)
    colnames(FiveEffect) <- c("SNP1","SNP2","SNP3","SNP4","SNP5","Value")
    
    FiveEffect <- FiveEffect[order(FiveEffect[,6],decreasing=T),]
    
    toc()
  }
  
  #############################
  # Other orders, return Error
  #############################
  if (MaxOrder>5)
  {
    stop("Error! 'MaxOrder' has exceeded the maximum value !\n")
  }
  
  #############################
  # Return results
  #############################
  list(SingleEffect=SingleEffect,TwoEffect=TwoEffect,ThreeEffect=ThreeEffect,
       FourEffect=FourEffect,FiveEffect=FiveEffect)
}
