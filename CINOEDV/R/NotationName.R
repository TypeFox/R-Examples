NotationName <-
function(MaxOrder,SingleEffect,TwoEffect,ThreeEffect,
                         FourEffect,FiveEffect,SNPNames=NA){
  # Notation of real SNP Name
  #
  # input
  #    MaxOrder: the specified maximum order, must be setted as 1,2,3,4 or 5.
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
  #    SNPName: Row -> 1, Column -> SNP Name
  # output
  #    SingleEffect: The notation of SingleEffect
  #    TwoEffect: The notation of TwoEffect
  #    ThreeEffect: The notation of ThreeEffect
  #    FourEffect: The notation of FourEffect
  #    FiveEffect: The notation of FiveEffect
  #
  # Junliang Shang
  # 4.20/2014
  
  if(length(SNPNames)>1){
    if(MaxOrder>=1){
      for (i in 1:nrow(SingleEffect)){
        SingleEffect[i,1] <- SNPNames[as.numeric(SingleEffect[i,1])]
      }
    }
    
    if(MaxOrder>=2){
      for (i in 1:nrow(TwoEffect)){
        TwoEffect[i,1] <- SNPNames[as.numeric(TwoEffect[i,1])]
        TwoEffect[i,2] <- SNPNames[as.numeric(TwoEffect[i,2])]
      }
    }
    
    if(MaxOrder>=3){
      for (i in 1:nrow(ThreeEffect)){
        ThreeEffect[i,1] <- SNPNames[as.numeric(ThreeEffect[i,1])]
        ThreeEffect[i,2] <- SNPNames[as.numeric(ThreeEffect[i,2])]
        ThreeEffect[i,3] <- SNPNames[as.numeric(ThreeEffect[i,3])]
      }
    }
    
    if(MaxOrder>=4){
      for (i in 1:nrow(FourEffect)){
        FourEffect[i,1] <- SNPNames[as.numeric(FourEffect[i,1])]
        FourEffect[i,2] <- SNPNames[as.numeric(FourEffect[i,2])]
        FourEffect[i,3] <- SNPNames[as.numeric(FourEffect[i,3])]
        FourEffect[i,4] <- SNPNames[as.numeric(FourEffect[i,4])]
      }
    }
    
    if(MaxOrder>=5){
      for (i in 1:nrow(FiveEffect)){
        FiveEffect[i,1] <- SNPNames[as.numeric(FiveEffect[i,1])]
        FiveEffect[i,2] <- SNPNames[as.numeric(FiveEffect[i,2])]
        FiveEffect[i,3] <- SNPNames[as.numeric(FiveEffect[i,3])]
        FiveEffect[i,4] <- SNPNames[as.numeric(FiveEffect[i,4])]
        FiveEffect[i,5] <- SNPNames[as.numeric(FiveEffect[i,5])]
      }
    }
  }
  
  list(SingleEffect=SingleEffect,TwoEffect=TwoEffect,ThreeEffect=ThreeEffect,
       FourEffect=FourEffect,FiveEffect=FiveEffect)
}
