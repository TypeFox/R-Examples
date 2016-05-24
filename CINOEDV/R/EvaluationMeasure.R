EvaluationMeasure <-
function(factor,pts,class,measure=1){
  # This function provides several evaluation measures for measuring main or
  # interaction effect of a SNP or SNP-combination.
  #
  # input
  #    factor: the considered SNP or SNP-combination, for example, factor <- 5
  #            or factor <- c(2,5)
  #    pts: Row -> Sample, Column -> SNP
  #         1 -> AA
  #         2 -> Aa
  #         3 -> aa
  #    class: Row -> 1, Column -> class label
  #           1 -> case
  #           2 -> control
  #    measure: the label of current used evaluation measure
  #            1 -> The classic co-information measure
  #            2 -> The Normalized co-information measure
  #            3 -> TingHu's Co-Information Measure
  #
  # output
  #    Value: the main or interaction effect value corresponding to a single SNP
  #    or a SNP combination.
  #
  # Junliang Shang
  # 3.26/2014
  
  # if the label of measure input error, then the first one is selected instead.
  if (!(measure %in% c(1,2,3))){
    measure <- 1
  }
  
  if (measure==1){
    #############################
    # Classic Co-Information Measure
    #############################
    
    Effect <- CoInformation(pts,class,factor)
    Value <- Effect$Co_Information_Value
  }
  
  if (measure==2){
    #############################
    # The Normalized Co-Information Measure
    #############################
    # It is normalized by dividing the entropy of the case-control status H(C), 
    # and thus give the percentage of predicting information on the phenotype status
    
    Effect <- CoInformation(pts,class,factor)
    Value <- Effect$Co_Information_Value
    
    DC <- t(class)
    Value <- Value/CombinationEntropy(DC)$Combination_Entropy_Value
  }
  
  if (measure==3){
    #############################
    # TingHu's Co-Information Measure
    #############################
    # reference: an information-gain approach to detecting three-way epistatic 
    # interactions in genetic association studies.
    
    if (length(factor)<3){
      Effect <- CoInformation(pts,class,factor)
      Value <- Effect$Co_Information_Value
      
      DC <- t(class)
      Value <- Value/CombinationEntropy(DC)$Combination_Entropy_Value
    }else
    {
      D123C <- cbind(pts[,factor],t(class))
      D123 <- pts[,factor]
      D1C <- cbind(pts[,factor[1]],t(class))
      D2C <- cbind(pts[,factor[2]],t(class))
      D3C <- cbind(pts[,factor[3]],t(class))
      D1 <- pts[,factor[1]]
      D2 <- pts[,factor[2]]
      D3 <- pts[,factor[3]]
      DC <- t(class)
      
      Value <- (-CombinationEntropy(D123C)$Combination_Entropy_Value+
                  CombinationEntropy(D123)$Combination_Entropy_Value+
                  CombinationEntropy(D1C)$Combination_Entropy_Value+
                  CombinationEntropy(D2C)$Combination_Entropy_Value+
                  CombinationEntropy(D3C)$Combination_Entropy_Value-
                  CombinationEntropy(D1)$Combination_Entropy_Value-
                  CombinationEntropy(D2)$Combination_Entropy_Value-
                  CombinationEntropy(D3)$Combination_Entropy_Value-
                  CombinationEntropy(DC)$Combination_Entropy_Value-
                  CombinationEntropy(DC)$Combination_Entropy_Value-
                  max(0,CoInformation(pts,class,c(factor[1],factor[2]))$Co_Information_Value)-
                  max(0,CoInformation(pts,class,c(factor[1],factor[3]))$Co_Information_Value)-
                  max(0,CoInformation(pts,class,c(factor[2],factor[3]))$Co_Information_Value))
      
      Value <- Value/CombinationEntropy(DC)$Combination_Entropy_Value
    }
  }
  
  list(Value=Value)
}
