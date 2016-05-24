CoInformation <-
function(pts,class,factor){
  # compute co-information value of a give SNP combination.
  # reference: http://en.wikipedia.org/wiki/Interaction_information
  #
  # input
  #    pts: Row -> Sample, Column -> SNP
  #         1 -> AA
  #         2 -> Aa
  #         3 -> aa
  #    class: Row -> 1, Column -> class label
  #           1 -> case
  #           2 -> control
  #    factor: the considered SNP or SNP-combination, for example, factor <- 5
  #            or factor <- c(2,5)
  #
  # output
  #    Co_Information_Value: co-information value of a give SNP combination.
  #
  # Junliang Shang
  # 3.12/2014
 
  data <- cbind(pts[,factor],t(class))
  FeatureNum <- ncol(data)
  value <- 0
  for (i in 1:FeatureNum){
    A <- combn(1:FeatureNum,i)
    for (j in 1:ncol(A)){
      value <- value-(-1)^(FeatureNum-i)*
        CombinationEntropy(data[,A[,j]])$Combination_Entropy_Value
    }
  }
  
  list(Co_Information_Value=value)
}
