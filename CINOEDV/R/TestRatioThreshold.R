TestRatioThreshold <-
function(MaxOrder,RatioThreshold){
  # Check the parameter 'RatioThreshold'.
  # 
  # input
  #   MaxOrder: the maximum order of epistatic interactions that CINOEDV considered.
  #   RatioThreshold:
  #    The length of RatioThreshold is equal to the parameter 'MaxOrder'.
  #    Each element is a decimal in [0,1] with a character format.
  #
  #    For example, RatioThreshold <- c("1.0","0.5","0.2"), MaxOrder <- 3
  #
  # Junliang Shang
  # 3.27/2014
   
  if (length(RatioThreshold)==MaxOrder){
    for (i in 1:MaxOrder){
      CharNum <- nchar(RatioThreshold[i])
      DotNum <- 0
      for (j in 1:CharNum){
        status <- (substr(RatioThreshold[i],j,j) %in% 
                     c("0","1","2","3","4","5","6","7","8","9","."))
        if (status==FALSE){
          stop("   Error! ",RatioThreshold[i]," should be a decimal.\n")
        }
        if (substr(RatioThreshold[i],j,j)=="."){
          DotNum <- DotNum+1
        }
      }
      if (DotNum>1){
        stop("   Error! ",RatioThreshold[i]," should be a decimal.\n")
      }
      if ((as.numeric(RatioThreshold[i])>1) || 
            (as.numeric(RatioThreshold[i])<0)){
        stop("   Error! ",RatioThreshold[i]," should be within [0,1]\n")
      } 
    }
    cat("   The RatioThreshold :", as.numeric(RatioThreshold),"\n\n")    
  }else
  {
    stop(" Error! ",MaxOrder," elements should be in RatioThreshold.\n")
  }
}
