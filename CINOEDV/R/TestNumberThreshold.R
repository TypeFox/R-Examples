TestNumberThreshold <-
function(MaxOrder,NumberThreshold){
  # Check the parameter 'NumberThreshold'.
  # 
  # input
  #   MaxOrder: the maximum order of epistatic interactions that CINOEDV considered.
  #   NumberThreshold:
  #    The length of NumberThreshold is equal to the parameter 'MaxOrder'.
  #    Each element is a integer with a character format.
  #
  #    For example, NumberThreshold <- c("5","20","10"), MaxOrder <- 3
  #
  # Junliang Shang
  # 3.28/2014
  
  
  if (length(NumberThreshold)==MaxOrder){
    for (i in 1:MaxOrder){
      CharNum <- nchar(NumberThreshold[i])
      DotNum <- 0
      for (j in 1:CharNum){
        status <- (substr(NumberThreshold[i],j,j) %in% 
                     c("0","1","2","3","4","5","6","7","8","9","."))
        if (status==FALSE){
          stop("   Error! ",NumberThreshold[i]," should be a number.\n")
        }
        if (substr(NumberThreshold[i],j,j)=="."){
          DotNum <- DotNum+1
        }
      }
      if (DotNum>1){
        stop("   Error! ",NumberThreshold[i]," should be a number.\n")
      }
      if (floor(as.numeric(NumberThreshold[i]))
          !=as.numeric(NumberThreshold[i])){
        stop("   Error! ",NumberThreshold[i]," should be an integer.\n")
      } 
    }
    cat("   The NumberThreshold :", as.numeric(NumberThreshold),"\n\n")
  }else
  {
    stop(" Error! ",MaxOrder," elements should be in NumberThreshold.\n")
  }
}
