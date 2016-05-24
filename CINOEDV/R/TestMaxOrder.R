TestMaxOrder <-
function(MaxOrder){
  # Check the parameter 'MaxOrder'.
  # 'MaxOrder' is the maximum order of epistatic interactions that CINOEDV
  # considered. In this version, it must be setted as 2,3,4,or 5.
  #
  # input
  #     MaxOrder <-  "2", "3", "4", or "5".
  #
  # Junliang Shang
  # 3.25/2014
  
  if (MaxOrder %in% c("2","3","4","5")){
    cat("    The Maximum order :", MaxOrder,"\n\n")
  }else
  {
    stop("    The Maximum order should be setted as 2,3,4 or 5 !\n")
  }
}
