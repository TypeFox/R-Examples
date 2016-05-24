rplau<-function(x){
  # Goal: Calculate a plausibility ratio of a proposition A,
  # namely Pl(A)/Pl(¬A)
  # 
  # x is the table resulting from the function belplau.
  # 
  # Example of an input table x
  #       Belief     Plausibility
  # [1,] 0.49743668  0.502461296
  # [2,] 0.00010202  0.005126633
  # [3,] 0.49743668  0.502461296
  # [4,] 1.00000000  1.000000000
  zr<-x
  zplau<-zr[,2]
  zbel<-zr[,1]
  zrplau<-zplau/(1-zbel)
  return(zrplau)
}