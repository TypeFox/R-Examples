### function to implement the soft-, hard, SCAD thresholding rule
# Input variables:
#	z: argument
#	type: thresholding rule
#		1 = (Adaptive) LASSO (default)
#		2 = hard thresholding
#		3 = SCAD
# 	delta: thresholding level
# 	a: default choice for SCAD penalty


thresh <- function(z,delta,  type=1, a=3.7){
  if(type==1){
    return(sign(z)*(abs(z)>=delta)*(abs(z)-delta))
  }
  
  if(type==2){
    return(z*(abs(z)>delta))
  }
  
  if(type==3){
    return(sign(z)*(abs(z)>=delta)*(abs(z)-delta)*(abs(z)<=2*delta)+
	((a-1)*z-sign(z)*a*delta)/(a-2)*(2*delta<abs(z))*(abs(z)<=a*delta)+z*(abs(z)>a*delta))
  }
}
