  .ROSthetacalc <- function(ROS,FROS,BROS,THETA){
  	 c1        <- cos(THETA)
  	 s1        <- sin(THETA)
     c1        <- ifelse(c1==0,cos(THETA+.001),c1)
     ROStheta  <- (ROS - BROS)/(2*c1) + (ROS + BROS)/(2*c1)*
  		(FROS*c1*sqrt(FROS*FROS * c1*c1 + (ROS*BROS)*s1*s1) - (ROS*ROS - BROS*BROS)*s1*s1)/
  		(FROS*FROS*c1*c1 + ((ROS+BROS)/2)*((ROS+BROS)/2)*s1*s1)		                                                                       # /* 94 - 2009 */
  	 ROStheta}
