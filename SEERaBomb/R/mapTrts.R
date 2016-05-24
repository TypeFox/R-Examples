mapTrts<-function(D){
  # adds a trt column to define cancer treatments of interest as a factor
  # this is how radiatn breaksdown. 0=none, 1=ExtBeam, 2+3=iso+implants, 4=extB+other
  # 5=NOS, 6=other rad (other=other than beam), 7=refused, 8+9=Unknown
  D$trt="noRad"  # will be left as 0 and 7. Do it this way to initialize the vector
  D$trt[D$radiatn%in%c(8,9)] ="unk"
  D$trt[D$radiatn>0 & D$radiatn<7] ="rad"
  D$trt=factor(D$trt,levels=c("rad","noRad","unk"))
  D
}  
