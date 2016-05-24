gconv<-function(g1,g2){
  #convolution of probability generating functions g1 and g2
  convolve(g1,rev(g2),type="open")
}
