owa<-function (x,y) 
{
  C5 = numeric()
  diff = y-x
  for(i in 1:length(x)){
    for(j in i:length(y)){
      C5 = c(C5, (diff[i]+diff[j])/2)
    }
  }
  owa=sort(C5)
  h.l = median(owa)
  list(owa=owa, h.l =h.l)
  
}
