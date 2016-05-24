colorhex <-
function(n,beg=1,end=256)
{ 
  if (beg < 1 || end < 1 || beg > 256 || end > 256) 
    stop("`beg' and `end' must be numbers in the interval [1,256]")
  c(gray(0.6),rev(rgb((1:n)/n, 0,0)))
}