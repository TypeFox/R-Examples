find.limits <-
function(y, u, A, b){
  z = y - sum(u*y)*u
  
  Au = A%*%u
  Az = A%*%z
  
  which.pos = Au > 0
  which.neg = Au < 0
  
  vals = (b - Az)/Au
  
  options(warn=-1)
  v.p = min (vals[which.pos])
  v.m = max (vals[which.neg])
  options(warn=0)
  
  c(v.m, v.p)
}
