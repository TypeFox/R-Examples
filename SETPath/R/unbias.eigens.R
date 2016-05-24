unbias.eigens <-
function(L,g,w,minalpha=NULL)
{
  # choose adjustment for alpha should L be too small to "de-bias"
  if(length(minalpha)==0) {minalpha = sqrt(max(g))}
  # get adjusted first eigens:
  a = c()
  for(k in 1:2)
  {
    #if(L[k]<(1+sqrt(g[k]))^2) {a[k] = 1+sqrt(g[k])} #old
    if(L[k]<(1+sqrt(g[k]))^2) {a[k] = 1+sqrt(max(g))+minalpha} #new
    if(L[k]>=(1+sqrt(g[k]))^2) {a[k] = ((1+L[k]-g[k]) + sqrt((1+L[k]-g[k])^2 - 4*L[k]))/2}
  }
  # get common 1st eigen, then bias correction	
  a0 = sum(a*w)
  a0 = max(c(a0,1+sqrt(g)))
  QLcorrection = (g[1]-g[2])*a0/(a0-1)	
  return(list(QLcorrection=QLcorrection,a0=a0,a=a))
}
