chao.sd <-
function(x) {
  f.one <- sum(x[x==1])
  f.two <- sum(x[x==2])
  gamma <- f.two*(0.25*(f.one/f.two)^4+(f.one/f.two)^3+0.5*(f.one/f.two)^2)
  return(sqrt(gamma))
}