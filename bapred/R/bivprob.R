bivprob <-
function (rho,lower,upper=-lower,mean=0) {
  ##require("mnormt")
  nu = 0
  low = rep(as.double((lower - mean)),2)
  upp = rep(as.double((upper - mean)),2)
  if (any(lower == upper)) return(0)
  infin = c(2, 2)
  infin = as.integer(infin)
  low = replace(low, low == -Inf, 0)
  upp = replace(upp, upp == Inf, 0)
  rho = as.double(rho)
  prob = as.double(0)
  a = lapply(rho,function(r,low,upp) mnormt::biv.nt.prob(df=Inf,lower=low,upper=upp,mean=rep(0,2),S=matrix(c(1,r,r,1),2,2)),
             low=low,upp=upp)
  return(unlist(a))
}
