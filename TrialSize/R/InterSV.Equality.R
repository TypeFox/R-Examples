InterSV.Equality <-
function(alpha, beta,vbt,vwt,vbr,vwr,m){

sigma=2*(vbt+vwt/m)^2+(vbr+vwr/m)^2+vwt^2/(m^2*(m-1))+vwr^2/(m^2*(m-1))

n=sigma*(qnorm(1-alpha/2)+qnorm(1-beta))^2/(vbt-vbr)^2
n
}
