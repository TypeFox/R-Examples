InterSV.NIS <-
function(alpha, beta,vbt,vwt,vbr,vwr,m,margin){

sigma=2*((vbt+vwt/m)^2+margin^4*(vbr+vwr/m)^2+vwt^2/(m^2*(m-1))+margin^4*vwr^2/(m^2*(m-1)))

n=sigma*(qnorm(1-alpha)+qnorm(1-beta))^2/(vbt-margin^2*vbr)^2
n
}
