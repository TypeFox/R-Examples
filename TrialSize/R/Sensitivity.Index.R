Sensitivity.Index <-
function(alpha,n,deltaT){
t=qt(1-alpha/2,n-2)
p=1-pt(t,n-2,deltaT)+pt(-t,n-2,deltaT)
}
