Vaccine.ELDI <-
function(alpha,beta,theta0,theta,pt,pc){
#theta=pt/(pt+pc)
n=(qnorm(1-alpha)*sqrt(theta0*(1-theta0))+qnorm(1-beta)*sqrt(theta*(1-theta)))^2/((pt+pc)*(theta-theta0)^2)
}
