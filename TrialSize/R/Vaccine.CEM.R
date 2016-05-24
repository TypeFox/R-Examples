Vaccine.CEM <-
function(alpha,beta,mu_t,mu_c,sigma_t,sigma_c,pt,pc){

mu_star=(mu_t+mu_c)/2
p_star=(pt+pc)/2

n=(1/(mu_t*pt-mu_c*pc)^2)*(qnorm(1-alpha/2)*sqrt(2*mu_star^2*(p_star*(1-p_star))+2*p_star*(sigma_t^2+sigma_c^2))
+qnorm(1-beta)*sqrt(pt*(sigma_t^2+mu_t^2*(1-pt))+pc*(sigma_c^2+mu_c^2*(1-pc))))^2

}
