pbd_loglik_rhs = function(t,x,pars)
{
   b = pars[[1]](t,as.numeric(pars[5:length(pars)]))
   mu_g = pars[[2]](t,as.numeric(pars[5:length(pars)]))
   la2 = pars[[3]](t,as.numeric(pars[5:length(pars)]))
   mu_i = pars[[4]](t,as.numeric(pars[5:length(pars)]))
   nu2 = la2 + mu_i
   H = x[1]
   p_i = x[2]
   q_i = x[3]
   q_g = x[4]
   dH = -H * b * (1 - p_i)
   dp_i = -(nu2 + b) * p_i + la2 * q_g + mu_i + b*p_i^2
   dq_i = -(nu2 + b) * q_i + la2 * q_g + mu_i + b*q_i^2
   dq_g = -(mu_g + b) * q_g + mu_g + b*q_i*q_g
   dx = c(dH,dp_i,dq_i,dq_g)
   return(list(dx))
}
