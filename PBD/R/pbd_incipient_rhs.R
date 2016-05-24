pbd_incipient_rhs = function(t,x,pars)
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
   alpha = x[5]
   B = x[6]
   C = x[7]
   D = x[8]
   E = x[9]
   dH = -H * b * (1 - p_i)
   dp_i = -(nu2 + b) * p_i + la2 * q_g + mu_i + b*p_i^2
   dq_i = -(nu2 + b) * q_i + la2 * q_g + mu_i + b*q_i^2
   dq_g = -(mu_g + b) * q_g + mu_g + b*q_i*q_g
   dalpha = (b * q_i - nu2) * alpha
   dB = -(p_i - q_i) * b * B
   dC = -(1 - q_i) * b * C
   dD = -(1 - q_i - alpha) * b * D
   dE = (1 - q_i - alpha) * b * C
   dx = c(dH,dp_i,dq_i,dq_g,dalpha,dB,dC,dD,dE)
   return(list(dx))
}
