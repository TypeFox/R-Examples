pbd_sim_cpp = function(pars,parsf = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},function(t,pars) {pars[3]},function(t,pars) {pars[4]}), age, soc = 2, plotltt = 1, methode = "lsoda")
{    

pbd_loglik_rhs_cpp = function(t,x,pars)
{
   b = pars[[1]](t,as.numeric(pars[5:(length(pars)-1)]))
   mu_g = pars[[2]](t,as.numeric(pars[5:(length(pars)-1)]))
   la2 = pars[[3]](t,as.numeric(pars[5:(length(pars)-1)]))
   mu_i = pars[[4]](t,as.numeric(pars[5:(length(pars)-1)]))
   nu2 = la2 + mu_i
   H = x[1]
   p_i = x[2]
   q_i = x[3]
   q_g = x[4]
   dH = -H * b * (1 - p_i)
   dp_i = -(nu2 + b) * p_i + la2 * q_g + mu_i + b*p_i^2
   dq_i = -(nu2 + b) * q_i + la2 * q_g + mu_i + b*q_i^2
   dq_g = -(mu_g + b) * q_g + mu_g + b*q_i*q_g
   dt = (x[1] >= pars[length(pars)])
   dx = c(dH,dp_i,dq_i,dq_g,dt)
   return(list(dx))
}

pars1 = c(parsf,pars)
abstol = 1e-16
reltol = 1e-10 
probs = c(1,1,0,0,0)
y = ode(probs,c(0,age),pbd_loglik_rhs_cpp,c(pars1,0),rtol = reltol,atol = abstol,method = methode)
pT = 1 - y[2,2]
#print(pT)
nd = sum(rgeom(soc,1 - pT))
#print(nd)
brts = rep(0,nd + 1)
brts[1] = age
i = 1
while(i <= nd)
{
   probs = c(1,1,0,0,0)
   y = ode(probs,c(0,age),pbd_loglik_rhs_cpp,c(pars1,1 - runif(1)*pT),rtol = reltol,atol = abstol,method = methode)
   brts[i + 1] = y[2,6]
   i = i + 1
}
brts = sort(brts, decreasing = T)
if(plotltt == 1)
{
    plot(c(-brts,0),c(soc:(length(brts)+soc-1),length(brts)+soc-1),log = 'y',type = 's',xlab = 'Time',ylab = 'Number of lineages')
}
return(brts)
}