pbd_brts_density = function(pars1,pars1f = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},function(t,pars) {pars[3]},function(t,pars) {pars[4]}), methode = "lsoda",brts)
{
# pbd_brts_density computes the density of node depth under the protracted speciation model

# pars1 contains model parameters
# In the simplest case where rates do not depend on time, we have
# - pars1[1] = b (= la_1 in ER2012) = speciation initiation rate
# - pars1[2] = mu_1 (= mu_g in ER2012) = extinction rate of good species
# - pars1[3] = la_1 (= la_2 in ER2012) = speciation completion rate
# - pars1[4] = mu_2 (= mu_i in ER2012) = extinction rate of incipient species
# When rates depend on time this time dependence should be specified in pars1f and pars1 then become the parameters used in pars1f
# pars1f contains the functions how the rates depend on time, default functions are constant functions of the parameters in pars1
# methode gives method of the numerical integration; see package deSolve for details
# brts = set of branching times for which the density needs to be computed

pars1 = c(pars1f,pars1)

brts = sort(abs(brts))
abstol = 1e-16
reltol = 1e-10 
b = pars1[[1]](brts,as.numeric(pars1[5:length(pars1)]))
S = (length(brts) + 1)

probs = c(1,1,0,0)
y = ode(probs,c(0,brts),pbd_loglik_rhs,c(pars1),rtol = reltol,atol = abstol,method = methode)
dens = b * y[2:S,2] * (1 - y[2:S,3])
}