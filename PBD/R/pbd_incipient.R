pbd_incipient = function(pars1,pars1f = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},function(t,pars) {pars[3]},function(t,pars) {pars[4]}), pars2 = c(1,1,2,1,"lsoda",1000),brts,missnumspec = 0)
{
# pbd_cryptic computes the distribution of incipient species under the protracted speciation model given a set of branching times and data

# pars1 contains model parameters
# In the simplest case where rates do not depend on time, we have
# - pars1[1] = b (= la_1 in ER2012) = speciation initiation rate
# - pars1[2] = mu_1 (= mu_g in ER2012) = extinction rate of good species
# - pars1[3] = la_1 (= la_2 in ER2012) = speciation completion rate
# - pars1[4] = mu_2 (= mu_i in ER2012) = extinction rate of incipient species
# When rates depend on time this time dependence should be specified in pars1f and pars1 then become the parameters used in pars1f
# pars1f contains the functions how the rates depend on time, default functions are constant functions of the parameters in pars1
# pars2 contains settings
# - pars2[1] = cond = conditioning on age (0), age and non-extinction of the clade (1) or age and number of extant taxa (2)
# - pars2[2] = btorph = likelihood for branching times (0) or phylogeny (1)
# - pars2[3] = soc = stem (1) or crown (2) age
# - pars2[4] = printing of parameters and likelihood (1) or not (0)
# - pars2[5] = method of the numerical integration; see package deSolve for details
# - pars2[6] = resolution of the probability distribution
# brts = set of branching times
# missnumspec = the number of species that belong to the same clade but are not in the phylogeny

# Example: pbd_incipient(pars1 = c(0.1,0.05,1,0.05), brts = 1:10, missnumspec = 4)

pars1 = c(pars1f,pars1)

brts = sort(abs(brts))
abstol = 1e-16
reltol = 1e-10 
b = pars1[[1]](brts,as.numeric(pars1[5:length(pars1)]))
methode = pars2[5]
cond = as.numeric(pars2[1])
btorph = as.numeric(pars2[2])
soc = as.numeric(pars2[3])
res = as.numeric(pars2[6])
S = length(brts) + (soc - 1)
m = missnumspec


probs = c(1,1,0,0,1,1,1,1,0)
y = ode(probs,c(0,brts),pbd_incipient_rhs,c(pars1),rtol = reltol,atol = abstol,method = methode)
#loglik = (btorph == 0) * lgamma(S) + sum(log(b) + log(y[2:S,2]) + log(1 - y[2:S,3])) - log(b[1]) + (soc == 2) * (log(y[S,2]) + log(1 - y[S,3])) - soc * (cond > 0) * (log(1 - y[S,3])) - (cond == 2) * ((soc == 2) * log(S + m - 1) + soc * log(y[S,2]) + (S + m - soc) * log(1 - y[S,2]))
print(y)
A = y[2:S,10]/(1 - y[2:S,9])
print(A)
B = y[2:S,7]
print(B)
betapar = 1 - y[2:S,4] - y[2:S,6]
gammapar = 1 - y[2:S,3]
#print(betapar/gammapar)
PN = rep(res + 1,0)
PN[1] = 1
for(j in 1:(S - 1))
{
    print(j)
    PLj = c(betapar[j]/gammapar[j],(1 - betapar[j]/gammapar[j]) * A[j] * (1 - A[j])^((1:res) - 1))
    #print(PLj[1:10])
    #print(sum(PLj[2:length(PLj)]))
    PRj = B[j] * (1 - B[j])^(0:res)
    #print(PLj[1:10])
    #print(sum(PRj))
    PNj = conv(PLj,PRj)[1:(res + 1)]
    #print(PNj[1:10])
    #print(sum(PNj))
    PN = conv(PN,PNj)[1:(res + 1)]
    print(PN[1:10])
}
PN = conv(PN,PNj)[1:(res + 1)]
print(PN[1:10])
expinc = t(PN) %*% (0:res)

#if(m > 0)
#{
#   if(soc == 1)
#   {
#      y2 = as.numeric(c(1-y[2:S,2]))
#   }
#   if(soc == 2)
#   {
#      y2 = as.numeric(c(1-y[2:S,2],1-y[S,2]))
#   }
#   x = rep(0,m + 1)
#   x[1] = 1
#   for(j in 1:S)
#   {
#       #x = convolve(x,rev((1:(m + 1)) * (y2[j]^(0:m))),type = 'open')[1:(m + 1)]
#       x = conv(x,(1:(m + 1)) * (y2[j]^(0:m)))[1:(m+1)]
#   }
#   loglik = loglik + lgamma(S + 1) + lgamma(m + 1) - lgamma(S + m + 1) + log(x[m + 1])
#}

if(as.numeric(pars2[4]) == 1)
{
    pastetxt = paste('Parameters:',pars1[[5]][1],sep = ' ')
    for(i in 6:length(pars1))
    {
        pastetxt = paste(pastetxt,pars1[[i]][1],sep = ', ')
    }
    s2 = sprintf(', Expected number of incipient speces: %f',expinc)
    cat(pastetxt,s2,"\n",sep = "")
    flush.console()
}

return(as.numeric(PN))
}