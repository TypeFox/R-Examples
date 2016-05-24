dd_lamuN = function(ddmodel,pars,N)
{
    la = pars[1]
    mu = pars[2]
    K = pars[3]
    n0 = (ddmodel == 2 | ddmodel == 4)
    if(length(pars) == 4)
    {
        r = pars[4]
    }
    if(ddmodel == 1)
    {
        # linear dependence in speciation rate
        laN = max(0,la - (la - mu) * N/K)
        muN = mu
    }
    if(ddmodel == 1.3)
    {
        # linear dependence in speciation rate
        laN = max(0,la * (1 - N/K))
        muN = mu
    }
    if(ddmodel == 2 | ddmodel == 2.1 | ddmodel == 2.2)
    {
        # exponential dependence in speciation rate
        al = (log(la/mu)/log(K+n0))^(ddmodel != 2.2)
        laN = la * (N + n0)^(-al)
        muN = mu
    }
    if(ddmodel == 2.3)
    {
        # exponential dependence in speciation rate
        al = K
        laN = la * (N + n0)^(-al)
        muN = mu
    }
    if(ddmodel == 3)
    {
        # linear dependence in extinction rate
        laN = la
        muN = mu + (la - mu) * N/K
    }
    if(ddmodel == 4 | ddmodel == 4.1 | ddmodel == 4.2)
    {
        # exponential dependence in extinction rate
        al = (log(la/mu)/log(K+n0))^(ddmodel != 4.2)
        laN = la
        muN = mu * (N + n0)^al
    }
    if(ddmodel == 5)
    {
        # linear dependence in speciation rate and extinction rate
        laN = max(0,la - 1/(r+1)*(la-mu) * N/K)
        muN = mu + r/(r+1)*(la-mu)/K * N
    }
    return(c(laN,muN))
}

dd_sim = function(pars,age,ddmodel = 1)
{
# Simulation of diversity-dependent process
#  . start from crown age
#  . no additional species at crown node
#  . no missing species in present
# pars = [la mu K]
#  . la = speciation rate per species
#  . mu = extinction rate per species
#  . K = diversification carrying capacity
#  . r = ratio of diversity-dependence in extinction rate over speciation rate
# age = crown age
# ddmodel = mode of diversity-dependence
#  . ddmodel == 1 : linear dependence in speciation rate with parameter K
#  . ddmodel == 1.3: linear dependence in speciation rate with parameter K'
#  . ddmodel == 2 : exponential dependence in speciation rate
#  . ddmodel == 2.1: variant with offset at infinity
#  . ddmodel == 2.2: 1/n dependence in speciation rate
#  . ddmodel == 2.3: exponential dependence in speciation rate with parameter x
#  . ddmodel == 3 : linear dependence in extinction rate
#  . ddmodel == 4 : exponential dependence in extinction rate
#  . ddmodel == 4.1: variant with offset at infinity
#  . ddmodel == 4.2: 1/n dependence in speciation rate
#  . ddmodel == 5 : linear dependence in speciation rate and in extinction rate

done = 0
while(done == 0)
{
    # number of species N at time t
    # i = index running through t and N
    t = rep(0,1)
    L = matrix(0,2,4)
    i = 1
    t[1] = 0
    N = 2
    # L = data structure for lineages,
    # . L[,1] = branching times
    # . L[,2] = index of parent species
    # . L[,3] = index of daughter species
    # . L[,4] = time of extinction
    # j = index running through L
    L[1,1:4] = c(0,0,-1,-1)
    L[2,1:4] = c(0,-1,2,-1)
    linlist = c(-1,2)
    newL = 2
    ff = dd_lamuN(ddmodel,pars,N[i])
    laN = ff[1]
    muN = ff[2]
    denom = (laN + muN) * N[i]
    t[i + 1] = t[i] - log(runif(1)) / denom
    while(t[i + 1] <= age)
    {
        i = i + 1
        ranL = sample2(linlist,1)
        if((laN * N[i - 1] / denom) >= runif(1))
        {
            # speciation event
            N[i] = N[i - 1] + 1
            newL = newL + 1
            L = rbind(L,c(t[i],ranL,sign(ranL) * newL,-1))
            linlist = c(linlist,sign(ranL) * newL)
        } else {
            # extinction event
            N[i] = N[i - 1] - 1
            L[abs(ranL),4] = t[i]
            w = which(linlist == ranL)
            linlist = linlist[-w]
            linlist = sort(linlist)
        }
        if(sum(linlist < 0) == 0 | sum(linlist > 0) == 0)
        {
            t[i + 1] = Inf
        } else {
            ff = dd_lamuN(ddmodel,pars,N[i])
            laN = ff[1]
            muN = ff[2]
            denom = (laN + muN) * N[i]
            t[i + 1] = t[i] - log(runif(1)) / denom
        } 
    }
    if(sum(linlist < 0) == 0 | sum(linlist > 0) == 0)
    {
       done = 0
    } else {
       done = 1
    }
}

L[,1] = age - c(L[,1])
notmin1 = which(L[,4] != -1)
L[notmin1,4] = age - c(L[notmin1,4])
L[which(L[,4] == age + 1),4] = -1
tes = L2phylo(L,dropextinct = T)
tas = L2phylo(L,dropextinct = F)
out = list(tes = tes,tas = tas,L = L)
return(out)

}