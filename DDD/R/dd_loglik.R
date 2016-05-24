# brts = branching times (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# - pars1[1] = la = (initial) speciation rate
# - pars1[2] = mu = extinction rate
# - pars1[3] = K = carrying capacity
# - pars1[4] = r = ratio of diversity-dependence in extinction rate over speciation rate
# - pars2[1] = lx = length of ODE variable x
# - pars2[2] = ddep = diversity-dependent model,mode of diversity-dependence
#  . ddep == 1 : linear dependence in speciation rate with parameter K
#  . ddep == 1.3 : linear dependence in speciation rate with parameter K'
#  . ddep == 2 : exponential dependence in speciation rate
#  . ddep == 2.1: variant with offset at infinity
#  . ddep == 2.2: 1/n dependence in speciation rate
#  . ddep == 2.3: exponential dependence in speciation rate with parameter x
#  . ddep == 3 : linear dependence in extinction rate
#  . ddep == 4 : exponential dependence in extinction rate
#  . ddep == 4.1: variant with offset at infinity
#  . ddep == 4.2: 1/n dependence in speciation rate
#  . ddep == 5 : linear dependence in speciation and extinction rate
# - pars2[3] = cond = conditioning
#  . cond == 0 : conditioning on stem or crown age
#  . cond == 1 : conditioning on stem or crown age and non-extinction of the phylogeny
#  . cond == 2 : conditioning on stem or crown age and on the total number of extant taxa (including missing species)
#  . cond == 3 : conditioning on the total number of extant taxa (including missing species)
# - pars2[4] = btorph = likelihood of branching times (0) or phylogeny (1), differ by a factor (S - 1)! where S is the number of extant species
# - pars2[5] = parameters and likelihood should be printed (1) or not (0)
# - pars2[6] = likelihood is for a tree with crown age (2) or stem age (1)
# missnumspec = number of missing species    
# methode = the method used in the numerical solving of the set of the ode's or 'analytical' which means matrix exponentiation is used

dd_loglik = function(pars1,pars2,brts,missnumspec,methode = 'analytical')
{
   if(methode == 'analytical')
   {
       out = dd_loglik2(pars1,pars2,brts,missnumspec)
   } else {
       out = dd_loglik1(pars1,pars2,brts,missnumspec,methode = methode)
   }
   return(out)
}

dd_loglik1 = function(pars1,pars2,brts,missnumspec,methode = 'lsoda')
{
if(length(pars2) == 4)
{
    pars2[5] = 0
    pars2[6] = 2
}
ddep = pars2[2]
cond = pars2[3]
btorph = pars2[4]
soc = pars2[6]
if(cond == 3) { soc = 2 }
la = pars1[1]
mu = pars1[2]
K = pars1[3]
if(ddep == 5) {r = pars1[4]} else {r = 0}
if(ddep == 1 | ddep == 5)
{ 
    lx = min(max(1 + missnumspec,1 + ceiling(la/(la - mu) * (r + 1) * K)),ceiling(pars2[1]))
} else {
if(ddep == 1.3)
{
    lx = min(ceiling(K),ceiling(pars2[1]))
} else {
    lx = round(pars2[1])
}}
n0 = (ddep == 2 | ddep == 4)
if((ddep == 1) & ((mu == 0 & missnumspec == 0 & floor(K) != ceiling(K) & la > 0.05) | K == Inf))
{
    loglik = bd_loglik(pars1[1:(2 + (K < Inf))],c(2*(mu == 0 & K < Inf),pars2[3:6]),brts,missnumspec)
} else {
abstol = 1e-16
reltol = 1e-10 
brts = -sort(abs(as.numeric(brts)),decreasing = TRUE)
if(sum(brts == 0) == 0)
{
   brts[length(brts) + 1] = 0
}

S = length(brts) + (soc - 2)
if(min(pars1) < 0)
{
    loglik = -Inf
} else {
if((mu == 0 & (ddep == 2 | ddep == 2.1 | ddep == 2.2)) | (la == 0 & (ddep == 4 | ddep == 4.1 | ddep == 4.2)) | (la <= mu))
{ 
    cat("These parameter values cannot satisfy lambda(N) = mu(N) for a positive and finite N.\n")
    loglik = -Inf
} else {
    if(((ddep == 1 | ddep == 5) & ceiling(la/(la - mu) * (r + 1) * K) < (S + missnumspec)) | ((ddep == 1.3) & (S + missnumspec > ceiling(K))))
    {
       loglik = -Inf
    } else {
       loglik = (btorph == 0) * lgamma(S)
       if(cond != 3)
       {
          probs = rep(0,lx)
          probs[1] = 1 # change if other species at stem/crown age 
          for(k in 2:(S + 2 - soc))
          {
             k1 = k + (soc - 2)
             y = ode(probs,brts[(k-1):k],dd_loglik_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
             probs = y[2,2:(lx+1)]
             if(k < (S + 2 - soc))
             {
                 probs = flavec(ddep,la,mu,K,r,lx,k1,n0) * probs # speciation event
                 if(sum(probs) <= 0)
                 {
                    loglik = -Inf
                    break
                 } else {
                    loglik = loglik + log(sum(probs))
                 }
                 probs = probs/sum(probs)
             }
          }    
       } else {
          probs = rep(0,lx + 1)
          probs[1 + missnumspec] = 1
          for(k in (S + 2 - soc):2)
          {
             k1 = k + (soc - 2)
             y = ode(probs,-brts[k:(k-1)],dd_loglik_bw_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
             probs = y[2,2:(lx+2)]
             if(k > soc)
             {
                 probs = c(flavec(ddep,la,mu,K,r,lx,k1-1,n0),1) * probs # speciation event
                 if(sum(probs[1:lx]) <= 0)
                 {
                    loglik = -Inf
                    break
                 } else {
                    loglik = loglik + log(sum(probs[1:lx]))
                 }
                 probs[1:lx] = probs[1:lx]/sum(probs[1:lx])
             }    
          }
       }
       if(probs[1 + missnumspec] <= 0 | loglik == -Inf)
       {
          loglik = -Inf
       } else  {        
          loglik = loglik + (cond != 3 | soc == 1) * log(probs[1 + (cond != 3) * missnumspec]) - lgamma(S + missnumspec + 1) + lgamma(S + 1) + lgamma(missnumspec + 1)
  
          logliknorm = 0
          if(cond == 1 | cond == 2)
          {
             probsn = rep(0,lx)
             probsn[1] = 1 # change if other species at stem or crown age
             k = soc
             t1 = brts[1] 
             t2 = brts[S + 2 - soc]
             y = ode(probsn,c(t1,t2),dd_loglik_rhs,c(pars1,k,ddep),rtol = reltol,atol = abstol,method = methode);
             probsn = y[2,2:(lx+1)]
             if(soc == 1) { aux = 1:lx }
             if(soc == 2) { aux = (2:(lx+1)) * (3:(lx+2))/6 }
             probsc = probsn/aux
             if(cond == 1) { logliknorm = log(sum(probsc)) }
             if(cond == 2) { logliknorm = log(probsc[S + missnumspec - soc + 1])}             
          }
          if(cond == 3)
          { 
             probsn = rep(0,lx + 1)
             probsn[S + missnumspec + 1] = 1
             TT = max(1,1/abs(la - mu)) * 1E+10 * max(abs(brts)) # make this more efficient later
             y = ode(probsn,c(0,TT),dd_loglik_bw_rhs,c(pars1,0,ddep),rtol = reltol,atol = abstol,method = methode)
             logliknorm = log(y[2,lx + 2])
             if(soc == 2)
             {
                probsn = rep(0,lx + 1)
                probsn[1:lx] = probs[1:lx]
                probsn = c(flavec(ddep,la,mu,K,r,lx,1,n0),1) * probsn # speciation event
                y = ode(probsn,c(max(abs(brts)),TT),dd_loglik_bw_rhs,c(pars1,1,ddep),rtol = reltol,atol = abstol,method = methode)
                logliknorm = logliknorm - log(y[2,lx + 2])
             }
          }
          loglik = loglik - logliknorm
       }
    }
}}
if(pars2[5] == 1)
{
    s1 = sprintf('Parameters: %f %f %f',pars1[1],pars1[2],pars1[3])
    if(ddep == 5) {s1 = sprintf('%s %f',s1,pars1[4])}
    s2 = sprintf(', Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    flush.console()
}
}
loglik = as.numeric(loglik)
if(is.nan(loglik) | is.na(loglik))
{
    loglik = -Inf
}
return(loglik)
}


dd_loglik2 = function(pars1,pars2,brts,missnumspec,methode = 'lsoda')
{
if(length(pars2) == 4)
{
    pars2[5] = 0
    pars2[6] = 2
}
ddep = pars2[2]
cond = pars2[3]
btorph = pars2[4]
soc = pars2[6]
if(cond == 3)
{ 
    soc = 2
}
la = pars1[1]
mu = pars1[2]
K = pars1[3]
if(ddep == 5)
{
    r = pars1[4]
} else
{
    r = 0
}
if(ddep == 1 | ddep == 5)
{ 
    lx = min(max(1 + missnumspec,1 + ceiling(la/(la - mu) * (r + 1) * K)),ceiling(pars2[1]))
} else if(ddep == 1.3)
{
    lx = min(ceiling(K),ceiling(pars2[1]))
} else {
    lx = round(pars2[1])
}
n0 = (ddep == 2 | ddep == 4)
if((ddep == 1) & ((mu == 0 & missnumspec == 0 & floor(K) != ceiling(K) & la > 0.05) | K == Inf))
{
    loglik = bd_loglik(pars1[1:(2 + (K < Inf))],c(2*(mu == 0 & K < Inf),pars2[3:6]),brts,missnumspec)
} else {
abstol = 1e-16
reltol = 1e-10 
brts = -sort(abs(as.numeric(brts)),decreasing = TRUE)
if(sum(brts == 0) == 0)
{
   brts[length(brts) + 1] = 0
}
S = length(brts) + (soc - 2)
if(min(pars1) < 0)
{
    loglik = -Inf
} else {
if((mu == 0 & (ddep == 2 | ddep == 2.1 | ddep == 2.2)) | (la == 0 & (ddep == 4 | ddep == 4.1 | ddep == 4.2)) | (la <= mu))
{ 
    cat("These parameter values cannot satisfy lambda(N) = mu(N) for a positive and finite N.\n")
    loglik = -Inf
} else {
    if(((ddep == 1 | ddep == 5) & ceiling(la/(la - mu) * (r + 1) * K) < (S + missnumspec)) | ((ddep == 1.3) & ((S + missnumspec) > ceiling(K))))
    {
       loglik = -Inf
    } else {
       loglik = (btorph == 0) * lgamma(S)
       if(cond != 3)
       {
          probs = rep(0,lx)
          probs[1] = 1 # change if other species at stem/crown age 
          for(k in 2:(S + 2 - soc))
          {
             k1 = k + (soc - 2)
             #y = ode(probs,brts[(k-1):k],dd_loglik_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
             #probs2 = y[2,2:(lx+1)]
             probs = dd_loglik_M(pars1,lx,k1,ddep,tt = abs(brts[k] - brts[k-1]),probs)
             #print(as.numeric(probs[1:10]))
             #print(as.numeric(probs2[1:10])) 
             if(k < (S + 2 - soc))
             {
                 #probs = flavec(ddep,la,mu,K,r,lx,k1,n0) * probs # speciation event
                 probs = lambdamu(0:(lx - 1) + k1,c(pars1[1:3],r),ddep)[[1]] * probs
                 if(sum(probs) <= 0 | sum(is.na(probs)) > 0 | sum(is.nan(probs)) > 0)
                 {
                    loglik = -Inf
                    break
                 } else {
                    loglik = loglik + log(sum(probs))
                 }
                 probs = probs/sum(probs)
             }
          }    
       } else {
          probs = rep(0,lx + 1)
          probs[1 + missnumspec] = 1
          for(k in (S + 2 - soc):2)
          {
             k1 = k + (soc - 2)
             #y = ode(probs,-brts[k:(k-1)],dd_loglik_bw_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
             #probs2 = y[2,2:(lx+2)]
             probs = dd_loglik_M_bw(pars1,lx,k1,ddep,tt = abs(brts[k] - brts[k-1]),probs[1:lx])
             probs = c(probs,0)
             if(k > soc)
             {
                 #probs = c(flavec(ddep,la,mu,K,r,lx,k1-1,n0),1) * probs # speciation event
                 probs = c(lambdamu(0:(lx - 1) + k1 - 1,pars1,ddep)[[1]],1) * probs
                 if(sum(probs[1:lx]) <= 0 | sum(is.na(probs[1:lx])) > 0 | sum(is.nan(probs[1:lx])) > 0)
                 {
                    loglik = -Inf
                    break
                 } else {
                    loglik = loglik + log(sum(probs[1:lx]))
                 }
                 probs[1:lx] = probs[1:lx]/sum(probs[1:lx])
             }    
          }
       }
       if(probs[1 + (cond != 3) * missnumspec] <= 0 | loglik == -Inf)
       {
          loglik = -Inf
       } else  {        
          loglik = loglik + (cond != 3 | soc == 1) * log(probs[1 + (cond != 3) * missnumspec]) - lgamma(S + missnumspec + 1) + lgamma(S + 1) + lgamma(missnumspec + 1)
  
          logliknorm = 0
          if(cond == 1 | cond == 2)
          {
             probsn = rep(0,lx)
             probsn[1] = 1 # change if other species at stem or crown age
             k = soc
             t1 = brts[1] 
             t2 = brts[S + 2 - soc]
             #y = ode(probsn,c(t1,t2),dd_loglik_rhs,c(pars1,k,ddep),rtol = reltol,atol = abstol,method = methode);
             #probsn = y[2,2:(lx+1)]
             probsn = dd_loglik_M(pars1,lx,k,ddep,tt = abs(t2 - t1),probsn)
             if(soc == 1) { aux = 1:lx }
             if(soc == 2) { aux = (2:(lx+1)) * (3:(lx+2))/6 }
             probsc = probsn/aux
             if(cond == 1) { logliknorm = log(sum(probsc)) }
             if(cond == 2) { logliknorm = log(probsc[S + missnumspec - soc + 1])}             
          }
          if(cond == 3)
          { 
             #probsn = rep(0,lx + 1)
             #probsn[S + missnumspec + 1] = 1 #/ (S + missnumspec)
             #TT = max(1,1/abs(la - mu)) * 100000000 * max(abs(brts)) # make this more efficient later
             #y = ode(probsn,c(0,TT),dd_loglik_bw_rhs,c(pars1,0,ddep),rtol = reltol,atol = abstol,method = methode)
             #logliknorm = log(y[2,lx + 2])
             probsn = rep(0,lx + 1)
             probsn[2] = 1
             MM = dd_loglik_M_aux(pars1,lx + 1,k = 0,ddep)
             MM = MM[-1,-1]
             #probsn = SparseM::solve(-MM,probsn[2:(lx + 1)])
             MMinv = SparseM::solve(MM)
             probsn = -MMinv %*% probsn[2:(lx + 1)]
             logliknorm = log(probsn[S + missnumspec])
             if(soc == 2)
             {
                #probsn = rep(0,lx + 1)
                #probsn[1:lx] = probs[1:lx]
                #probsn = c(flavec(ddep,la,mu,K,r,lx,1,n0),1) * probsn # speciation event
                #probsn = c(lambdamu(0:(lx - 1) + 1,pars1,ddep)[[1]],1) * probsn # speciation event
                #y = ode(probsn,c(max(abs(brts)),TT),dd_loglik_bw_rhs,c(pars1,1,ddep),rtol = reltol,atol = abstol,method = methode)
                #logliknorm = logliknorm - log(y[2,lx + 2])
                probsn2 = rep(0,lx)
                probsn2 = lambdamu(0:(lx - 1) + 1,pars1,ddep)[[1]] * probs[1:lx]
                MM = dd_loglik_M_bw_aux(pars1,lx,k = 1,ddep)
                MMinv = SparseM::solve(MM)
                #probsn2 = SparseM::solve(-MM,probsn2[1:lx])
                probsn2 = -MMinv %*% probsn2[1:lx]
                logliknorm = logliknorm - log(probsn2[1])            
             }
          }
          loglik = loglik - logliknorm
       }
    }
}}
if(pars2[5] == 1)
{
    s1 = sprintf('Parameters: %f %f %f',pars1[1],pars1[2],pars1[3])
    if(ddep == 5) {s1 = sprintf('%s %f',s1,pars1[4])}
    s2 = sprintf(', Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    flush.console()
}
}
loglik = as.numeric(loglik)
if(is.nan(loglik) | is.na(loglik) | loglik == Inf)
{
    loglik = -Inf
}
return(loglik)
}