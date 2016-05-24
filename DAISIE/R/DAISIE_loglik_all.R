DAISIE_loglik_rhs = function(t,x,pars)
{
lx = (length(x) - 1)/2
lac = pars[1]
mu = pars[2]
K = pars[3]
gam = pars[4]
laa = pars[5]
kk = pars[6]
ddep = pars[7]

nn = -2:(lx+2*kk+1)
lnn = length(nn)
nn = pmax(rep(0,lnn),nn)

if(ddep == 0)
{
    laavec = laa * rep(1,lnn)
    lacvec = lac * rep(1,lnn)
    muvec = mu * rep(1,lnn)
    gamvec = gam * rep(1,lnn)
} else {
if(ddep == 1)
{
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * (1 - nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = gam * rep(1,lnn)
} else {
if(ddep == 2)
{
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * exp(-nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = gam * rep(1,lnn)
} else {
if(ddep == 11)
{
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * (1 - nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = pmax(rep(0,lnn),gam * (1 - nn/K))
} else {
if(ddep == 21)
{
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * exp(-nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = pmax(rep(0,lnn),gam * exp(-nn/K))
} else {
if(ddep == 3)
{
    laavec = laa * rep(1,lnn)
    lacvec = lac * rep(1,lnn)
    muvec = mu * (1 + nn/K)
    gamvec = gam * rep(1,lnn)
}
}}}}}

#x = x * (x > 0)

xx1 = c(0,0,x[1:lx],0)
xx2 = c(0,0,x[(lx + 1):(2 * lx)],0)
xx3 = x[2 * lx + 1]

nil2lx = 3:(lx + 2)

il1 = nil2lx+kk-1
il2 = nil2lx+kk+1
il3 = nil2lx+kk
il4 = nil2lx+kk-2

in1 = nil2lx+2*kk-1
in2 = nil2lx+1
in3 = nil2lx+kk

ix1 = nil2lx-1
ix2 = nil2lx+1
ix3 = nil2lx
ix4 = nil2lx-2

dx1 = laavec[il1 + 1] * xx2[ix1] + lacvec[il4 + 1] * xx2[ix4] + muvec[il2 + 1] * xx2[ix3] +
      lacvec[il1] * nn[in1] * xx1[ix1] + muvec[il2] * nn[in2] * xx1[ix2] +
      -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
      -gamvec[il3] * xx1[ix3]
dx1[1] = dx1[1] + laavec[il3[1]] * xx3 * (kk == 1)
dx1[2] = dx1[2] + 2 * lacvec[il3[1]] * xx3 * (kk == 1)

dx2 = gamvec[il3] * xx1[ix3] +
      lacvec[il1 + 1] * nn[in1] * xx2[ix1] + muvec[il2 + 1] * nn[in2] * xx2[ix2] +
      -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
      -laavec[il3 + 1] * xx2[ix3]

dx3 = -(laavec[il3[1]] + lacvec[il3[1]] + gamvec[il3[1]] + muvec[il3[1]]) * xx3

return(list(c(dx1,dx2,dx3)))
}

DAISIE_loglik_rhs2 = function(t,x,pars)
{
lx = length(x)/3
lac = pars[1]
mu = pars[2]
K = pars[3]
gam = pars[4]
laa = pars[5]
kk = pars[6]
ddep = pars[7]

nn = -2:(lx+2*kk+1)
lnn = length(nn)
nn = pmax(rep(0,lnn),nn)

if(ddep == 0)
{
    laavec = laa * rep(1,lnn)
    lacvec = lac * rep(1,lnn)
    muvec = mu * rep(1,lnn)
    gamvec = gam * rep(1,lnn)
} else {
if(ddep == 1)
{
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * (1 - nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = gam * rep(1,lnn)
} else {
if(ddep == 2)
{
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * exp(-nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = gam * rep(1,lnn)
} else {
if(ddep == 11)
{
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * (1 - nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = pmax(rep(0,lnn),gam * (1 - nn/K))
} else {
if(ddep == 21)
{
    laavec = laa * rep(1,lnn)
    lacvec = pmax(rep(0,lnn),lac * exp(-nn/K))
    muvec = mu * rep(1,lnn)
    gamvec = pmax(rep(0,lnn),gam * exp(-nn/K))
} else {
if(ddep == 3)
{
    laavec = laa * rep(1,lnn)
    lacvec = lac * rep(1,lnn)
    muvec = mu * (1 + nn/K)
    gamvec = gam * rep(1,lnn)
}
}}}}}

#x = x * (x > 0)

xx1 = c(0,0,x[1:lx],0)
xx2 = c(0,0,x[(lx + 1):(2 * lx)],0)
xx3 = c(0,0,x[(2 * lx + 1):(3 * lx)],0)

nil2lx = 3:(lx + 2)

il1 = nil2lx+kk-1
il2 = nil2lx+kk+1
il3 = nil2lx+kk
il4 = nil2lx+kk-2

in1 = nil2lx+2*kk-1
in2 = nil2lx+1
in3 = nil2lx+kk
in4 = nil2lx-1

ix1 = nil2lx-1
ix2 = nil2lx+1
ix3 = nil2lx
ix4 = nil2lx-2

# inflow:
# anagenesis in colonist when k = 1: Q_M,n -> Q^1_n; n+k species present
# cladogenesis in colonist when k = 1: Q_M,n-1 -> Q^1_n; n+k-1 species present; rate twice
# anagenesis of reimmigrant: Q^M,k_n-1 -> Q^k,n; n+k-1+1 species present
# cladogenesis of reimmigrant: Q^M,k_n-2 -> Q^k,n; n+k-2+1 species present; rate once
# extinction of reimmigrant: Q^M,k_n -> Q^k,n; n+k+1 species present
# cladogenesis in one of the n+k-1 species: Q^k_n-1 -> Q^k_n; n+k-1 species present; rate twice for k species
# extinction in one of the n+1 species: Q^k_n+1 -> Q^k_n; n+k+1 species present
# outflow:
# all events with n+k species present
dx1 = (laavec[il3] * xx3[ix3] + 2 * lacvec[il1] * xx3[ix1]) * (kk == 1) + 
      laavec[il1 + 1] * xx2[ix1] + lacvec[il4 + 1] * xx2[ix4] + muvec[il2 + 1] * xx2[ix3] +
      lacvec[il1] * nn[in1] * xx1[ix1] + muvec[il2] * nn[in2] * xx1[ix2] +
      -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] - gamvec[il3] * xx1[ix3]

# inflow:
# immigration when there are n+k species: Q^k,n -> Q^M,k_n; n+k species present
# cladogenesis in n+k-1 species: Q^M,k_n-1 -> Q^M,k_n; n+k-1+1 species present; rate twice for k species
# extinction in n+1 species: Q^M,k_n+1 -> Q^M,k_n; n+k+1+1 species present
# outflow:
# all events with n+k+1 species present
dx2 = gamvec[il3] * xx1[ix3] +
      lacvec[il1 + 1] * nn[in1] * xx2[ix1] + muvec[il2 + 1] * nn[in2] * xx2[ix2] +
      -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
      -laavec[il3 + 1] * xx2[ix3]
      
# only when k = 1         
# inflow:
# cladogenesis in one of the n-1 species: Q_M,n-1 -> Q_M,n; n+k-1 species present; rate once
# extinction in one of the n+1 species: Q_M,n+1 -> Q_M,n; n+k+1 species present
# outflow:
# all events with n+k species present
dx3 = lacvec[il1] * nn[in4] * xx3[ix1] + muvec[il2] * nn[in2] * xx3[ix2] +
      -(lacvec[il3] + muvec[il3]) * nn[in3] * xx3[ix3] +
      -(laavec[il3] + gamvec[il3]) * xx3[ix3]
      
return(list(c(dx1,dx2,dx3)))
}

checkprobs = function(lx,loglik,probs)
{
   probs = probs * (probs > 0)
   if(sum(probs[1:(2 * lx)]) <= 0)
      {
         loglik = -Inf
      } else {
         loglik = loglik + log(sum(probs[1:(2 * lx)]))
         probs[1:(2 * lx)] = probs[1:(2 * lx)]/sum(probs[1:(2 * lx)])
      }
   return(list(loglik,probs))
}

checkprobs2 = function(lx,loglik,probs)
{
   probs = probs * (probs > 0)
   if(sum(probs) <= 0)
      {
         loglik = -Inf
      } else {
         loglik = loglik + log(sum(probs))
         probs = probs/sum(probs)
      }   
   return(list(loglik,probs))
}

divdepvec = function(lacgam,K,lx,k1,ddep)
{
   if(ddep == 1 | ddep == 11)
   {
	    vec = pmax(rep(0,lx + 1),lacgam * (1 - ((0:lx)+k1) / K))
   } else {
      if(ddep == 2 | ddep == 21)
      {
         vec = pmax(rep(0,lx + 1),lacgam * exp(-((0:lx)+k1) / K))
      } else {
         if(ddep == 0 | ddep == 3)
  		   {
		        vec = lacgam * rep(1,lx + 1)
         }
      }
   }
   return(vec)
}        

DAISIE_loglik = function(pars1,pars2,brts,stac,missnumspec,methode = "lsodes")
{
# brts = branching times (positive, from present to past)
# - max(brts) = age of the island
# - next largest brts = stem age / time of divergence from the mainland
# The interpretation of this depends on stac (see below)
# For stac = 0, there is no other value.
# For stac = 1, this is the time since divergence from the immigrant's sister on the mainland.
# The immigrant must have immigrated at some point since then.
# For stac = 2 and stac = 3, this is the time since divergence from the mainland.
# The immigrant that established the clade on the island must have immigrated precisely at this point.
# For stac = 3, it must have reimmigrated, but only after the first immigrant had undergone speciation.
# - min(brts) = most recent branching time (only for stac = 2, or stac = 3)
# pars1 = model parameters
# - pars1[1] = lac = (initial) cladogenesis rate
# - pars1[2] = mu = extinction rate
# - pars1[3] = K = maximum number of species possible in the clade
# - pars1[4] = gam = (initial) immigration rate
# - pars1[5] = laa = (initial) anagenesis rate
# pars2 = model settings
# - pars2[1] = lx = length of ODE variable x
# - pars2[2] = ddep = diversity-dependent model,mode of diversity-dependence
#  . ddep == 0 : no diversity-dependence
#  . ddep == 1 : linear dependence in speciation rate (anagenesis and cladogenesis)
#  . ddep == 11 : linear dependence in speciation rate and immigration rate
#  . ddep == 3 : linear dependence in extinction rate
# - pars2[3] = cond = conditioning
#  . cond == 0 : no conditioning
#  . cond == 1 : conditioning on presence on the island (not used in this single loglikelihood)
# - pars2[4] = parameters and likelihood should be printed (1) or not (0)
# stac = status of the clade formed by the immigrant
#  . stac == 0 : immigrant is not present and has not formed an extant clade
#  . stac == 1 : immigrant is present but has not formed an extant clade
#  . stac == 2 : immigrant is not present but has formed an extant clade
#  . stac == 3 : immigrant is present and has formed an extant clade
#  . stac == 4 : immigrant is present but has not formed an extant clade, and it is known when it immigrated.
# missnumspec = number of missing species

ddep = pars2[2]
cond = pars2[3]
if(cond > 0)
{
   cat("Conditioning has not been implemented and may not make sense. Cond is set to 0.\n")
}

lac = pars1[1]
mu = pars1[2]
K = pars1[3]
if(ddep == 0)
{
   K = Inf
}
gam = pars1[4]
laa = pars1[5]

abstol = 1e-16
reltol = 1e-10
brts = -sort(abs(as.numeric(brts)),decreasing = TRUE)
if(sum(brts == 0) == 0)
{
   brts[length(brts) + 1] = 0
}
# for stac = 0 and stac = 1, brts will contain origin of island and 0; length = 2; no. species should be 0
# for stac = 1, brts will contain origin of island and 0; length = 2; no. species should be 1
# for stac = 4, brts will contain origin of island, colonization event and 0; length = 3; no. species should be 1
# for stac = 2, brts with contain origin of island, colonization event, branching times, 0; no. species should be no. branching times + 1
# for stac = 3, brts with contain origin of island, colonization event, branching times, 0; no. species should be no. branching times + 2
#S = length(brts) - 2 * (stac == 0) - (stac == 1) - 2 * (stac == 4) - 2 * (stac == 2) - (stac == 3)
S = length(brts) - (stac %% 2 == 1) - 2 * (stac %% 2 == 0)
S2 = S - (stac == 1) - (stac == 3) - (stac == 4)
loglik = -lgamma(S2 + missnumspec + 1) + lgamma(S2 + 1) + lgamma(missnumspec + 1)
if(min(pars1) < 0)
{
   loglik = -Inf
} else {
   if((ddep == 1 | ddep == 11) & ceiling(K) < (S + missnumspec))
   {
      loglik = -Inf
   } else {
      if(lac == Inf & mu != Inf & missnumspec == 0)
      {
        loglik = DAISIE_loglik_high_lambda(pars1,-brts,stac)        
      } else {
        if(ddep == 1 | ddep == 11)
        {
           lx = min(max(1 + missnumspec,1 + ceiling(K)),round(pars2[1]) + missnumspec)
        } else {
           lx = roundn(pars2[1]) + missnumspec
        }
        if(loglik > -Inf)
        { 
           # in all cases we integrate from the origin of the island to the first branching point (stac > 1) or to the present (stac <= 1)
           probs = rep(0,2 * lx + 1)
           probs[1] = 1
           k1 = 0
           y = ode(probs,brts[1:2],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
           probs = y[2,2:(2 * lx + 2)]
           cp = checkprobs(lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]      
           if(stac == 0)
           # for stac = 0, the integration is from the origin of the island until the present
           # and we evaluate the probability of no clade being present and no immigrant species,
           # but there can be missing species
           {     
              loglik = loglik + log(probs[1 + missnumspec])
           } else {
             if(stac == 1)
             # for stac = 1, the integration is from the maximum colonization time (usually the
             # island age + tiny time unit) until the present, where we set all probabilities where
             # the immigrant is already present to 0
             # and we evaluate the probability of the immigrant species being present,
             # but there can be missing species           
             {         
                probs[(lx + 1):(2 * lx)] = 0
                y = ode(probs,brts[2:3],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
                probs = y[2,2:(2 * lx + 2)]
                cp = checkprobs(lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]               
                loglik = loglik + log(probs[lx + 1 + missnumspec])
             } else {
             # for stac > 1, integration is then from the colonization event until the first branching time (stac = 2 and 3) or the present (stac = 4). We add a set of equations for Q_M,n, the probability that the process is compatible with the data, and speciation has not happened; during this time immigration is not allowed because it would alter the colonization time. After speciation, colonization is allowed again (re-immigration)
             # all probabilities of states with the immigrant present are set to zero and all probabilities of states with endemics present are transported to the state with the colonist present waiting for speciation to happen. We also multiply by the (possibly diversity-dependent) immigration rate
                gamvec = divdepvec(gam,K,lx,k1,ddep * (ddep == 11 | ddep == 21))
                probs[(2 * lx + 1):(3 * lx)] = gamvec[1:lx] * probs[1:lx]
                probs[1:(2 * lx)] = 0        
                k1 = 1
                y = ode(probs,c(brts[2:3]),DAISIE_loglik_rhs2,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
                probs = y[2,2:(3 * lx + 1)]
                cp = checkprobs2(lx,loglik,probs); loglik = cp[[1]]; probs = cp[[2]]
                if(stac == 4)
                # if stac = 4, we're done and we take an element from Q_M,n
               	{
                 	loglik = loglik + log(probs[2 * lx + 1 + missnumspec])
                } else {         
                # for stac = 2 and 3, at the first branching point all probabilities of states Q_M,n are transferred to probabilities where only endemics are present. Then go through the branching points.
                  S1 = length(brts) - 1
                	if(S1 >= 3)
                  {
                     lacvec = divdepvec(lac,K,lx,k1,ddep)
                     probs[1:lx] = lacvec[1:lx] * (probs[1:lx] + probs[(2 * lx + 1):(3 * lx)])
                     probs[(lx + 1):(2 * lx)] = lacvec[2:(lx + 1)] * probs[(lx + 1):(2 * lx)]
                     probs = probs[-c((2 * lx + 2):(3 * lx))]
                     probs[2 * lx + 1] = 0
                     for(k in 3:S1)
             	       {
              	        k1 = k - 1
             		        y = ode(probs,brts[k:(k+1)],DAISIE_loglik_rhs,c(pars1,k1,ddep),rtol = reltol,atol = abstol,method = methode)
                        probs = y[2,2:(2 * lx + 2)]
             		        if(k < S1)
              	        {
                           # speciation event      
                           lacvec = divdepvec(lac,K,lx,k1,ddep)
          		             probs[1:(2 * lx)] = c(lacvec[1:lx],lacvec[2:(lx + 1)]) * probs[1:(2 * lx)]
             		        }
                     }            
                  }
                  # we evaluate the probability of the phylogeny with any missing species at the present without (stac = 2) or with (stac = 3) the immigrant species; there can be no missing species for stac = 4
                 	loglik = loglik + log(probs[(stac == 3) * lx + 1 + missnumspec])
                }   
              }     
           }           
        }
     }
   }
}

if(pars2[4] == 1)
{
   s1 = sprintf('Status of colonist: %d, Parameters: %f %f %f %f %f',stac,pars1[1],pars1[2],pars1[3],pars1[4],pars1[5])
   s2 = sprintf(', Loglikelihood: %f',loglik)
   cat(s1,s2,"\n",sep = "")
   flush.console()
}

return(as.numeric(loglik))
}

DAISIE_loglik_all = function(pars1,pars2,datalist,methode = "lsodes")
{
# datalist = list of all data: branching times, status of clade, and numnber of missing species
# datalist[[,]][1] = list of branching times (positive, from present to past)
# - max(brts) = age of the island
# - next largest brts = stem age / time of divergence from the mainland
# The interpretation of this depends on stac (see below)
# For stac = 0, this needs to be specified only once.
# For stac = 1, this is the time since divergence from the immigrant's sister on the mainland.
# The immigrant must have immigrated at some point since then.
# For stac = 2 and stac = 3, this is the time since divergence from the mainland.
# The immigrant that established the clade on the island must have immigrated precisely at this point.
# For stac = 3, it must have reimmigrated, but only after the first immigrant had undergone speciation.
# - min(brts) = most recent branching time (only for stac = 2, or stac = 3)
# datalist[[,]][2] = list of status of the clades formed by the immigrant
#  . stac == 0 : immigrant is not present and has not formed an extant clade
# Instead of a list of zeros, here a number must be given with the number of clades having stac = 0
#  . stac == 1 : immigrant is present but has not formed an extant clade
#  . stac == 2 : immigrant is not present but has formed an extant clade
#  . stac == 3 : immigrant is present and has formed an extant clade
#  . stac == 4 : immigrant is present but has not formed an extant clade, and it is known when it immigrated.
# datalist[[,]][3] = list with number of missing species in clades for stac = 2 and stac = 3;
# for stac = 0 and stac = 1, this number equals 0.
# pars1 = model parameters
# - pars1[1] = lac = (initial) cladogenesis rate
# - pars1[2] = mu = extinction rate
# - pars1[3] = K = maximum number of species possible in the clade
# - pars1[4] = gam = (initial) immigration rate
# - pars1[5] = laa = (initial) anagenesis rate
# - pars1[6]...pars1[10] = same as pars1[1]...pars1[5], but for a second type of immigrant
# - pars1[11] = proportion of type 2 immigrants in species pool
# pars2 = model settings
# - pars2[1] = lx = length of ODE variable x
# - pars2[2] = ddep = diversity-dependent model,mode of diversity-dependence
#  . ddep == 0 : no diversity-dependence
#  . ddep == 1 : linear dependence in speciation rate (anagenesis and cladogenesis)
#  . ddep == 11 : linear dependence in speciation rate and immigration rate
#  . ddep == 3 : linear dependence in extinction rate
# - pars2[3] = cond = conditioning
#  . cond == 0 : no conditioning
#  . cond == 1 : conditioning on presence on the island (not used in this single loglikelihood)
# - pars2[4] = parameters and likelihood should be printed (1) or not (0)

cond = pars2[3]
if(length(pars1) == 5)
{
   logp0 = DAISIE_loglik(pars1[1:5],pars2,datalist[[1]]$island_age,0,0,methode)
   if(is.null(datalist[[1]]$not_present))
   {
      loglik = (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) * logp0
      numimm = (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) + length(datalist) - 1
   } else {
      loglik = datalist[[1]]$not_present * logp0
      numimm = datalist[[1]]$not_present + length(datalist) - 1
   }
   logcond = (cond == 1) * log(1 - exp(numimm * logp0))
   for(i in 2:length(datalist))
   {
       datalist[[i]]$type1or2 = 1
   }
} else {
   numimm = datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2 + length(datalist) - 1
   numimm_type2 = length(which(unlist(datalist)[which(names(unlist(datalist)) == "type1or2")] == 2))
   numimm_type1 = length(datalist) - 1 - numimm_type2
   if(is.na(pars1[11]) == FALSE)
   {
       if(pars1[11] < numimm_type2/numimm | pars1[11] > (1 - numimm_type1 /numimm)) { return(-Inf) }
       datalist[[1]]$not_present_type2 = max(0,round(pars1[11] * numimm) - numimm_type2)
       datalist[[1]]$not_present_type1 = numimm - (length(datalist) - 1) - datalist[[1]]$not_present_type2
   }
   logp0_type1 = DAISIE_loglik(pars1[1:5],pars2,datalist[[1]]$island_age,0,0,methode)
   logp0_type2 = DAISIE_loglik(pars1[6:10],pars2,datalist[[1]]$island_age,0,0,methode)
   loglik = datalist[[1]]$not_present_type1 * logp0_type1 + datalist[[1]]$not_present_type2 * logp0_type2
   logcond = (cond == 1) * log(1 - exp((datalist[[1]]$not_present_type1 + numimm_type1) * logp0_type1 + (datalist[[1]]$not_present_type2 + numimm_type2) * logp0_type2))
}
loglik = loglik - logcond

if(length(datalist) > 1)
{
   for(i in 2:length(datalist))
   {
      if(datalist[[i]]$type1or2 == 1)
      {
         pars = pars1[1:5]
      } else {
         pars = pars1[6:10]
      }
      loglik = loglik + DAISIE_loglik(pars,pars2,datalist[[i]]$branching_times,datalist[[i]]$stac,datalist[[i]]$missing_species,methode)
   }
}
return(loglik)
}