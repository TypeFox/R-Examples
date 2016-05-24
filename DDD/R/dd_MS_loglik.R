# brtsM = branching times of main clade M (positive, from present to past)
# brtsS = branching times of subclade S (positive, from present to past)
# - max(brts) = crown age
# - min(brts) = most recent branching time
# - pars1[1] = laM = (initial) speciation rate of main clade
# - pars1[2] = muM = extinction rate of main clade
# - pars1[3] = K = carrying capacity
# - pars1[4] = laS = (initial) speciation rate of subclade
# - pars1[5] = muS = extinction rate of subclade
# - pars1[6] = tinn = time of key innovation
# - pars2[1] = lx = length of ODE variable x
# - pars2[2] = ddep = diversity-dependent model, mode of diversity-dependence
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
# - pars2[3] = cond = conditioning
#  . cond == 0 : no conditioning
#  . cond == 1 : conditioning on non-extinction of the phylogeny
# - pars2[4] = tsplit = time of split of innovative branch
# - pars2[5] = printing of parameters and likelihood (1) or not (0)
# - pars2[6] = likelihood is for a tree with crown age (2) or stem age (1)
# missnumspec = number of missing species in main clade M and subclade S
# methode = the method used in the numerical solving of the set of the ode's

dd_MS_loglik = function(pars1,pars2,brtsM,brtsS,missnumspec,methode = 'analytical')
{
   if(methode == 'analytical')
   {
       out = dd_MS_loglik2(pars1,pars2,brtsM,brtsS,missnumspec)
   } else {
       out = dd_MS_loglik1(pars1,pars2,brtsM,brtsS,missnumspec,methode = methode)
   }
   return(out)
}

dd_MS_loglik1 = function(pars1,pars2,brtsM,brtsS,missnumspec,methode = 'analytical')
{
if(length(pars2) == 4)
{
   pars2[5] = 0
   pars2[6] = 2
}
brtsM = -sort(abs(brtsM),decreasing = TRUE)
maxbrtsS = 0
if(!is.null(brtsS))
{
   brtsS = -sort(abs(brtsS),decreasing = TRUE)
   maxbrtsS = max(brtsS)
}
tinn = 0
if(!is.na(pars1[6]))
{
   tinn = -abs(pars1[6])
}
if(min(pars1[1:5]) < 0 | tinn <= min(brtsM) | tinn > maxbrtsS)
{
   loglik = -Inf         
} else {
   if(((pars1[2] == 0 | pars1[4] == 0) & pars2[2] == 2) | ((pars1[1] == 0 | pars1[3] == 0) & pars2[2] == 4) | pars1[1] <= pars1[2] | pars1[4] <= pars1[5])
   { 
      cat("These parameter values cannot satisfy lambda(N) = mu(N) for some finite N.\n")
      loglik = -Inf
   } else {
      abstol = 1e-16
      reltol = 1e-14
      m = missnumspec
      summ = sum(m)
      # order branching times
      soc = pars2[6]
      S1 = length(brtsM) + (soc - 2)
      S2 = length(brtsS) + 1
      S = S1 + S2 
      if(is.na(pars1[6]))
      {
         tinn = NULL
      } 
      # avoid coincidence of branching time and key innovation time
      if(!is.null(tinn))
      {
         if(sum(abs(c(brtsM,brtsS) - tinn) < 1E-14) == 1)
         { 
            tinn = tinn - 1E-8
         }
      }
      brtsM = rbind(brtsM,1)
      if(length(brtsS) > 0)
      {
         brtsS = rbind(brtsS,2)
      }
      tinn1 = 1E-20
      if(length(tinn) > 0)
      {
         tinn = rbind(tinn,3)
         tinn1 = tinn[1]
      }
      brts = abs(cbind(brtsM,brtsS,tinn))
      brtsorder = order(abs(brts[1,]),decreasing = T)
      brts = brts[,brtsorder]
      brts[1,] = -abs(brts[1,])
      if(sum(brts[1,] == 0) == 0)
      { 
         brts = cbind(brts,rbind(0,0))
      }
      laM = pars1[1]
      muM = pars1[2]
      K = pars1[3]
      laS = pars1[4]
      muS = pars1[5]   
      lmax = pars2[1]
      ddep = pars2[2]
      cond = pars2[3]
      tsplit = -abs(pars2[4])
      soc = pars2[6]
      if(ddep != 1.3 & ddep != 2.3)
      {
          cat("This only works for ddmodel = 1.2 or ddmodel == 1.3.\n")
          loglik = -Inf
      } else {
         lx = ceiling(min(max(1 + m[1],1 + K),lmax))
         lx2 = lx * lx  
         n0 = (ddep == 2 | ddep == 4)
         if(ddep == 1.3 & ceiling(K) < S1 + S2 + summ)
         { 
            cat("The carrying capacity is lower than the actual number of species.\n")
            loglik = -Inf
         } else {
            loglik = 0
            # compute likelihood
            #probs = matrix(0,lx,lx)
            #probs[1,1] = 1
            #dim(probs) = c(lx2,1)
            probs = rep(0,lx * lx)
            probs[1] = 1 # change if other species at crown age
            #loglik2 = 0
            #probs2 = rep(0,lx)
            #probs2[1] = 1
            kM = 2
            kS = 0
            nx1 = rep(-1:lx,lx + 2)
            dim(nx1) = c(lx + 2,lx + 2) # row index = number of species in first group 
            nx2 = t(nx1) # column index = number of species in second group
            nxt = nx1 + nx2
            for(i in 2:length(brts[1,]))
            {
               t1 = brts[1,i - 1]
               t2 = brts[1,i]                           
               nxt2 = nxt + kM + kS
               if(ddep == 1.3) 
               { 
                   lavec = pmax(matrix(0,lx + 2,lx + 2),(1 - nxt2/K))
                   muvec = matrix(1,lx + 2,lx + 2)
               } 
               if(ddep == 2.3)
               { 
                   x = K
                   lavec = pmax(matrix(0,lx + 2,lx + 2),(nxt2 + n0)^x)
                   muvec = matrix(1,lx + 2,lx + 2)                   
               }    
               laMvec = laM * lavec
               laSvec = laS * lavec
               muMvec = muM * muvec
               muSvec = muS * muvec
               mm = list()
               mm[[1]] = laMvec[1:lx,2:(lx+1)] * (nx1[1:lx,2:(lx+1)] + 2 * kM)
               mm[[2]] = muMvec[3:(lx+2),2:(lx+1)] * nx1[3:(lx+2),2:(lx+1)]
               mm[[3]] = (laMvec[2:(lx+1),2:(lx+1)] + muMvec[2:(lx+1),2:(lx+1)]) * (nx1[2:(lx+1),2:(lx+1)] + kM)
               mm[[4]] = laSvec[2:(lx+1),1:lx] * (nx2[2:(lx+1),1:lx] + 2 * kS)
               mm[[5]] = muSvec[2:(lx+1),3:(lx+2)] * nx2[2:(lx+1),3:(lx+2)]
               mm[[6]] = (laSvec[2:(lx+1),2:(lx+1)] + muSvec[2:(lx+1),2:(lx+1)]) * (nx2[2:(lx+1),2:(lx+1)] + kS)
               y = ode(probs,c(t1,t2),dd_logliknorm_rhs2,mm,rtol = reltol,atol = abstol,method = methode)
               probs = y[2,2:(lx2 + 1)]
               dim(probs) = c(lx,lx)
               
               #k1 = i + (soc - 2)
               #y = lsoda(probs2,c(t1,t2),dd_loglik_rhs,c(pars1[1:3],k1,ddep),rtol = reltol,atol = abstol)
               #probs2 = y[2,2:(lx+1)]              
               
               if(t2 < 0 & t2 != tinn1)
               {
                  if(brts[2,i] == 1)
                  { 
                     lavec = laMvec
                  }
                  if(brts[2,i] == 2)
                  { 
                     lavec = laSvec
                  }
                  probs = lavec[2:(lx+1),2:(lx+1)] * probs # speciation event
                  sumprobs = sum(probs)
                  if(sumprobs <= 0)
                  { 
                     loglik = -Inf
                     break
                  } else {
                     loglik = loglik + log(sumprobs)
                  }
                  probs = probs/sumprobs
                  
                  #probs2 = flavec(ddep,laM,muM,K,0,lx,k1,n0) * probs2 # speciation event
                  #loglik2 = loglik2 + log(sum(probs2))
                  #probs2 = probs2/sum(probs2)   
               }
               kM = kM + (brts[2,i] == 1) - (brts[2,i] == 3)
               kS = kS + (brts[2,i] >= 2)
            }    
            #loglik2 = loglik2 + log(probs2[1])
            if(length(m) == 1)
            { 
               probstot = 0
               for(i in 0:m)
               {
                  probstot = probstot + probs[1 + i,1 + m - i]  
               }
               loglik = loglik + log(probstot)
            } else {
               loglik = loglik + log(probs[1 + m[1],1 + m[2]])   
            }          
            if(is.nan(loglik) | is.na(loglik))
            {
               loglik = -Inf
            }       
            if(cond == 0 | loglik == -Inf)
            {
               logliknorm = 0
            } else {   
               cat("Cond > 0 is not yet implemented.\n")
               # COMPUTE NORMALIZATION
               tcrown = brts[1]
               tpres = 0
               logliknorm = 0 ###
            }
            if(length(m) > 1)
            {
               Sv = c(S1,S2)
            } else {
               Sv = S
            }
            loglik = loglik - logliknorm - sum(lgamma(Sv + m + 1) - lgamma(Sv + 1) - lgamma(m + 1))           
         }          
      }
   } 
}
if(pars2[5] == 1)
{
    s1 = sprintf('Parameters: %f %f %f %f %f %f, ',pars1[1],pars1[2],pars1[3],pars1[4],pars1[5],pars1[6])
    s2 = sprintf('Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    flush.console()
}
loglik = as.numeric(loglik)
if(is.nan(loglik) | is.na(loglik))
{
    loglik = -Inf
}
return(loglik)
}   


dd_MS_loglik2 = function(pars1,pars2,brtsM,brtsS,missnumspec,methode = 'lsoda')
{
if(length(pars2) == 4)
{
   pars2[5] = 0
   pars2[6] = 2
}
brtsM = -sort(abs(brtsM),decreasing = TRUE)
maxbrtsS = 0
if(!is.null(brtsS))
{
   brtsS = -sort(abs(brtsS),decreasing = TRUE)
   maxbrtsS = max(brtsS)
}
tinn = 0
if(!is.na(pars1[6]))
{
   tinn = -abs(pars1[6])
}
if(min(pars1[1:5]) < 0 | tinn <= min(brtsM) | tinn > maxbrtsS)
{
   loglik = -Inf         
} else {
   if(((pars1[2] == 0 | pars1[4] == 0) & pars2[2] == 2) | ((pars1[1] == 0 | pars1[3] == 0) & pars2[2] == 4) | pars1[1] <= pars1[2] | pars1[4] <= pars1[5])
   { 
      cat("These parameter values cannot satisfy lambda(N) = mu(N) for some finite N.\n")
      loglik = -Inf
   } else {
      abstol = 1e-16
      reltol = 1e-14
      m = missnumspec
      summ = sum(m)
      # order branching times
      soc = pars2[6]
      S1 = length(brtsM) + (soc - 2)
      S2 = length(brtsS) + 1
      S = S1 + S2 
      if(is.na(pars1[6]))
      {
         tinn = NULL
      } 
      # avoid coincidence of branching time and key innovation time
      if(!is.null(tinn))
      {
         if(sum(abs(c(brtsM,brtsS) - tinn) < 1E-14) == 1)
         { 
            tinn = tinn - 1E-8
         }
      }
      brtsM = rbind(brtsM,1)
      if(length(brtsS) > 0)
      {
         brtsS = rbind(brtsS,2)
      }
      tinn1 = 1E-20
      if(length(tinn) > 0)
      {
         tinn = rbind(tinn,3)
         tinn1 = tinn[1]
      }
      brts = abs(cbind(brtsM,brtsS,tinn))
      brtsorder = order(abs(brts[1,]),decreasing = T)
      brts = brts[,brtsorder]
      brts[1,] = -abs(brts[1,])
      if(sum(brts[1,] == 0) == 0)
      { 
         brts = cbind(brts,rbind(0,0))
      }
      laM = pars1[1]
      muM = pars1[2]
      K = pars1[3]
      laS = pars1[4]
      muS = pars1[5]   
      lmax = pars2[1]
      ddep = pars2[2]
      cond = pars2[3]
      tsplit = -abs(pars2[4])
      soc = pars2[6]
      if(ddep != 1.3 & ddep != 2.3)
      {
          cat("This only works for ddmodel = 1.2 or ddmodel == 1.3.\n")
          loglik = -Inf
      } else {
         lx = ceiling(min(max(1 + m[1],1 + K),lmax))
         lx2 = lx * lx  
         n0 = (ddep == 2 | ddep == 4)
         if(ddep == 1.3 & ceiling(K) < S1 + S2 + summ)
         { 
            cat("The carrying capacity is lower than the actual number of species.\n")
            loglik = -Inf
         } else {
            loglik = 0
            # compute likelihood
            probs = rep(0,lx2)
            probs[1] = 1 # change if other species at crown age
            #probs2 = probs
            kM = 2
            kS = 0
            nx1 = rep(-1:lx,lx + 2)
            dim(nx1) = c(lx + 2,lx + 2) # row index = number of species in first group 
            nx2 = t(nx1) # column index = number of species in second group
            nxt = nx1 + nx2
            for(i in 2:length(brts[1,]))
            {
               t1 = brts[1,i - 1]
               t2 = brts[1,i]                           
               nxt2 = nxt + kM + kS
               if(ddep == 1.3) 
               { 
                   lavec = pmax(matrix(0,lx + 2,lx + 2),(1 - nxt2/K))
                   muvec = matrix(1,lx + 2,lx + 2)
               } 
               if(ddep == 2.3)
               { 
                   x = K
                   lavec = pmax(matrix(0,lx + 2,lx + 2),(nxt2 + n0)^x)
                   muvec = matrix(1,lx + 2,lx + 2)                   
               }    
               laMvec = laM * lavec
               laSvec = laS * lavec
               #muMvec = muM * muvec
               #muSvec = muS * muvec
               #mm = list()
               #mm[[1]] = laMvec[1:lx,2:(lx+1)] * (nx1[1:lx,2:(lx+1)] + 2 * kM)
               #mm[[2]] = muMvec[3:(lx+2),2:(lx+1)] * nx1[3:(lx+2),2:(lx+1)]
               #mm[[3]] = (laMvec[2:(lx+1),2:(lx+1)] + muMvec[2:(lx+1),2:(lx+1)]) * (nx1[2:(lx+1),2:(lx+1)] + kM)           
               #mm[[4]] = laSvec[2:(lx+1),1:lx] * (nx2[2:(lx+1),1:lx] + 2 * kS)
               #mm[[5]] = muSvec[2:(lx+1),3:(lx+2)] * nx2[2:(lx+1),3:(lx+2)]
               #mm[[6]] = (laSvec[2:(lx+1),2:(lx+1)] + muSvec[2:(lx+1),2:(lx+1)]) * (nx2[2:(lx+1),2:(lx+1)] + kS)
               #y = ode(probs2,c(t1,t2),dd_logliknorm_rhs2,mm,rtol = reltol,atol = abstol,method = methode)
               #probs2 = y[2,2:(lx2 + 1)]
               #print(as.numeric(probs2[1:10]))
               probs = dd_loglik_M3(pars1,lx,ddep,tt = abs(t2 - t1),p = probs,kM,kS)
               #print(as.numeric(probs[seq(1,10*lx,by = lx)]))
               dim(probs) = c(lx,lx)
               #dim(probs2) = c(lx,lx)
               if(t2 < 0 & t2 != tinn1)
               {
                  if(brts[2,i] == 1)
                  { 
                     lavec = laMvec
                  } else if(brts[2,i] == 2)
                  { 
                     lavec = laSvec
                  }
                  probs = lavec[2:(lx+1),2:(lx+1)] * probs # speciation event
                  #probs2 = lavec[2:(lx+1),2:(lx+1)] * probs2
                  dim(probs) = c(lx2,1)
                  #print(as.numeric(probs2[1:10]))  
                  #print(as.numeric(probs[seq(1,10*lx,by = lx)]))            
                  sumprobs = sum(probs)
                  #sumprobs2 = sum(probs2)
                  if(sumprobs <= 0)
                  { 
                     loglik = -Inf
                     break
                  } else {
                     loglik = loglik + log(sumprobs)
                     #loglik = loglik + log(sumprobs2)
                  }
                  probs = probs/sumprobs   
                  #probs2 = probs2/sumprobs2
                  #print(as.numeric(probs2[1:10]))  
                  #print(as.numeric(probs[seq(1,10*lx,by = lx)]))   
               }
               kM = kM + (brts[2,i] == 1) - (brts[2,i] == 3)
               kS = kS + (brts[2,i] >= 2)
            }    
            dim(probs) = c(lx,lx)
            probs = t(probs)
            if(length(m) == 1)
            { 
               probstot = 0
               for(i in 0:m)
               {
                  probstot = probstot + probs[1 + i,1 + m - i]  
               }
               loglik = loglik + log(probstot)
            } else {
               loglik = loglik + log(probs[1 + m[1],1 + m[2]])   
            }          
            if(is.nan(loglik) | is.na(loglik))
            {
               loglik = -Inf
            }       
            if(cond == 0 | loglik == -Inf)
            {
               logliknorm = 0
            } else {   
               cat("Cond > 0 is not yet implemented.\n")
               # COMPUTE NORMALIZATION
               tcrown = brts[1]
               tpres = 0
               logliknorm = 0 ###
            }
            if(length(m) > 1)
            {
               Sv = c(S1,S2)
            } else {
               Sv = S
            }
            loglik = loglik - logliknorm - sum(lgamma(Sv + m + 1) - lgamma(Sv + 1) - lgamma(m + 1))           
         }          
      }
   } 
}
if(pars2[5] == 1)
{
    s1 = sprintf('Parameters: %f %f %f %f %f %f, ',pars1[1],pars1[2],pars1[3],pars1[4],pars1[5],pars1[6])
    s2 = sprintf('Loglikelihood: %f',loglik)
    cat(s1,s2,"\n",sep = "")
    flush.console()
}
loglik = as.numeric(loglik)
if(is.nan(loglik) | is.na(loglik))
{
    loglik = -Inf
}
return(loglik)
}   
