DAISIE_loglik_all_choosepar = function(trparsopt,trparsfix,idparsopt,idparsfix,idparsnoshift,idparseq,pars2,datalist,methode)
{
   if(sum(idparsnoshift == (6:10)) != 5)
   {
       trpars1 = rep(0,10)
   } else {
       trpars1 = rep(0,5)
   }
   trpars1[idparsopt] = trparsopt
   if(length(idparsfix) != 0)
   {
      trpars1[idparsfix] = trparsfix
   }
   if(sum(idparsnoshift == (6:10)) != 5)
   {
      trpars1[idparsnoshift] = trpars1[idparsnoshift - 5]
   }
   if(max(trpars1) > 1 | min(trpars1) < 0)
   {
      loglik = -Inf
   } else {
      pars1 = trpars1/(1 - trpars1)
      if(pars2[5] > 0)
      {
         pars1 = DAISIE_eq(datalist,pars1,pars2)
         if(sum(idparsnoshift == (6:10)) != 5)
         {
             pars1[idparsnoshift] = pars1[idparsnoshift - 5]
         }
      }
      if(min(pars1) < 0)
      {
         loglik = -Inf
      } else {
         loglik = DAISIE_loglik_all(pars1,pars2,datalist,methode)
      }
      if(is.nan(loglik) || is.na(loglik))
      {
         cat("There are parameter values used which cause numerical problems.\n")
         loglik = -Inf
      }
   }
   return(loglik)
}

DAISIE_simplex = function(trparsopt,idparsopt,trparsfix,idparsfix,idparsnoshift,idparseq,pars2,datalist,methode)
{
numpar = length(trparsopt)
## Setting up initial simplex
v = t(matrix(rep(trparsopt,each = numpar + 1),nrow = numpar + 1))
for(i in 1:numpar)
{
    parsoptff = 1.05 * trparsopt[i]/(1 - trparsopt[i])
    trparsoptff = parsoptff/(1 + parsoptff)
    fac = trparsoptff/trparsopt[i]
    if(v[i,i + 1] == 0)
    {
       v[i,i + 1] = 0.00025
    } else {
       v[i,i + 1] = v[i,i + 1] * min(1.05,fac)
    }
}

fv = rep(0,numpar + 1)
for(i in 1:(numpar + 1))
{
   fv[i] = -DAISIE_loglik_all_choosepar(v[,i],trparsfix,idparsopt,idparsfix,idparsnoshift,idparseq,pars2,datalist,methode)
}

how = "initial"
itercount = 1
string = itercount
for(i in 1:numpar)
{
   string = paste(string, v[i,1]/(1 - v[i,1]), sep = " ")
}
string = paste(string, -fv[1], how, "\n", sep = " ")
cat(string)
flush.console()

tmp = order(fv)
if(numpar == 1)
{
   v = matrix(v[tmp],nrow = 1,ncol = 2)
} else {
   v = v[,tmp]
}
fv = fv[tmp]

## Iterate until stopping criterion is reached
reltolx = pars2[6]
reltolf = pars2[7]
abstolx = pars2[8]
maxiter = pars2[9]
rh = 1
ch = 2
ps = 0.5
si = 0.5

v2 = t(matrix(rep(v[,1],each = numpar + 1),nrow = numpar + 1))

while(itercount <= maxiter & ( ( is.nan(max(abs(fv - fv[1]))) | (max(abs(fv - fv[1])) - reltolf * abs(fv[1]) > 0) ) + ( (max(abs(v - v2) - reltolx * abs(v2)) > 0) | (max(abs(v - v2)) - abstolx > 0) ) ) )
{
   ## Calculate reflection point

   if(numpar == 1)
   {
       xbar = v[1]
   } else {
       xbar = rowSums(v[,1:numpar])/numpar
   }
   xr = (1 + rh) * xbar - rh * v[,numpar + 1]
   fxr = -DAISIE_loglik_all_choosepar(xr,trparsfix,idparsopt,idparsfix,idparsnoshift,idparseq,pars2,datalist,methode)

   if(fxr < fv[1])
   {
       ## Calculate expansion point
       xe = (1 + rh * ch) * xbar - rh * ch * v[,numpar + 1]
       fxe = -DAISIE_loglik_all_choosepar(xe,trparsfix,idparsopt,idparsfix,idparsnoshift,idparseq,pars2,datalist,methode)
       if(fxe < fxr)
       {
           v[,numpar + 1] = xe
           fv[numpar + 1] = fxe
           how = "expand"
       } else {
           v[,numpar + 1] = xr
           fv[numpar + 1] = fxr
           how = "reflect"
       }
   } else {
       if(fxr < fv[numpar])
       {
           v[,numpar + 1] = xr
           fv[numpar + 1] = fxr
           how = "reflect"
       } else {
           if(fxr < fv[numpar + 1])
           {
              ## Calculate outside contraction point
              xco = (1 + ps * rh) * xbar - ps * rh * v[,numpar + 1]
              fxco = -DAISIE_loglik_all_choosepar(xco,trparsfix,idparsopt,idparsfix,idparsnoshift,idparseq,pars2,datalist,methode)
              if(fxco <= fxr)
              {
                 v[,numpar + 1] = xco
                 fv[numpar + 1] = fxco
                 how = "contract outside"
              } else {
                 how = "shrink"
              }
           } else {
              ## Calculate inside contraction point
              xci = (1 - ps) * xbar + ps * v[,numpar + 1]
              fxci = -DAISIE_loglik_all_choosepar(xci,trparsfix,idparsopt,idparsfix,idparsnoshift,idparseq,pars2,datalist,methode)
              if(fxci < fv[numpar + 1])
              {
                 v[,numpar + 1] = xci
                 fv[numpar + 1] = fxci
                 how = "contract inside"
              } else {
                 how = "shrink"
              }
           }
           if(how == "shrink")
           {
               for(j in 2:(numpar + 1))
               {

                   v[,j] = v[,1] + si * (v[,j] - v[,1])
                   fv[j] = -DAISIE_loglik_all_choosepar(v[,j],trparsfix,idparsopt,idparsfix,idparsnoshift,idparseq,pars2,datalist,methode)
               }
           }
       }
   }
   tmp = order(fv)
   if(numpar == 1)
   {
      v = matrix(v[tmp],nrow = 1,ncol = 2)
   } else {
      v = v[,tmp]
   }
   fv = fv[tmp]
   itercount = itercount + 1
   string = itercount;
   for(i in 1:numpar)
   {
       string = paste(string, v[i,1]/(1 - v[i,1]), sep = " ")
   }
   string = paste(string, -fv[1], how, "\n", sep = " ")
   cat(string)
   flush.console()
   v2 = t(matrix(rep(v[,1],each = numpar + 1),nrow = numpar + 1))
}
if(itercount < maxiter)
{
   cat("Optimization has terminated successfully.","\n")
} else {
   cat("Maximum number of iterations has been exceeded.","\n")
}
out = list(par = v[,1], fvalues = -fv[1], conv = -as.numeric(itercount > maxiter))
invisible(out)
}

DAISIE_ML = function(datalist,initparsopt,idparsopt,parsfix,idparsfix,idparsnoshift = 6:10, res = 100, ddmodel = 0, cond = 0, eqmodel = 0, x_E = 0.95, x_I = 0.98, tol = c(1E-4, 1E-5, 1E-7), maxiter = 1000 * round((1.25)^length(idparsopt)), methode = "lsodes")
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
# datalist[[,]][3] = list with number of missing species in clades for stac = 2 and stac = 3;
# for stac = 0 and stac = 1, this number equals 0.
# initparsopt, parsfix = optimized and fixed model parameters
# - pars1[1] = lac = (initial) cladogenesis rate
# - pars1[2] = mu = extinction rate
# - pars1[3] = K = maximum number of species possible in the clade
# - pars1[4] = gam = (initial) immigration rate
# - pars1[5] = laa = (initial) anagenesis rate
# - pars1[6]...pars1[10] = same as pars1[1]...pars1[5], but for a second type of immigrant
# - pars1[11] = proportion of type 2 immigrants in species pool
# idparsopt, idparsfix = ids of optimized and fixed model parameters
# - res = pars2[1] = lx = length of ODE variable x
# - ddmodel = pars2[2] = diversity-dependent model,mode of diversity-dependence
#  . ddmodel == 0 : no diversity-dependence
#  . ddmodel == 1 : linear dependence in speciation rate (anagenesis and cladogenesis)
#  . ddmodel == 11 : linear dependence in speciation rate and immigration rate
#  . ddmodel == 3 : linear dependence in extinction rate
# - cond = conditioning; currently only cond = 0 is possible
#  . cond == 0 : no conditioning
#  . cond == 1 : conditioning on presence on the island
# - eqmodel = equilibrium model
#  . eqmodel = 0 : no equilibrium is assumed
#  . eqmodel = 1 : equilibrium is assumed on deterministic equation for total number of species
#  . eqmodel = 2 : equilibrium is assumed on total number of species using deterministic equation for endemics and immigrants
#  . eqmodel = 3 : equilibrium is assumed on endemics using deterministic equation for endemics and immigrants
#  . eqmodel = 4 : equilibrium is assumed on immigrants using deterministic equation for endemics and immigrants
#  . eqmodel = 5 : equilibrium is assumed on endemics and immigrants using deterministic equation for endemics and immigrants

options(warn=-1)
idparseq = c()
if(eqmodel == 1 | eqmodel == 3 | eqmodel == 13)
{
   idparseq = 2
}
if(eqmodel == 2 | eqmodel == 4)
{
   idparseq = 4
}
if(eqmodel == 5 | eqmodel == 15)
{
   idparseq = c(2,4)
}

idpars = sort(c(idparsopt,idparsfix,idparsnoshift,idparseq))
#print(idpars)
missnumspec = unlist(lapply(datalist,function(list) {list$missing_species}))
if(sum(missnumspec) > (res - 1))
{
   cat("The number of missing species is too large relative to the resolution of the ODE.\n")
   out2 = data.frame(lambda_c = -1, mu = -1,K = -1, gamma = -1, lambda_a = -1, loglik = -1, df = -1, conv = -1)
} else {
  if((sum(idpars == (1:10)) != 10) || (length(initparsopt) != length(idparsopt)) || (length(parsfix) != length(idparsfix)))
  {
     cat("The parameters to be optimized and/or fixed are incoherent.\n")
     out2 = data.frame(lambda_c = -1, mu = -1,K = -1, gamma = -1, lambda_a = -1, loglik = -1, df = -1, conv = -1)
  } else {
    if(length(idparsopt) > 11)
    {
       cat("The number of parameters to be optimized is too high.\n")
       out2 = data.frame(lambda_c = -1, mu = -1,K = -1, gamma = -1, lambda_a = -1, loglik = -1, df = -1, conv = -1)
    } else {
      namepars = c("lambda_c","mu","K","gamma","lambda_a","lambda_c2","mu2","K2","gamma2","lambda_a2","prop_type2")
      if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
      cat("You are optimizing",optstr,"\n")
      if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
      cat("You are fixing",fixstr,"\n")
      if(sum(idparsnoshift == (6:10)) != 5)
      {
         noshiftstring = namepars[idparsnoshift]
         cat("You are not shifting",noshiftstring,"\n")
      }
      if(length(idparseq) == 0)
      {
         #cat("You are not assuming equilibrium\n")
      } else {
         if(ddmodel == 3)
         {
            cat("Equilibrium optimization is not implemented for ddmodel = 3\n")
         } else {
            cat("You are assuming equilibrium. Extinction and/or immigration will be considered a function of the other parameters, the species pool size, the number of endemics, and/or the number of non-endemics\n")
         }
      }
      cat("Calculating the likelihood for the initial parameters.","\n")
      flush.console()
      trparsopt = initparsopt/(1 + initparsopt)
      trparsopt[which(initparsopt == Inf)] = 1
      trparsfix = parsfix/(1 + parsfix)
      trparsfix[which(parsfix == Inf)] = 1
      pars2 = c(res,ddmodel,cond,0,eqmodel,tol,maxiter,x_E,x_I)
      initloglik = DAISIE_loglik_all_choosepar(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,idparsnoshift = idparsnoshift,idparseq = idparseq, pars2 = pars2,datalist = datalist,methode)
      cat("The loglikelihood for the initial parameter values is",initloglik,"\n")
      if(initloglik == -Inf)
      {
         cat("The initial parameter values have a likelihood that is equal to 0 or below machine precision. Try again with different initial values.\n")
         out2 = data.frame(lambda_c = -1, mu = -1,K = -1, gamma = -1, lambda_a = -1, loglik = -1, df = -1, conv = -1)
      } else {
        cat("Optimizing the likelihood - this may take a while.","\n")
        flush.console()
        out = DAISIE_simplex(trparsopt,idparsopt,trparsfix,idparsfix,idparsnoshift,idparseq,pars2,datalist,methode)
        if(out$conv > 0)
        {
           cat("Optimization has not converged. Try again with different initial values.\n")
           out2 = data.frame(lambda_c = -1, mu = -1,K = -1, gamma = -1, lambda_a = -1, loglik = -1, df = -1, conv = -1)
        } else {
          MLtrpars = as.numeric(unlist(out$par))
          MLpars = MLtrpars/(1-MLtrpars)
          ML = as.numeric(unlist(out$fvalues))
          if(sum(idparsnoshift == (6:10)) != 5)
          {
              MLpars1 = rep(0,10)
          } else {
              MLpars1 = rep(0,5)
          }
          MLpars1[idparsopt] = MLpars
          if(length(idparsfix) != 0) { MLpars1[idparsfix] = parsfix }
          if(eqmodel > 0)
          {
              MLpars1 = DAISIE_eq(datalist,MLpars1,pars2)
          }
          if(MLpars1[3] > 10^7){ MLpars1[3] = Inf }
          if(sum(idparsnoshift == (6:10)) != 5)
          {
              if(length(idparsnoshift) != 0) { MLpars1[idparsnoshift] = MLpars1[idparsnoshift - 5] }
              if(MLpars1[8] > 10^7){ MLpars1[8] = Inf }
              out2 = data.frame(lambda_c = MLpars1[1], mu = MLpars1[2], K = MLpars1[3], gamma = MLpars1[4], lambda_a = MLpars1[5], lambda_c2 = MLpars1[6], mu2 = MLpars1[7], K2 = MLpars1[8], gamma2 = MLpars1[9], lambda_a2 = MLpars1[10], prop_type2 = MLpars1[11], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
              s1 = sprintf('Maximum likelihood parameter estimates: lambda_c: %f, mu: %f, K: %f, gamma: %f, lambda_a: %f, lambda_c2: %f, mu2: %f, K2: %f, gamma2: %f, lambda_a2: %f, prop_type2: %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5],MLpars1[6],MLpars1[7],MLpars1[8],MLpars1[9],MLpars1[10],MLpars1[11])
          } else {
              out2 = data.frame(lambda_c = MLpars1[1], mu = MLpars1[2], K = MLpars1[3], gamma = MLpars1[4], lambda_a = MLpars1[5], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
              s1 = sprintf('Maximum likelihood parameter estimates: lambda_c: %f, mu: %f, K: %f, gamma: %f, lambda_a: %f',MLpars1[1],MLpars1[2],MLpars1[3],MLpars1[4],MLpars1[5])
          }
          s2 = sprintf('Maximum loglikelihood: %f',ML)
          cat("\n",s1,"\n",s2,"\n")
          if(eqmodel > 0)
          {
              M = calcMN(datalist,MLpars1)
              ExpEIN = DAISIE_ExpEIN(datalist[[1]]$island_age,MLpars1,M)
              cat("The expected number of endemics, non-endemics, and the total at these parameters is: ", ExpEIN[[1]], ExpEIN[[2]],ExpEIN[[3]])
          }
        }
      }
    }
  }
}
invisible(out2)
}