DAMOCLES_loglik_choosepar = function(trparsopt,trparsfix,idparsopt,idparsfix,phy,pa1,pa2,pchoice,edgeTList)
{
   pa = cbind(pa1,pa2)
   trpars1 = rep(0,3)
   trpars1[idparsopt] = trparsopt
   if(length(idparsfix) != 0) { trpars1[idparsfix] = trparsfix }
   if(max(trpars1) > 1 || min(trpars1) < 0)
   {
      loglik = -Inf
   } else {
      pars1 = trpars1/(1 - trpars1)
      loglik = DAMOCLES_loglik(phy,pa,pars1,pchoice,edgeTList)
   }
   return(loglik)
}

DAMOCLES_simplex = function(trparsopt,trparsfix,idparsopt,idparsfix,phy,pa1,pa2,pars2,edgeTList)
{
  pchoice = pars2[5]
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
     fv[i] = -DAMOCLES_loglik_choosepar(v[,i],trparsfix,idparsopt,idparsfix,phy,pa1,pa2,pchoice,edgeTList)
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
  reltolx = pars2[1]
  reltolf = pars2[2]
  abstolx = pars2[3]
  maxiter = pars2[4]
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
     fxr = -DAMOCLES_loglik_choosepar(xr,trparsfix,idparsopt,idparsfix,phy,pa1,pa2,pchoice,edgeTList)
   
     if(fxr < fv[1])
     {
         ## Calculate expansion point
         xe = (1 + rh * ch) * xbar - rh * ch * v[,numpar + 1]
         fxe = -DAMOCLES_loglik_choosepar(xe,trparsfix,idparsopt,idparsfix,phy,pa1,pa2,pchoice,edgeTList)
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
                fxco = -DAMOCLES_loglik_choosepar(xco,trparsfix,idparsopt,idparsfix,phy,pa1,pa2,pchoice,edgeTList)
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
                fxci = -DAMOCLES_loglik_choosepar(xci,trparsfix,idparsopt,idparsfix,phy,pa1,pa2,pchoice,edgeTList)
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
                     fv[j] = -DAMOCLES_loglik_choosepar(v[,j],trparsfix,idparsopt,idparsfix,phy,pa1,pa2,pchoice,edgeTList)
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
  out = list(par = v[,1], fvalues = -fv[1], conv = as.numeric(itercount > maxiter))
  invisible(out)
}

DAMOCLES_ML = function(
   phy = rcoal(10),
   pa = matrix(c(phy$tip.label,sample(c(0,1),Ntip(phy),replace = T)),nrow = Ntip(phy),ncol = 2),
   initparsopt = c(0.1,0.1),
   idparsopt = 1:length(initparsopt),
   parsfix = 0,
   idparsfix = (1:3)[-idparsopt],
   pars2 = c(1E-3,1E-4,1E-5,1000),
   pchoice = 0)
{
  if(is.matrix(pa) == 0){pa = matrix(c(phy$tip.label,pa),nrow = length(pa),ncol = 2)}
  options(warn=-1)
  idpars = sort(c(idparsopt,idparsfix))
  if(sum(idpars == (1:3)) != 3)
  {
    cat("The parameters to be optimized and fixed are incoherent.")
  } else {
    namepars = c("mu","gamma_0","gamma_1")
    if(length(namepars[idparsopt]) == 0) { optstr = "nothing" } else { optstr = namepars[idparsopt] }
    cat("You are optimizing",optstr,"\n")
    if(length(namepars[idparsfix]) == 0) { fixstr = "nothing" } else { fixstr = namepars[idparsfix] }
    cat("You are fixing",fixstr,"\n")
    trparsopt = initparsopt/(1 + initparsopt)
    trparsfix = parsfix/(1 + parsfix)
    trparsfix[parsfix == Inf] = 1
    flush.console()
    pars2[5] = pchoice
    edgeTList = compute_edgeTList(phy)
    out = DAMOCLES_simplex(trparsopt = trparsopt,trparsfix = trparsfix,idparsopt = idparsopt,idparsfix = idparsfix,phy = phy,pa1 = pa[,1],pa2 = pa[,2],pars2,edgeTList)
    if(out$conv > 0)
    {
      cat("Optimization has not converged. Try again with different starting values.\n")
    } else {
      MLtrpars = unlist(out$par)
      MLpars = MLtrpars/(1 - MLtrpars)
      out$par = list(MLpars)
      MLpars1 = rep(0,3)
      MLpars1[idparsopt] = MLpars
      if(length(idparsfix) != 0) { MLpars1[idparsfix] = parsfix }
      ML = as.numeric(unlist(out$fvalues))
      out2 = data.frame(mu = MLpars1[1], gamma_0 = MLpars1[2], gamma_1 = MLpars1[3], loglik = ML, df = length(initparsopt), conv = unlist(out$conv))
      s1 = sprintf('Maximum likelihood parameter estimates: mu: %f, gamma_0: %f, gamma_1: %f',MLpars1[1],MLpars1[2],MLpars1[3])
      s2 = sprintf('Maximum loglikelihood: %f',out$fvalues)
      cat("\n",s1,"\n",s2,"\n")
      invisible(out2)
    }
  }  
}

