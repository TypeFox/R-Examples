dd_loglik_choosepar = function(trparsopt,trparsfix,idparsopt,idparsfix,pars2,brts,missnumspec,methode)
{
   trpars1 = rep(0,3)
   if(pars2[2] == 5)
   {   
      trpars1 = rep(0,4)
   }
   trpars1[idparsopt] = trparsopt
   if(length(idparsfix) != 0)
   {
      trpars1[idparsfix] = trparsfix
   }
   if(max(trpars1) > 1 || max(trpars1[1:2]) == 1 || min(trpars1[1:3]) < 0 || trpars1[1] <= trpars1[2])
   {
      loglik = -Inf
   } else {
      pars1 = trpars1/(1 - trpars1)
      loglik = dd_loglik(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec, methode = methode)
      if(is.nan(loglik) || is.na(loglik) || loglik == Inf)
      {
         cat("There are parameter values used which cause numerical problems.\n")
         loglik = -Inf
      }
   }
   return(loglik)
}