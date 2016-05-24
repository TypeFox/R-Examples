bd_loglik_choosepar = function(trparsopt,trparsfix,idparsopt,idparsfix,pars2,brts,missnumspec,methode)
{
   trpars1 = rep(0,4)
   trpars1[idparsopt] = trparsopt
   if(length(idparsfix) != 0)
   {
      trpars1[idparsfix] = trparsfix
   }
   if(max(trpars1) > 1 | min(trpars1) < 0)
   {
      loglik = -Inf
   } else {
      pars1 = trpars1/(1 - trpars1)
      loglik = bd_loglik(pars1 = pars1,pars2 = pars2[1:6],brts = brts,missnumspec = missnumspec, methode = methode)
      if(is.nan(loglik) || is.na(loglik))
      {
         cat("There are parameter values used which cause numerical problems.\n")
         loglik = -Inf
      }
   }
   return(loglik)
}