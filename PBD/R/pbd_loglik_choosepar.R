pbd_loglik_choosepar = function(trparsopt,trparsfix,idparsopt,idparsfix,exteq,parsfunc,pars2,brts,missnumspec)
{
   trpars1 = rep(0,4)  
   trpars1[idparsopt] = trparsopt
   if(length(idparsfix) != 0)
   {
      trpars1[idparsfix] = trparsfix
   }
   if(exteq == 1)
   {
      trpars1[4] = trpars1[2]
   }
   if(max(trpars1) > 1 || min(trpars1) < 0 || max(trpars1[c(1,2,4)]) == 1)
   {
      loglik = -Inf
   } else {
      pars1 = trpars1/(1 - trpars1)
      loglik = pbd_loglik(pars1,parsfunc,pars2,brts,missnumspec)
      if(is.nan(loglik) || is.na(loglik))
      {
         cat("There are parameter values used which cause numerical problems.\n")
         loglik = -Inf
      }
   }
   return(loglik)
}