dd_SR_loglik_choosepar = function(trparsopt,trparsfix,idparsopt,idparsfix,idparsnoshift,pars2,brts,missnumspec,methode)
{
   trpars1 = rep(0,7)
   trpars1[idparsopt] = trparsopt
   if(length(idparsfix) != 0)
   {
      trpars1[idparsfix] = trparsfix
   }
   if(length(idparsnoshift) != 0)
   {
      trpars1[idparsnoshift] = trpars1[idparsnoshift - 3]
   }
   brts = -sort(abs(as.numeric(brts)),decreasing = TRUE)
   pars1 = trpars1/(1 - trpars1)
   if(max(trpars1) > 1 || min(trpars1) < 0 || trpars1[1] <= trpars1[2] || trpars1[4] <= trpars1[5] || -pars1[7] <= min(brts))
   {
      loglik = -Inf
   } else {
      loglik = dd_SR_loglik(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec, methode = methode)
      if(is.nan(loglik) || is.na(loglik))
      {
          cat("There are parameter values used which cause numerical problems.\n")
          loglik = -Inf
      }
   }
   return(loglik)
}