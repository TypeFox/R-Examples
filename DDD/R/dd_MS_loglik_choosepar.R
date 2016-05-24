dd_MS_loglik_choosepar = function(trparsopt,trparsfix,idparsopt,idparsfix,idparsnoshift,pars2,brtsM,brtsS,missnumspec,methode)
{   
   trpars1 = rep(0,6)
   trpars1[idparsopt] = trparsopt
   if(length(idparsfix) != 0)
   {
      trpars1[idparsfix] = trparsfix
   }
   if(length(idparsnoshift) != 0)
   {
      trpars1[idparsnoshift] = trpars1[idparsnoshift - 3]
   }
   brtsM = -sort(abs(as.numeric(brtsM)),decreasing = TRUE)
   brtsS = -sort(abs(as.numeric(brtsS)),decreasing = TRUE)
   pars1 = trpars1/(1 - trpars1)
   if(max(trpars1) > 1 || min(trpars1) < 0 || trpars1[1] <= trpars1[2] || trpars1[4] <= trpars1[5] || -abs(pars1[6]) <= -abs(pars2[4])  || -abs(pars1[6]) >= min(brtsS) )
   {
      loglik = -Inf
   } else {
      loglik = dd_MS_loglik(pars1 = pars1,pars2 = pars2,brtsM = brtsM,brtsS = brtsS,missnumspec = missnumspec, methode = methode)
      if(is.nan(loglik) || is.na(loglik))
      {
         cat("There are parameter values used which cause numerical problems.\n")
         loglik = -Inf
      }
   }
   return(loglik)
}