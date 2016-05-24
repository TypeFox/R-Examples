DAISIE_loglik_high_lambda = function(pars1,brts,stac)
{
   lbrts = length(brts)
   if(brts[lbrts] == 0)
   {
       brts = brts[-lbrts]
       lbrts = length(brts)
   }
   N = lbrts - 1
   mu = pars1[2]
   gam = pars1[4]
   brtsdiff = brts - c(brts[2:(N+1)],0)   
   if(stac == 0)
   {
      out = -gam * brts[1]
   }
   if(stac == 2)
   {
      out = -gam * brtsdiff[1] +
        log(gam) +
        log(N) +
        (N - 1) * log(mu) +
        lgamma(N) +
        - (N - 1) * log(N - 1) +
        - mu/(N - 1) * sum((1:N)*(0:(N-1)) * brtsdiff[2:(N+1)])
   }
   if(stac == 1 | stac == 3 | stac == 4)
   {
      out = -Inf
   }   
   return(out)
}
