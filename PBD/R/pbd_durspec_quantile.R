pbd_durspec_quantile = function(pars,p)
{
   expdurspec = pbd_durspec_mean(pars)
   if(expdurspec < 1E-7)
   {
      q = 0
   } else {
      found = 0
      uptau = 100 * expdurspec
      while(found == 0)
      {
          if(pbd_durspec_cumdensity(pars,uptau) > p)
          {
              found = 1
          } else {
              uptau = 10*uptau
          }
      }
      q = uniroot(function(x) pbd_durspec_cumdensity(pars,x) - p,c(0,uptau))$root
   }
   return(q)
}