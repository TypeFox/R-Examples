bisection.update <- function(ModelMx, ObsTbl, tolerance)
{
   obs.s <- sum(ObsTbl);
   b <- suff.stat(ModelMx, ObsTbl/obs.s);
   gamma.one <- 1/(sum(b));
   gamma.two <- min(1/b);
   gamma.mid <- (gamma.one+gamma.two)/2;

   F.gamma.mid <- ipf.gamma(ModelMx, ObsTbl, gamma.mid, tolerance, "probabilities")
   p.mid <- F.gamma.mid$fitted.values;
   
  
   while(abs(sum(p.mid)/obs.s -1)> tolerance)
   {
      F.gamma.two <- ipf.gamma(ModelMx, ObsTbl, gamma.two, tolerance, "probabilities")
      p.two <- F.gamma.two$fitted.values;

      if( sign(sum(p.mid)/obs.s -1) == sign(sum(p.two)/obs.s -1))
      {
         gamma.two <- gamma.mid;
      }   else { gamma.one <- gamma.mid };

      gamma.mid <- (gamma.one + gamma.two)/2;
      F.gamma.mid <- ipf.gamma(ModelMx, ObsTbl, gamma.mid, tolerance, "probabilities")
      p.mid <- F.gamma.mid$fitted.values;
        
    }
    bisection.result <- list(gamma.tilde =  gamma.mid,
                             model.tilde =  F.gamma.mid); 
    return(bisection.result);
}