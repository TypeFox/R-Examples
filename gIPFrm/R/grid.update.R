grid.update <- function(ModelMx, ObsTbl, tolerance)
{
   
   T <- 1;
   #T <- length(ObsTbl);
   m <- 0;
   s.star <- 0;
   p.star <- NULL;
   F.star <- NULL;
   gamma.star <- NULL;
   obs.s <- sum(ObsTbl);
   b <- suff.stat(ModelMx, ObsTbl/obs.s);
           
           gamma.left <- 1/(sum(b));
           gamma.right <- min(1/b);

   
   while(s.star == 0)
   {
        gamma.m <- gamma.left + m*(gamma.right - gamma.left)/T;
        tolerance.T <- max(c(2.220446e-16, tolerance/{T}));

        F.gamma.m <- ipf.gamma(ModelMx, ObsTbl, gamma.m, tolerance.T, "probabilities")
        p.gamma.m <- F.gamma.m$fitted.values; 
        
        if(abs(sum(p.gamma.m)/obs.s -1)<tolerance)
        {
           p.star <- p.gamma.m;
           F.star <- F.gamma.m;
           gamma.star <- gamma.m;
           s.star <- 1
        }
        else
        {
            if(m < T)   
            {
                m <- m + 1;
            }  
            else
            {
               T <- T + 1; 
               m <- 0
            }        
         }
    }

    grid.result <- list(gamma.tilde =  gamma.star,
                             model.tilde =  F.star); 
    return(grid.result);
}