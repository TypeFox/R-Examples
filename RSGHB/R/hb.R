hb <- function(a, b, d, f, env = parent.frame())
{
    
     env$starttime <- Sys.time()
     
     # sets random seed.	  
     if (env$gSeed == 0) env$gSeed <- ceiling(runif(1) * 1000000)
     	       
     set.seed(env$gSeed, kind = "default", normal.kind = "default")
     
     p <- env$likelihood(f, b, env)
     
     # the burn-in iterations
     for (r in 1:env$gNCREP)
     {

          if (env$gNIV > 0) {
          
               out    <- nextB(a, b, d, p, f, env)
               b      <- out[[1]]
               env$acceptanceRatePerc <- out[[2]]
               p      <- out[[3]]
               
               a <- nextA(b, d, env)
               
               # this constrains the means of the random parameters to be user-specified
               # this is necessary for such things as error components logit.
               if (!is.null(env$fixedA))
               {
                    for (rp in 1:env$gNIV)
                    {
                         if (!is.na(env$fixedA[rp])) a[rp] <- env$fixedA[rp]
                    }
               }     
                    
               if (env$gFULLCV) d <- nextD(a, b, env)
               if (!env$gFULLCV) d <- nextDind(a, b, env)
               
               if (!is.null(env$fixedD))
               {
                    for (rp in 1:env$gNIV)
                    {
                         if (!is.na(env$fixedD[rp])) d[rp, rp] <- env$fixedD[rp]
                    }
               } 
               
               
          }

          # drawing a new set of fixed coefficients
          if (env$gFIV > 0) {
               
               out <- nextF(p, f, b, env)
               
               if (sum(out[[1]] == f) != env$gFIV) env$acceptanceRateF <- env$acceptanceRateF + 1
               
               # targeting an acceptance rate of 0.25
               if ((r %% 100) == 0)
               {
                    env$acceptanceRateFPerc <- env$acceptanceRateF / 100
                    
                    if (env$acceptanceRateFPerc < env$targetAcceptanceFixed) env$rhoF <- env$rhoF - env$rhoF/50
                    if (env$acceptanceRateFPerc > env$targetAcceptanceFixed) env$rhoF <- env$rhoF + env$rhoF/50
                    
                    env$acceptanceRateF <- 0
               }     
                    
               f  <- out[[1]]
               p  <- out[[2]]     
          }
                   
          if (r %% env$gINFOSKIP == 0 | r == 1) {
               progreport(r, p, a, b, d, f, env)
          }
          
     }
     
     # for storing draws
     if (env$gStoreDraws){     
          for (i in 1:env$gNP)
          {
               env$storedDraws[[i]] <- matrix(0, env$gNEREP, env$gNIV)
               dimnames(env$storedDraws[[i]]) <- list(NULL, env$gVarNamesNormal)
          }    
     }
     
     # Iterate after convergence collecting averages 
     n <- env$gNEREP * env$gNSKIP
     
     for (r in 1:n)
     {
          
          if (env$gNIV >0){
                    
               out    <- nextB(a, b, d, p, f, env)
               b      <- out[[1]]
               env$acceptanceRatePerc <- out[[2]]
               p      <- out[[3]]     
               
               a <- nextA(b, d, env)
               
               # this constrains the means of the random parameters to be user-specified
               # this is necessary for such things as error components logit.
               if (!is.null(env$fixedA))
               {
                    for (rp in 1:env$gNIV)
                    {
                         if (!is.na(env$fixedA[rp])) a[rp] <- env$fixedA[rp]
                    }
               } 
               
               if (env$gFULLCV) d <- nextD(a, b, env)
               if (!env$gFULLCV) d <- nextDind(a, b, env)
               
               if (!is.null(env$fixedD))
               {
                    for (rp in 1:env$gNIV)
                    {
                         if (!is.na(env$fixedD[rp])) d[rp,rp] <- env$fixedD[rp]    
                    }
               }   
               
               if (r %% env$gNSKIP == 0) {
                    C <- trans(b,env)
                    env$ma[, r / env$gNSKIP] <- a
                    env$md[, , r / env$gNSKIP] <- d
                    env$mp[, r / env$gNSKIP] <- p^(1/env$TIMES)    #RLH Calculation
                    env$mb <- env$mb + b
                    env$mb.squared <- env$mb.squared + b^2
                    env$mc <- env$mc + C
                    env$mc.squared <- env$mc.squared + C^2
                    if (env$gStoreDraws)
                    {
                         for (i in 1:env$gNP) env$storedDraws[[i]][(r / env$gNSKIP),] <- C[i,]
                    }
               } 
               
          }
          
          if (env$gFIV > 0) {
               # drawing a new set of fixed coefficients
               out <- nextF(p, f, b, env)
               
               if (sum(out[[1]] == f) != env$gFIV) env$acceptanceRateF <- env$acceptanceRateF + 1
               
               # targeting an acceptance rate of 0.25
               if ((r %% 100) == 0)
               {
                    env$acceptanceRateFPerc <- env$acceptanceRateF / 100
                    
                    if (env$acceptanceRateFPerc < env$targetAcceptanceFixed) env$rhoF <- env$rhoF - env$rhoF/50
                    if (env$acceptanceRateFPerc > env$targetAcceptanceFixed) env$rhoF <- env$rhoF + env$rhoF/50
                    
                    env$acceptanceRateF <- 0
               }                          

               f   <- out[[1]]
               p   <- out[[2]]
               
               if (r %% env$gNSKIP == 0) {
                    env$mf[, (r / env$gNSKIP)] <- f
               }
               
          }              
          
          if (r %% env$gINFOSKIP == 0) {
               progreport(env$gNCREP + r, p, a, b, d, f, env)
          }
     }
     
     return(TRUE)
}
