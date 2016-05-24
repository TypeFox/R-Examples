`Kmm` <-
function (mippp, r = 1:10, nsim=NULL) 
{
    
    dataname <- deparse(substitute(mippp))
    mippp1 = setmarks(mippp, 1)
    lambda = mippp$n/area.owin(mippp$window)
    mu = mean(mippp$marks)
    Kmm = NULL
    Kmm1 = NULL
    Kmmsim <- NULL
 

  if(!is.null(nsim)){
         pepesim <- list(mippp)
         pepesim <- rep(pepesim, nsim) 
         pepesim <- lapply(pepesim, function(x)  x = x %mark% sample(x$marks))
         # class(pepesim) <- c(class(pepesim), "splitppp")
        
         m0sim <- sapply(pepesim, markstat, sum, R=0)
 }
       m0 = markstat(mippp, sum, R = 0)
       m01 = markstat(mippp1, sum, R = 0)
       
    


    cat("computing Kmm for distance class ")
    for (i in 1:length(r)) {
    progressreport(i, length(r))
        mr <- markstat(mippp, sum, R = r[i])
        mr1 <- markstat(mippp1, sum, R = r[i])

        sumatorio <- (mr - m0) * m0
        sumatorio1 <- (mr1 - m01) * m01
            
E0 <- mean(sumatorio)
        E01 <- mean(sumatorio1)

        Kmm = c(Kmm, E0/(lambda * mu^2))
        Kmm1 = c(Kmm1, E01/(lambda))

if(!is.null(nsim)){
   mrsim <- sapply(pepesim,markstat,sum, R=r[i])
   sumatoriosim <- (mrsim - m0sim) * m0sim
   E0sim <- apply(sumatoriosim,2, mean)
   Kmmsim <- rbind(Kmmsim,E0sim/(lambda*mu^2))
}

   }
    
   if(!is.null(nsim)) Kmmsim.n <- apply(Kmmsim,2, function(x) x/Kmm1) else Kmmsim.n <-NULL
   
   result=list(dataname=dataname, r = r, nsim=nsim, kmm = Kmm, 
                   kmm.n = Kmm/Kmm1, kmmsim=Kmmsim, 
                   kmmsim.n=Kmmsim.n)

   class(result) =c("ecespa.kmm",class(result))
return(result)   
}

