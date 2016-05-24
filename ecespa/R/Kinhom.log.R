`Kinhom.log` <-
function(A, lambda=NULL, mod=NULL, lifemark="0", prob=NULL,
                                r=NULL, nsim=99, correction="trans", ngrid=200){
  
   
    dataname <- deparse(substitute(A))
    if (is.null(prob)) probname <- "null" else probname <- deparse(substitute(prob))
   Avivos <- unmark(A[A$marks == lifemark])
   Atotal <- unmark(A)
   
   # if there is not probability vector, give everyone the same probability
   # (get it from the frequency of "live" cases):
   
   if (is.null(prob)) prob <- Avivos$n / Atotal$n
   
   # if there is no lambda nor mod, generate a constant lambda
   # this is useful for homogeneous patterns:
   if (is.null(lambda) & is.null (mod)){
         mod <- ppm(Avivos, ~1, Poisson())
         lambda <- predict.ppm(mod, type="trend", ngrid=ngrid)
   }
   
   # if there is not lambda and there is mod, predict lambda from mod
   if (is.null(lambda) & !is.null(mod)){
         mod <- update(mod, Q=Avivos)
         lambda <- predict.ppm(mod, type="trend", ngrid=ngrid)
    }
   
   # if there are both lambda and mod, use mod and forget lambda
    if (!is.null(lambda) & !is.null(mod)){
         mod <- update(mod, Q=Avivos)
         lambda <- predict.ppm(mod, type="trend", ngrid=ngrid)
    }

# if there is lambda but not model, approximate the model with lambda as a covariate
    if (!is.null(lambda) & is.null(mod)){
         mod <- ppm(Avivos,~1 + lam, covariates=list(lam=lambda))
 lambda <- predict.ppm(mod, type="trend", ngrid=ngrid)
    }

    ## rec<e1>lculo de los que est<e1>n dentro de la ventana que genera density.ppp
    ## son los problemas t<ed>picos de decordinaci<f3>n con inside.owin en la versi<f3>n
    ## 1.11.1 de spatstat

    dentro <- !is.na(lambda[Avivos,drop=FALSE])
    Avivos <- Avivos[dentro]

    ## FIN DEL REC<c1>LCULO/CORRECCI<d3>N

    Kao <- Kinhom(Avivos , lambda, correction=correction, r=r)
    mi.r <- Kao$r
    Kao <- Kao[[3]] # 3 es el lugar de la lista resultante de Kinhom que ocupa la K observada
    Kas <- NULL
    
    ##simulations
    for(i in 1:nsim){
      progressreport(i,nsim)
      Avivos <- rthin(Atotal, prob)  
      mod <- update(mod, Q=Avivos)
      lambda <- predict.ppm(mod, type="trend")
      
      dentro <- !is.na(lambda[Avivos,drop=FALSE])
      Avivos <- Avivos[dentro]
      
      Kas=cbind(Kas, Kinhom(Avivos , lambda,  correction=correction, 
                             r=mi.r)[[3]]) 

   }
   result <- (list(r=mi.r, kia = Kao, kia.s=Kas,  datanamea=dataname, type="Kinhom.log", probname=probname, modtrend=mod$trend, nsim=nsim))
   class(result)<-c("ecespa.kci", class(result))
   return(result)
   
}

