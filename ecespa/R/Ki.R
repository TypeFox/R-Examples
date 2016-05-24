`Ki` <-
function(mod1,correction="trans", nsim=99, ngrid=200,
                         nrep=1e5, r=NULL, spctype=1)
{
## datos b<e1>sicos
 modnamea <- deparse(substitute(mod1))
 
if(inherits(mod1,"ppm")){  
   I.ppp= mod1$Q$data
   lambdaI = predict(ppm(mod1$Q$data, mod1$trend), type="trend", ngrid=ngrid)
   Isim="ppm"
   dataname.a <- mod1$call[[2]]
}
if(inherits(mod1,"ecespa.minconfit")){  
   I.ppp= mod1$data
   lambdaI = mod1$lambda
   Isim="spc"
   dataname.a <- mod1$dataname
}


## correcci<f3>n de los problemas de descoordinaci<f3>n entre inside.owin y density en
## la versi<f3>n 1.11.1 de spatstat( hay que asegurarse de que no hay puntos fuera de 
## las ventanas de intensidad
 #      dentro = !is.na(lambdaI[I.ppp, drop=FALSE])
 #      I.ppp=I.ppp[dentro]
 #      dentro = !is.na(lambdaJ[J.ppp, drop=FALSE])
 #      J.ppp=J.ppp[dentro]
##### FIN DE LA CORRECCI<d3>N

  ## C<e1>lculo de los Ki de cada patr<f3>n
   Kia=Kinhom(I.ppp, lambdaI, correction=correction, r=r)
   mi.r = Kia$r
   Kia = Kia[[3]]
   Kia.s= NULL
   Isim.ppp = I.ppp
      for (i in 1: nsim){ 
           progressreport(i,nsim)
        ## aseguramos que no haya NAs en el vector simulado de lambdas
            if(Isim=="ppm"){      Isim.ppp = rmh(mod1, start=list(x.start=I.ppp), 
                                                control=list(p=1, nrep=nrep), verbose=FALSE) }
            else if (Isim=="spc") Isim.ppp = rIPCP (mod1, type=spctype)

            dentro = !is.na(lambdaI[Isim.ppp, drop=FALSE]) 
            Isim.ppp = Isim.ppp[dentro]
            ## K simuladas
            Kia.s= cbind(Kia.s, Kinhom(Isim.ppp, lambdaI, 
                           correction=correction, r = mi.r, nlarge=Inf)[[3]])
    }
    
   result=(list(r=mi.r, kia = Kia, kia.s=Kia.s, lambda=lambdaI, datanamea=dataname.a, modnamea=modnamea, type="Ki"))
   class(result)<-c("ecespa.kci", class(result))
   return(result)
}

