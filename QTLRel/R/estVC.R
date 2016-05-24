
machineEps<- (.Machine$double.eps)^(2/4)
inf<- max(1e+38,sqrt(.Machine$double.xmax))

# extract info from specified variance components
fv <- function(vv){
# vv: list of variance components -- list(AA,DD,HH,AD,MH,EE,...) (Abney 200 pp635)
   nms<- c("AA","DD","HH","AD","MH","EE")
   if(any(!is.element(nms,names(vv)))){
      cat("Assume components are in the order AA,DD,HH,AD,MH,EE...\n")
      names(vv)[1:6]<- nms
   }
   vvTmp<- vv
      vvTmp$AA<- vvTmp$DD<- vvTmp$HH<- vvTmp$AD<- vvTmp$MH<- vvTmp$EE<- NULL
   vv<- list(AA=vv$AA,
             DD=vv$DD,
             HH=vv$HH,
             AD=vv$AD,
             MH=vv$MH,
             EE=vv$EE)
      vv<- c(vv,vvTmp)
   if(is.null(vv[[nms[2]]])){
      if(!is.null(vv[[nms[3]]])){
         vv[nms[3]]<- list(NULL)
         cat(nms[3], "is set to null because", nms[2], "is null.\n")
      }
      if(!is.null(vv[[nms[5]]])){
         vv[nms[5]]<- list(NULL)
         cat(nms[5], "is set to null because", nms[2], "is null.\n")
      }
   }
   if(is.null(vv[[nms[1]]])){
      if(!is.null(vv[[nms[4]]])){
        vv[nms[4]]<- list(NULL)
        cat(nms[4], "is set to null because", nms[1], "is null.\n")
      }
   }
   if(is.null(vv[[nms[3]]])){
      if(!is.null(vv[[nms[4]]])){
         vv[nms[4]]<- list(NULL)
         cat(nms[4], "is set to null because", nms[3], "is null.\n")
      }
   }
   nv<- length(vv)
   if(nv>0){
      nnl<- NULL
      for(i in 1:nv){
         nnl<- c(nnl,!is.null(vv[[i]]))
         if(!is.null(vv[[i]])){
            if(!all(is.finite(vv[[i]])))
               stop("Only finite numerical elements are allowed in variance matrices!")
         }
      }
      nn<- cumsum(nnl)
   }else stop("At least the environmental variance component should be included.\n")

   list(v=vv,nv=nv,nnl=nnl,nn=nn)
}

# one of the main functions, Nelder-Mead method
estVC <-
   function(y,
            x,
            v = vector("list",6),
            initpar,
            nit = 25,
            method = c("Nelder-Mead", "BFGS", "CG", "SANN"),
            control = list(),
            hessian = FALSE)
{
   UseMethod("estVC")
}

estVC.default <-
   function(y,
            x,
            v = vector("list",6),
            initpar,
            nit = 25,
            method = c("Nelder-Mead", "BFGS", "CG", "SANN"),
            control = list(),
            hessian = FALSE)
{
   if(!all(is.finite(y)))
      stop("y: non-numeric or infinite data points not allowed.")
   if(!missing(x))
      if(any(sapply(x,is.infinite) | sapply(x,is.na)))
         stop("x: missing or infinite data points not allowed.")

   if(!is.null(dim(y))){
      if(length(dim(y))>2) stop("y: may be wrong.\n")
      if(dim(y)[2]>1)
         warning("y: only the fisrt column will be analyzed.")
      y<- y[,1]
   }
   if(!missing(x)){
      oTmp<- data.frame(y=y,x)
   }else oTmp<- data.frame(y=y)
   oTmp<- model.frame(y~.,oTmp)
   y<- model.response(oTmp)
   x<- model.matrix(y~.,oTmp)
   method<-  match.arg(method)

   estVC.4(y = y,
           x = x,
           v = v,
           initpar = initpar,
           nit = nit,
           method = method,
           control = control,
           hessian = hessian)
}

estVC.1 <-
   function(y,
            x,
            v,
            initpar,
            nit,
            method,
            control,
            hessian)
{
# estimate all background genetic variance (bgv)
# y: vector, response
# x: covariates
# v: list of variance components -- list(AA,DD,HH,AD,MH,EE,...) (Abney 200 pp635)
# initpar: initial parameters, will be initilized automatically if missing
# nit: number of iterations to call optim()
# method: the method to be used
# control: A list of control parameters
# hessian: logical. should a numerically differentiated Hessian matrix be returned?
   fs<- control$fnscale
      if(!is.null(fs) && fs<0) fs<- -1 else fs<- 1

   ny<- length(y)
   nb<- ncol(x)
   ov<- fv(v)
   if(missing(initpar)) initpar<- c(rep(mean(y),ncol(x)),rep(var(y),sum(ov$nnl)))
   for(i in 1:ov$nv){
      if(ov$nnl[i]){
         initpar[nb+ov$nn[i]]<- initpar[nb+ov$nn[i]]/mean(diag(ov$v[[i]]))
         if(i!=4) initpar[nb+ov$nn[i]]<- log(initpar[nb+ov$nn[i]])
      }
   }

   optfct<- function(par,a=list(nb=nb,ny=ny,ov=ov)){
      b<- par[1:a$nb]
      if(a$ov$nnl[4]){
         tmp<- sqrt(exp(par[a$nb+a$ov$nn[1]]+par[a$nb+a$ov$nn[3]])/2)
         if(abs(par[a$nb+a$ov$nn[4]])>tmp) par[a$nb+a$ov$nn[4]]<- sign(par[a$nb+a$ov$nn[4]])*tmp
      }

      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(par[a$nb+a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*par[a$nb+a$ov$nn[i]]
         }
      }
      if(!all(is.finite(S))) return(fs*inf)

      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(fs*inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(fs*inf)

      fs*tmp
   }

   oo<- list(par=initpar)
   val1<- val2<- inf
   while(nit>0){
      oo<- optim(oo$par,optfct,gr=NULL,a=list(nb=nb,ny=ny,ov=ov),method=method,control=control,hessian=FALSE)
      if(ov$nnl[4]){
         tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
         if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
      }

      nit<- nit-1
      val1<- val2
      val2<- oo$value
      if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
         break
   }

   oo<- optim(oo$par,optfct,gr=NULL,a=list(nb=nb,ny=ny,ov=ov),method="Nelder-Mead",control=control,hessian=hessian)
      if(ov$nnl[4]){
         tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
         if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
      }
      for(i in 1:ov$nv){
         if(ov$nnl[i]){
            if(i!=4) oo$par[nb+ov$nn[i]]<- exp(oo$par[nb+ov$nn[i]])
         }
      }
   oo$value<- as.numeric(oo$value)

   oo$value<- -oo$value #-fs*oo$value
   attributes(oo$value)<- NULL
   names(oo$par)[1:nb]<- colnames(x)
   names(oo$par)[nb+1:sum(ov$nnl)]<- names(ov$v[ov$nnl])
   oo$y<- as.matrix(y)
   oo$x<- as.matrix(x)
   oo$v<- v
   oo$nv<- ov$nv
   oo$nnl<- ov$nnl
   oo$nn<- ov$nn
   class(oo)<- "bgv"
   oo
}

# NOTES: as accurate as but about 2.5 times as fast as estVC.1
estVC.2 <-
   function(y,
            x,
            v,
            initpar,
            nit,
            method,
            control,
            hessian)
{
# estimate all background genetic variance (bgv)
# y: vector, response
# x: covariates
# v: list of variance components -- list(AA,DD,HH,AD,MH,EE,...) (Abney 200 pp635)
# initpar: initial parameters, will be initilized automatically if missing
# nit: number of iterations to call optim()
# method: the method to be used
# control: A list of control parameters
# hessian: logical. should a numerically differentiated Hessian matrix be returned?
   fs<- control$fnscale
      if(!is.null(fs) && fs<0) fs<- -1 else fs<- 1

   ny<- length(y)
   nb<- ncol(x)
   ov<- fv(v)
   if(missing(initpar)) initpar<- c(rep(mean(y),ncol(x)),rep(var(y),sum(ov$nnl)))
   for(i in 1:ov$nv){
      if(ov$nnl[i]){
         initpar[nb+ov$nn[i]]<- initpar[nb+ov$nn[i]]/mean(diag(ov$v[[i]]))
         if(i!=4) initpar[nb+ov$nn[i]]<- log(initpar[nb+ov$nn[i]])
      }
   }

   optfct<- function(par,a=list(nb=nb,ny=ny,ov=ov)){
      b<- par[1:a$nb]
      if(a$ov$nnl[4]){
         tmp<- sqrt(exp(par[a$nb+a$ov$nn[1]]+par[a$nb+a$ov$nn[3]])/2)
         if(abs(par[a$nb+a$ov$nn[4]])>tmp) par[a$nb+a$ov$nn[4]]<- sign(par[a$nb+a$ov$nn[4]])*tmp
      }

      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(par[a$nb+a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*par[a$nb+a$ov$nn[i]]
         }
      }
      if(!all(is.finite(S))) return(fs*inf)

      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(fs*inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(fs*inf)

      fs*tmp
   }
   optfct.b<- function(a=list(par=initpar,nb=nb,ny=ny,ov=ov)){
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(a$par[a$nb+a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*a$par[a$nb+a$ov$nn[i]]
         }
      }
      if(!all(is.finite(S))) return(a$par[1:a$nb])

      dd<- eigen(S,symmetric=T)
         uu<- dd$vec
         dd<- abs(dd$val)
         dd[dd<machineEps]<- machineEps
      yy<- t(uu)%*%y
         yy<- as.matrix(yy)
      xx<- t(uu)%*%x
         xx<- as.matrix(xx)
      b<- lm.wfit(x=xx,y=yy,w=1/dd)$coef
         b[!is.finite(b)]<- 0
      b
   }
   optfct.v<- function(par,a=list(par=initpar,nb=nb,ny=ny,ov=ov)){
      b<- a$par[1:a$nb]
      if(a$ov$nnl[4]){
         tmp<- sqrt(exp(par[a$ov$nn[1]]+par[a$ov$nn[3]])/2)
         if(abs(par[a$ov$nn[4]])>tmp) par[a$ov$nn[4]]<- sign(par[a$ov$nn[4]])*tmp
      }

      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(par[a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*par[a$ov$nn[i]]
         }
      }
      if(!all(is.finite(S))) return(fs*inf)

      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(fs*inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(fs*inf)

      fs*tmp
   }

   oo<- list(par=initpar)
   val1<- val2<- inf
   pparTmp<- initpar
   if(length(initpar) < nb+2){
      while(nit>0){
         pparTmp[1:nb]<- optfct.b(a=list(par=pparTmp,nb=nb,ny=ny,ov=ov))
         oo<- optimize(optfct.v,interval=c(-1000,1000),a=list(par=pparTmp,nb=nb,ny=ny,ov=ov),maximum=ifelse(fs<0,TRUE,FALSE))
            pparTmp[-c(1:nb)]<- oo$objective
            oo<- list(par=pparTmp, value=oo$minimum)
            if(ov$nnl[4]){
               tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
               if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
            }
         pparTmp<- oo$par

         val1<- val2
         val2<- oo$value
         if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
            break
         nit<- nit-1
      }
   }else{
      while(nit>0){
         pparTmp[1:nb]<- optfct.b(a=list(par=pparTmp,nb=nb,ny=ny,ov=ov))
         oo<- optim(pparTmp[-c(1:nb)],optfct.v,gr=NULL,a=list(par=pparTmp,nb=nb,ny=ny,ov=ov),
                    method=method,control=control,hessian=FALSE)
            pparTmp[-c(1:nb)]<- oo$par
            oo$par<- pparTmp
            if(ov$nnl[4]){
               tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
               if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
            }
         pparTmp<- oo$par

         val1<- val2
         val2<- oo$value
         if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
            break
         nit<- nit-1
      }
   }
   oo<- optim(oo$par,optfct,gr=NULL,a=list(nb=nb,ny=ny,ov=ov),method="Nelder-Mead",control=control,hessian=hessian)
      if(ov$nnl[4]){
         tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
         if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
      }
      for(i in 1:ov$nv){
         if(ov$nnl[i]){
            if(i!=4) oo$par[nb+ov$nn[i]]<- exp(oo$par[nb+ov$nn[i]])
         }
      }

   oo$value<- as.numeric(oo$value)
   oo$value<- -oo$value
   attributes(oo$value)<- NULL
   names(oo$par)[1:nb]<- colnames(x)
   names(oo$par)[nb+1:sum(ov$nnl)]<- names(ov$v[ov$nnl])
   oo$y<- as.matrix(y)
   oo$x<- as.matrix(x)
   oo$v<- ov$v
   oo$nv<- ov$nv
   oo$nnl<- ov$nnl
   oo$nn<- ov$nn
   class(oo)<- "bgv"
   oo
}

# NOTES: as accurate as but about 5 times as fast as estVC.1
# as accurate as but about 2 times as fast as estVC.2; may Not stable!!!
estVC.3 <-
   function(y,
            x,
            v,
            initpar,
            nit,
            method,
            control,
            hessian)
{
# estimate all background genetic variance (bgv)
# y: vector, response
# x: desig matrix including overall mean !!!
# v: list of variance components -- list(AA,DD,HH,AD,MH,EE,...) (Abney 200 pp635)
# initpar: initial parameters, will be initilized automatically if missing
# nit: number of iterations
   control$fnscale<- 1
   ny<- length(y)
   nb<- ncol(x)
   ov<- fv(v)
   if(missing(initpar)) initpar<- c(rep(mean(y),ncol(x)),rep(var(y),sum(ov$nnl)))
   for(i in 1:ov$nv){
      if(ov$nnl[i]){
         initpar[nb+ov$nn[i]]<- initpar[nb+ov$nn[i]]/mean(diag(ov$v[[i]]))
         if(i!=4) initpar[nb+ov$nn[i]]<- log(initpar[nb+ov$nn[i]])
      }
   }

   optfct<- function(par,a=list(nb=nb,ny=ny,ov=ov)){
      b<- par[1:a$nb]
      if(a$ov$nnl[4]){
         tmp<- sqrt(exp(par[a$nb+a$ov$nn[1]]+par[a$nb+a$ov$nn[3]])/2)
         if(abs(par[a$nb+a$ov$nn[4]])>tmp) par[a$nb+a$ov$nn[4]]<- sign(par[a$nb+a$ov$nn[4]])*tmp
      }

      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(par[a$nb+a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*par[a$nb+a$ov$nn[i]]
         }
      }
      if(!all(is.finite(S))) return(inf)

      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(inf)

      tmp
   }
   optfct.b<- function(a=list(par=par,nb=nb,ny=ny,ov=ov)){
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(a$par[a$nb+a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*a$par[a$nb+a$ov$nn[i]]
         }
      }
      if(!all(is.finite(S))) return(a$par[1:a$nb])

      dd<- eigen(S,symmetric=T)
         uu<- dd$vec
         dd<- abs(dd$val)
         dd[dd<machineEps]<- machineEps
      yy<- t(uu)%*%y
         yy<- as.matrix(yy)
      xx<- t(uu)%*%x
         xx<- as.matrix(xx)
      b<- lm.wfit(x=xx,y=yy,w=1/dd)$coef
         b[!is.finite(b)]<- 0
      b
   }
   optfct.v<- function(par,a=list(par=par,nb=nb,ny=ny,ov=ov)){
      b<- a$par[1:a$nb]
      if(a$ov$nnl[4]){
         tmp<- sqrt(exp(par[a$ov$nn[1]]+par[a$ov$nn[3]])/2)
         if(abs(par[a$ov$nn[4]])>tmp) par[a$ov$nn[4]]<- sign(par[a$ov$nn[4]])*tmp
      }

      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(par[a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*par[a$ov$nn[i]]
         }
      }
      if(!all(is.finite(S))) return(inf)

      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(inf)

      tmp
   }

   val1<- val2<- inf
   pparTmp<- initpar
   while(nit>0){
      pparTmp[1:nb]<- optfct.b(a=list(par=pparTmp,nb=nb,ny=ny,ov=ov))
      oo<- nlm(optfct.v,pparTmp[-c(1:nb)],a=list(par=pparTmp,nb=nb,ny=ny,ov=ov))
         pparTmp[-c(1:nb)]<- oo$estimate; oo$value<- oo$minimum
         oo$par<- pparTmp
         if(ov$nnl[4]){
            tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
            if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
         }
      pparTmp<- oo$par
      val1<- val2
      val2<- oo$value
      if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
         break
      nit<- nit-1
   }
   oo<- optim(oo$par,optfct,gr=NULL,a=list(nb=nb,ny=ny,ov=ov),method="Nelder-Mead",control=control,hessian=hessian)
      if(ov$nnl[4]){
         tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
         if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
      }
      for(i in 1:ov$nv){
         if(ov$nnl[i]){
            if(i!=4) oo$par[nb+ov$nn[i]]<- exp(oo$par[nb+ov$nn[i]])
         }
      }

   oo$value<- as.numeric(oo$value)
   oo$value<- -oo$value
   attributes(oo$value)<- NULL
   names(oo$par)[1:nb]<- colnames(x)
   names(oo$par)[nb+1:sum(ov$nnl)]<- names(ov$v[ov$nnl])
   oo$y<- as.matrix(y)
   oo$x<- as.matrix(x)
   oo$v<- ov$v
   oo$nv<- ov$nv
   oo$nnl<- ov$nnl
   oo$nn<- ov$nn
   class(oo)<- "bgv"
   oo
}

# the above ones give more optional but less sensible solutions
estVC.4 <-
   function(y,
            x,
            v,
            initpar,
            nit,
            method,
            control,
            hessian)
{
# estimate all background genetic variance (bgv)
# y: vector, response
# x: covariates
# v: list of variance components -- list(AA,DD,HH,AD,MH,EE,...) (Abney 200 pp635)
# initpar: initial parameters, will be initilized automatically if missing
# nit: number of iterations to call optim()
# method: the method to be used
# control: A list of control parameters
# hessian: logical. should a numerically differentiated Hessian matrix be returned?
   fs<- control$fnscale
      if(!is.null(fs) && fs<0) fs<- -1 else fs<- 1

   ny<- length(y)
   nb<- ncol(x)
   ov<- fv(v)
   if(missing(initpar)) initpar<- c(rep(mean(y),ncol(x)),rep(var(y),sum(ov$nnl)))
   for(i in 1:ov$nv){
      if(ov$nnl[i]){
         initpar[nb+ov$nn[i]]<- initpar[nb+ov$nn[i]]/mean(diag(ov$v[[i]]))
         if(i!=4) initpar[nb+ov$nn[i]]<- log(initpar[nb+ov$nn[i]])
      }
   }

   optfct<- function(par,a=list(nb=nb,ny=ny,ov=ov)){
      b<- par[1:a$nb]
      if(a$ov$nnl[4]){
         tmp<- sqrt(exp(par[a$nb+a$ov$nn[1]]+par[a$nb+a$ov$nn[3]])/2)
         if(abs(par[a$nb+a$ov$nn[4]])>tmp) par[a$nb+a$ov$nn[4]]<- sign(par[a$nb+a$ov$nn[4]])*tmp
      }

      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(par[a$nb+a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*par[a$nb+a$ov$nn[i]]
         }
      }
      if(!all(is.finite(S))) return(fs*inf)

      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(fs*inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(fs*inf)

      fs*tmp
   }

   optfct.b<- function(a=list(par=par,nb=nb,ny=ny,ov=ov)){
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(a$par[a$nb+a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*a$par[a$nb+a$ov$nn[i]]
         }
      }
      if(!all(is.finite(S))) return(a$par[1:a$nb])

      dd<- eigen(S,symmetric=T)
         uu<- dd$vec
         dd<- abs(dd$val)
         dd[dd<machineEps]<- machineEps
      yy<- t(uu)%*%y
         yy<- as.matrix(yy)
      xx<- t(uu)%*%x
         xx<- as.matrix(xx)
      b<- lm.wfit(x=xx,y=yy,w=1/dd)$coef
         b[!is.finite(b)]<- 0
      b
   }
   optfct.v<- function(par,a=list(par=par,nb=nb,ny=ny,ov=ov)){
      b<- a$par[1:a$nb]
      if(a$ov$nnl[4]){
         tmp<- sqrt(exp(par[a$ov$nn[1]]+par[a$ov$nn[3]])/2)
         if(abs(par[a$ov$nn[4]])>tmp) par[a$ov$nn[4]]<- sign(par[a$ov$nn[4]])*tmp
      }

      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
         if(a$ov$nnl[i]){
            if(i!=4){
               S<- S + a$ov$v[[i]]*exp(par[a$ov$nn[i]])
            }else S<- S + a$ov$v[[i]]*par[a$ov$nn[i]]
         }
      }
      if(!all(is.finite(S))) return(inf)

      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(inf)

      tmp
   }

   upper<- rep(25,length(ov$v))
      names(upper)<- names(ov$v)
      upper["AD"]<- Inf
      upper<- upper[ov$nnl]
      upper<- c(rep(Inf,nb),upper)
   val1<- val2<- inf
   pparTmp<- initpar
   while(nit>0){
      pparTmp[1:nb]<- optfct.b(a=list(par=pparTmp,nb=nb,ny=ny,ov=ov))
      oo<- nlminb(pparTmp[-c(1:nb)],optfct.v,a=list(par=pparTmp,nb=nb,ny=ny,ov=ov),upper=upper[-c(1:nb)])
         pparTmp[-c(1:nb)]<- oo$par; oo$value<- oo$objective
         oo$par<- pparTmp
         if(ov$nnl[4]){
            tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
            if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
         }
#         for(i in 1:ov$nv){
#            if(ov$nnl[i]){
#               if(i!=4) oo$par[nb+ov$nn[i]]<- min(10,exp(oo$par[nb+ov$nn[i]]))
#            }
#         }
      pparTmp<- oo$par
      val1<- val2
      val2<- oo$value
      if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
         break
      nit<- nit-1
   }

   oo<- nlminb(oo$par,optfct,gradient=NULL,a=list(nb=nb,ny=ny,ov=ov),upper=upper)
      oo$value<- oo$objective
      if(ov$nnl[4]){
         tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
         if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
      }
   oo<- optim(oo$par,optfct,gr=NULL,a=list(nb=nb,ny=ny,ov=ov),
              method="Nelder-Mead",control=list(maxit=1),hessian=hessian)
      for(i in 1:ov$nv){
         if(ov$nnl[i]){
            if(i!=4) oo$par[nb+ov$nn[i]]<- exp(oo$par[nb+ov$nn[i]])
         }
      }

   oo$value<- -fs*oo$value
   attributes(oo$value)<- NULL
   names(oo$par)[1:nb]<- colnames(x)
   names(oo$par)[nb+1:sum(ov$nnl)]<- names(ov$v[ov$nnl])
   oo$y<- as.matrix(y)
   oo$x<- as.matrix(x)
   oo$v<- ov$v
   oo$nv<- ov$nv
   oo$nnl<- ov$nnl
   oo$nn<- ov$nn
   class(oo)<- "bgv"
   oo
}

blup<- function(object)
{
   UseMethod("blup")
}

blup.bgv<- function(object){
# best linear unbiased prediction (BLUP) for all effects
# y: ny by 1 matrix
# x: design matrix including overall mean !!!
# object: object from estBVG
# v: list of variance components corresponding object$par[-c(1:nb)]
   ny<- nrow(object$y)
   nb<- ncol(object$x)
   nv<- length(object$v)
   idx<- NULL
   for(i in 1:nv)
      if(!is.null(object$v[[i]])) idx<- c(idx,i)
   object$v<- object$v[idx]
   nv<- length(object$v)

   S<- matrix(0,nrow=ny,ncol=ny)
   for(i in 1:nv)
      S<- S + object$v[[i]]*object$par[nb+i]

   out<- vector("list",nv)
   rr<- object$y-object$x%*%object$par[1:nb]
   rd<- solve(S)%*%rr
   out[[1]]<- sweep(object$x,2,object$par[1:nb],"*")
      names(out)[1]<- "fixed"
   for(i in 1:nv){
      out[[1+i]]<- object$v[[i]]%*%rd*object$par[nb+i]
      names(out)[1+i]<- names(object$v)[i]
   }

   out
}

print.bgv<- function(x,...){
   cat("\nvalue:\n"); print(x$value)
   cat("\nparameters:\n"); print(x$par)
}

