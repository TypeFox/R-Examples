hessianLinear<-function(beta,sigma,y,w,X,variant){

  n<-NROW(w)
  J<-NCOL(w)
  K<-NCOL(X)

  y<-sapply(1:J,function(j) y)

  eta<-sapply(1:J,function(j){
      Xj<-cbind(X[,,j])
      betaj<-beta[,j,drop=FALSE]
      Xj%*%betaj
    })

  h<-(w/sigma)*exp(-0.5*(y-eta)^2/sigma^2)

  g<-rowSums(h)

  dlogP<-function(k,which.param){
    vari<-if (which.param=="beta") variant[k] else FALSE
    dhk<-dh(k,which.param=which.param)
    if (vari) colSums(dhk/g) else sum(rowSums(dhk)/g)
  }

  dh<-function(k,which.param){
    if (which.param=="beta") res <- h*(y-eta)*X[,k,]/sigma^2
    if (which.param=="sigma") res <- -h/sigma+h*(y-eta)^2/sigma^3
    res
  }

  S<-NULL
  for (k in 1:K){
    Sk<-dlogP(k,"beta")
    names(Sk)<-rep(paste("var",k,sep=""),length(Sk))
    S<-c(S,Sk)
  }
  S<-c(S,"sigma"=dlogP(1,"sigma"))

  d2h<-function(k,kprime,which.param1,which.param2){
    if (which.param1=="beta" & which.param2=="beta") res<-X[,k,]*dh(kprime,"beta")*(y-eta)/sigma^2-X[,k,]*X[,kprime,]*h/sigma^2
    if (which.param1=="beta" & which.param2=="sigma") res<-X[,k,]*(dh(1,"sigma")/sigma^2-2*h/sigma^3)*(y-eta)
    if (which.param1=="sigma" & which.param2=="sigma") res<- -(dh(1,"sigma")*sigma-h)/sigma^2+(y-eta)^2*(dh(1,"sigma")*sigma^3-3*h*sigma^2)/sigma^6
     res
  }

  d2logP<-function(k,kprime,which.param1,which.param2){

    variantk<-variant[k]
    variantkprime<-variant[kprime]

    if (which.param1=="beta" & which.param2=="beta"){
      if (variantk & variantkprime){
        res<-matrix(0,J,J)
        for (j in 1:J){
          for (jprime in 1:J){
            dhk<-dh(k,"beta")[,j]
            dhkprime<-dh(kprime,"beta")[,jprime]
            d2hkkprime<-d2h(k,kprime,"beta","beta")[,j]
            if (j!=jprime) res[j,jprime]<- sum(-dhk*dhkprime/g^2)
            if (j==jprime) res[j,jprime]<- sum((d2hkkprime*g-dhk*dhkprime)/g^2)
          }
        }
      }
      if (variantk & !variantkprime){
        res<-matrix(NA,J,1)
        for (j in 1:J){
          dhk<-dh(k,"beta")[,j]
          dhkprime<-rowSums(dh(kprime,"beta"))
          d2hkkprime<-d2h(k,kprime,"beta","beta")[,j]
          res[j,1]<-sum((d2hkkprime*g-dhk*dhkprime)/g^2)
        }
        res<-cbind(res)
      }
      if (!variantk & variantkprime){
        res<-matrix(NA,1,J)
        for (j in 1:J){
          dhk<-rowSums(dh(k,"beta"))
          dhkprime<-dh(kprime,"beta")[,j]
          d2hkkprime<-d2h(k,kprime,"beta","beta")[,j]
          res[1,j]<-sum((d2hkkprime*g-dhk*dhkprime)/g^2)
        }
        res<-rbind(res)
      }
      if (!variantk & !variantkprime){
        dhk<-rowSums(dh(k,"beta"))
        dhkprime<-rowSums(dh(kprime,"beta"))
        d2hkkprime<-rowSums(d2h(k,kprime,"beta","beta"))
        res<-cbind(sum((d2hkkprime*g-dhk*dhkprime)/g^2))
      }
    }

    if (which.param1=="beta" & which.param2=="sigma"){
      if (variantk){
        res<-matrix(0,J,1)
        for (j in 1:J){
          dhbeta<-dh(k,"beta")[,j]
          dhsigma<-rowSums(dh(1,"sigma"))
          d2hbetasigma<-d2h(k,1,"beta","sigma")[,j]
          res[j,1]<- sum((d2hbetasigma*g-dhbeta*dhsigma)/g^2)
        }
      }
      if (!variantk){
        dhbeta<-rowSums(dh(k,"beta"))
        dhsigma<-rowSums(dh(1,"sigma"))
        d2hbetasigma<-rowSums(d2h(k,kprime,"beta","sigma"))
        res<- sum((d2hbetasigma*g-dhbeta*dhsigma)/g^2)
      }
    }

    if (which.param1=="sigma" & which.param2=="sigma"){
      dhsigma<-rowSums(dh(1,"sigma"))
      d2hsigma2<-rowSums(d2h(k,kprime,"sigma","sigma"))
      res<- sum((d2hsigma2*g-dhsigma^2)/g^2)
    }

    res

  }

  Hbeta<-NULL
  for (k in 1:K){
    Hfilak<-NULL
    for (kprime in 1:K){
      res<-d2logP(k,kprime,"beta","beta")
      colnames(res)<-rep(paste("var",kprime,sep=""),NCOL(res))
      rownames(res)<-rep(paste("var",k,sep=""),NROW(res))
      Hfilak<-cbind(Hfilak,res)
    }
    Hbeta<-rbind(Hbeta,Hfilak)
  }

  Hbetasigma<-NULL
  for (k in 1:K) Hbetasigma<-rbind(Hbetasigma,d2logP(k,1,"beta","sigma"))

  Hsigma2<-d2logP(1,1,"sigma","sigma")

  H<-rbind(cbind(Hbeta,Hbetasigma),cbind(t(Hbetasigma),Hsigma2))
  colnames(H)[ncol(H)]<-"sigma"
  rownames(H)[nrow(H)]<-"sigma"

  return(list(S=S,H=H))

}
