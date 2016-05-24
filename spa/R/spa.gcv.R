"spa.gcv"<-
  function(D,y,type,lev,kernel,control,L,U,...){
    ## squared error cv -1
    m=length(L)
    n=dim(D)[1]    
    fit=rep(NA,n)
    p1<-length(levels(as.factor(y)))==2
    p2<-control$pce
    eps=control$adjust
    adjust=p1&p2
    ldepth=control$ldepth
    a=density(D[D!=Inf])
    lqmin=control$lqmin
    lqmax=control$lqmax
    lgrid=control$lgrid
    lmin<-max(quantile(a$x,lqmin),control$ltmin)
    lmax<-quantile(a$x,lqmax)
    if(!is.null(lgrid)){
      lam=seq(lmin,lmax,length=lgrid)
    }else{          
      gcvs=rep(0,ldepth)
      lam=rep(0,ldepth)
    }
    sim=control$dissimilar
    if(type=="hard"){
      K=length(lev)
      probs<-matrix(0,nrow=n,ncol=K)
    }
    if(control$gcv=="tGCV"){
      gcv<-function(i){
        W=kernel(D,i)+eps
        Sj=W/apply(W,1,sum)
        vals=ftf.fit(y,Sj[c(L,U),c(L,U)],"soft",control=control,GCV=TRUE)
        f=vals[[1]][1:m]
        tdf<-vals[[2]]
        if(adjust){
          f[f<0.001]=0.001
          f[f>0.999]=0.999
        }
        err=sum( (y-f)^2)
        err/(1-tdf/m)^2
      }
    }
    if(control$gcv=="aGCV"){
      gcv<-function(i){
        W=kernel(D,i)+eps
        ind<-apply(W[U,L],1,sum)!=0
        flag=FALSE
        if(sum(ind)>0){
          flag=TRUE
        }
        if(flag){
          U1=U[ind]
          vec1=apply(W,1,sum)
          vec2=apply(W[U1,L],1,sum)
          vec3=apply(W[L,U1],1,sum)
          V=t(W[U1,L]/vec2)%*%(W[U1,L]/vec1[U1])
          ##          V=diag(vec3/vec1[L])%*%V/apply(V,1,sum)
          V=(vec3/vec1[L])*(V/apply(V,1,sum))
          V[is.nan(V)]=0
          gr<-W[L,L]/vec1[L]+V
          tdf=sum(diag(gr))
          
          if(tdf/m==1)
            flag=FALSE
        }
        if(!flag){
          W=W[L,L]
          vec1<-apply(W,1,sum)
          gr=W/vec1
          tdf=sum(1/vec1)
        }    
        f=as.numeric(gr%*%y)
        if(adjust){
          f[f<0.001]=0.001
          f[f>0.999]=0.999
        }
        err=sum( (y-f)^2)
        err/(1-tdf/m)^2
      }
    }
    if(control$gcv=="fGCV"){
      gcv<-function(i){
        W=kernel(D,i)+eps
        ind<-apply(W[U,L],1,sum)!=0
        if(sum(ind)>0){
          U1=U[ind]
          vec1=apply(W,1,sum)
          vec2=apply(W[U1,L],1,sum)
          vec3=apply(W[L,U1],1,sum)
          V=t(W[U1,L]/vec2)%*%(W[U1,L]/vec1[U1])
          V=(vec3/vec1[L])*(V/apply(V,1,sum)) 
          V[is.nan(V)]=0
          gr<-W[L,L]/vec1[L]+V
          tdf=sum(1/vec1)
        }else{
          W=W[L,L]
          vec1=apply(W,1,sum)
          gr=W/vec1
          tdf=sum(1/vec1)
          n=m
        }
        f=as.numeric(gr%*%y)
        if(adjust){
          f[f<0.001]=0.001
          f[f>0.999]=0.999
        }
        err=sum( (y-f)^2)
        err/(1-tdf/n)^2
      }
    }
    if(control$gcv=="lGCV"){
      gcv<-function(i){
        W=kernel(D[L,L],i)+eps
        gr=W/apply(W,1,sum)  
        tdf=sum(diag(gr))
        err=sum( (y-gr%*%y)^2)
        err/(1-tdf/m)^2
      }
    }
    if(sim){
      if(!is.null(lgrid)){
        gcvs=sapply(lam,function(i)gcv(i))
      }else{
        gcvs[2]=gcv(lmax)
        gcvs[1]=gcv(lmin)

        l1=lmin
        l2=lmax
        g1=gcvs[1]
        g2=gcvs[2]
        if(is.nan(g1)|is.nan(g2)){
          val=sum(diag(kernel(D[L,L],lmin)))-m==0
          warning("Cv-estimation problems .... diagnose ...\n")
          if(val){
           warning(" W[L,L]==diag(m) which means lGCV and aGCV are undefined, use tGCV or fGCV in this case\n")
          }else{
            warning("problem not known, try other GCV type...\n")
          }
          gcv.str=list(cvlam=(lmin+lmax)/2,gcv=NaN,lambda=NaN,msg="try tGCV")
          return(gcv.str)
        }
        lam[1]=lmin
        lam[2]=lmax
        for(i in 3:(ldepth)){
          if(g1<g2){
            l2=(l1+l2)/2
            g2=gcv(l2)
            lam[i]=l2
            gcvs[i]=g2
          }else{
            l1=(l1+l2)/2
            g1=gcv(l1)
            lam[i]=l1
            gcvs[i]=g1
          }
        }
      }
      vec=!is.na(gcvs)
      if(sum(vec)==0){
        warning("All GCV(lam) are NA using max(lam) ...",
                "other problems may occur...")
        cvlam=lam[length(lam)]
      }else{
        tlam=lam[vec]
        a2=gcvs[vec]
        m1=(1:length(a2)) [a2==min(a2)]
        if(length(m1)>1)
          m1=sample(m1,1)
        cvlam=tlam[m1]
      }
    }else{
      cvlam=lmin
    }
    
    gcv.str=list(cvlam=cvlam,gcv=gcvs,lambda=lam)
    return(gcv.str)
  }
