
####################################
#  functions for model selection:  #
# multivariate multiple regression #
################################################################################

##################
# part I: mtcmim #
##################
Array<- function(k){
  A<- matrix(0,nrow=2^k,ncol=k)
  for(j in 1:k){
    e<- 0
    i<- 0
    for(m in 1:2^j){
      for(n in 1:2^(k-j)){
        i<- i+1
        A[i,j]<- (-1)^e
      }
      e<- e+1
    }
  }
  (A+1)/2
}

# mymtcmim: estimate parameters
mymtcmim<- function(y,W,nws,ws,a,sigma,P,G,ngs,gs,b,init=1,
  iter=10000,tol=1e-12){
# y: n by p matrix
# W: n by nW matrix
# nws: vector of length p, nws_j covariates for y_j, specified by ws
# a: covariate effects
# sigma: residual matrix
# P: n by np matrix, mixing proportions
# G: np by nG matrix, genetic matrix
# ngs: vector of length p, ngs_j QTL for y_j, specified by gs
# b: QTL effects
  if(!is.matrix(y)) stop("\a  y should be a matrix...")
  if(!is.matrix(W)) stop("\a  W should be a matrix...")
  if(!is.matrix(sigma)) stop("\a  sigma should be a matrix...")
  if(missing(P)){
    P<- matrix(1,nrow=n,ncol=1)
    G<- numeric(0)
    ngs<- numeric(0)
    gs<- numeric(0)
    b<- numeric(0)
  }else if(!is.matrix(P)){
    stop("\a  P should be a matrix...")
  }else if(ncol(P)>1){
    if(!is.matrix(G)) stop("\a  G should be a matrix or nothing...")
  }

  n<- nrow(y)
  p<- ncol(y)
  np<- ncol(P)
  nG<- ncol(G)
  nW<- ncol(W)
  loglik<- -1e+308
  out<- .C("mtcmimEstc",
    as.double(t(y)),
    as.integer(n),
    as.integer(p),
    as.double(t(P)),
    as.integer(np),
    as.double(t(G)),
    as.integer(nG),
    as.integer(ngs),
    as.integer(gs),
    as.double(t(W)),
    as.integer(nW),
    as.integer(nws),
    as.integer(ws),
    a=as.double(a),
    b=as.double(b),
    sigma=as.double(t(sigma)),
    loglik=as.double(loglik),
    as.integer(init),
    as.integer(iter),
    as.double(tol))
  list(loglik=out$loglik,a=out$a,b=out$b,
       sigma=matrix(out$sigma,nrow=p,byrow=TRUE))
}

fw<- function(xid){
  np<- length(xid)
  nws<- rep(0,np)
  ws<- numeric(0)
  for(i in 1:np){
    nws[i]<- length(xid[[i]])
    ws<- c(ws,xid[[i]])
  }
  list(nws=nws,ws=ws)
}

fgg<- function(qtl,eps,nQ){
  neps<- nrow(eps)
  gg<- qtl
  if(!is.null(neps)){
    n0<- nQ
    for(i in 1:neps){
      n0<- n0+1
      j<- eps$y[i]
      gg[[j]]<- c(gg[[j]],n0)
    }
  }
  gg
}

fgs<- function(gg){
  np<- length(gg)
  ngs<- rep(0,np)
  gs<- numeric(0)
  for(i in 1:np){
    ngs[i]<- length(gg[[i]])
    gs<- c(gs,gg[[i]])
  }
  list(ngs=ngs,gs=gs)
}

fdist<- function(mpos){
# creat dists from mpos if mpos$dist is not cumulative
  mid<- mpos$id; o<- order(mid,decreasing=F)
  mid<- mid[o]
  ch<- mpos$ch[o]
  d<- mpos$dist[o]
  chrs<- sort(unique(ch),decreasing=F)
  for(i in chrs){
    w<- ch==i
    d[w]<- cumsum(d[w])
  }
  data.frame(mid=mid,ch=ch,d=d)
}

forder<- function(dists){
# the i-th element is the ord[i]-th smallest
  ord<- order(dists$ch,dists$d)
  order(ord)
}

fo<- function(olddist,newdist){
  o1<- forder(olddist)
  o2<- forder(newdist)
  o<- order(o1,decreasing=F)
  o2[o]
}

feps<- function(eps,o){
  neps<- nrow(eps)
  ep<- eps
  if(!is.null(neps)&&neps>0)for(i in 1:neps){
    ep$q1[i]<- o[eps$q1[i]]
    ep$q2[i]<- o[eps$q2[i]]
  }
  ep
}

fqtl<- function(qtl,o){
  p<- length(qtl)
  qt<- vector("list",p)
  for(i in 1:p){
    if(length(qtl[[i]])>0)
      qt[[i]]<- o[qtl[[i]]]
  }
  qt
}

fG<- function(A,eps){#should be consistent with fgg()
  nr<- nrow(A)
  nQ<- ncol(A)
  neps<- nrow(eps)
  if(is.null(neps)){
    G<- A-1/2
  }else if(neps>0){
    nc<- nQ+neps
    G<- matrix(0.0,nrow=nr,ncol=nc)
    G[,1:nQ]<- as.matrix(A-1/2)
    j0<- nQ
    for(j in 1:neps){
      j0<- j0+1
      G[,j0]<- G[,eps$q1[j]]*G[,eps$q2[j]]
    }
  }else stop("\a  eps: wrong input...")
  G
}

fP<- function(A,mpos,mdat,dists,pp=1){
  n<- nrow(mdat)
  nP<- nrow(A)
  nQ<- ncol(A)
  nm<- ncol(mdat)
  P<- matrix(1.0,nrow=n,ncol=nP)
  mid<- sort(unique(dists$mid))
  nmid<- length(mid)
  if(nQ!=nrow(dists))stop("\a  fP: wrong input...")

  out<- .C("fPc",
    as.integer(t(A)),
    as.integer(nP),
    as.integer(nQ),
    as.integer(t(mdat)),
    as.integer(n),
    as.integer(nm),
    as.double(t(mpos)),
    as.integer(dists$ch),
    as.integer(dists$mid),
    as.double(dists$d),
    as.integer(mid),
    as.integer(nmid),
    P=as.double(t(P)),
    as.integer(pp))
  matrix(out$P,nrow=n,byrow=T)
}

fs<- function(dv,win,dists,k,range=0){
# k: k-th row of dists
# range: genome-wide (0), the same chromosome (-1)
  dv<- rbind(dv,c(mid=dists$mid[k],ch=dists$ch[k],dist=dists$d[k]))
  k0<- dists$ch[k]
  d0<- dists$d[k]
  if(!is.numeric(range)) stop("range: wrong...")
  if(length(range)==1){
    if(range==-1){
      ch<- (dv$ch==k0) & (dv$dist>=d0-win)&(dv$dist<=d0+win)
    }else if(range==0){
      ch<- rep(TRUE,nrow(dv))
    }else stop("range: -1 or 0 only.")
  }else stop("range: -1 or 0 only.")
  if(sum(ch)<1) stop("win or range: wrong...")

  dv[ch,]
}

fqtdst<- function(qt,dst,ep=NULL){
   np<- length(qt)
   xx<- unique(unlist(qt))
   if(length(xx)<1){
      list(qtl=NULL,eps=ep,dists=dst)
   }else{
      dst<- dst[order(dst$ch,dst$d),]
      xx<- sort(xx,decreasing=FALSE)
      oo<- order(xx,decreasing=F)
      xxx<- rep(Inf,max(xx))
         xxx[xx]<- oo
      for(j in 1:np){
         if(length(qt[[j]])>0)
            qt[[j]]<- sort(xxx[qt[[j]]])
      }
      if(!is.null(ep)){
         ep$q1<- xxx[ep$q1]
         ep$q2<- xxx[ep$q2]
      }
      list(qtl=qt,eps=ep,dists=dst[xx,])
   }
}

# mtcmim: multiple-trait composite-interval-mapping 
mtcmim<- function(y,mpos,mdat,x,xid,dists,a,b,sigma,
  qtl=NULL,eps=NULL,win=Inf,range=0,pp=1,len=2,init=1,iter=10000,tol=1e-12){
# y: n by p matrix, traits
# x: n by nW matrix, covarites; excluding intercept
# xid: list of length p, xid[[j]] specifies columns of x as covariates for y_j 
# dists: data frame (ch=chrom id,mid=marker id,d=dist on the chrom), specifies QTLs
# a: covarite effects
# b: QTL effects
# sigma: residual variance matrix
# qtl: list of length p, qtl[[j]]--which qtls for y_j
# eps: data frame (y=which trait,q1=,q2=)
# win: window width to search
# range: genome-wide (0), the same chromosome (-1)
# pp: mapping population, BC-1, RIL-selfing-2, RIL-brother-sister-mating-3
# len: move step length in search
  if(is.data.frame(mdat)){
    mdat<- as.matrix(mdat)
  }else if(!is.matrix(mdat))
     stop("mdat: should be a matrix.")
  tmp<- unique(as.matrix(mdat))
    tmp<- setdiff(tmp,c(0,1))
  if(length(tmp)>0) stop("mdat elements: 0 or 1 only.")
  if(!is.null(dists$ch) && !is.numeric(dists$ch))
    stop("dists: chromosome IDs should be integers.")

  if(is.data.frame(y)){
    y<- as.matrix(y)
  }else if(!is.matrix(y))
    stop("\a  y: should be a matrix...")
  ny<- nrow(y); np<- ncol(y)

  if(missing(x)){
    xx<- matrix(1,nrow=ny,ncol=1)
    xid<- vector("list",np)
    for(n in 1:np) xid[[n]]<- 1; rm(n)
  }else{
    if(is.data.frame(x)){
      x<- as.matrix(x)
    }else if(!is.matrix(x))
      stop("\a  x should be a matrix...")
    if(!is.list(xid)) stop("\a  xid: should be a list...")
    xx<- cbind(1,x); xx<- as.matrix(xx)
    for(n in 1:np) xid[[n]]<- c(0,xid[[n]])+1; rm(n)
  }
  if(missing(a)){
    na<- NULL; for(j in 1:np) na<- c(na,length(xid[[j]]))
    a<- rep(0,sum(na))
    na<- c(0,cumsum(na))+1
      na<- na[1:np]
    a[na]<- colMeans(y)
  }
  if(missing(b)){
    nb<- 0; for(j in 1:np) nb<- nb+length(qtl[[j]])
    if(!is.null(eps)) nb<- nb+nrow(eps)
    b<- rep(0,nb)
  }
  if(missing(sigma))
    sigma<- var(y) #diag(np)

  if(!is.null(qtl)){
    if(!is.list(qtl))stop("\a  qtl: should be a list...")
    if(length(qtl)!=np)stop("\a  qtl: wrong input...")
  }
  qtdstTmp<- fqtdst(qt=qtl,dst=dists,ep=eps)
  qtl<- qtdstTmp$qtl
  eps<- qtdstTmp$eps
  dists<- qtdstTmp$dists

  nW<- ncol(xx)
  ob<- fw(xid); nws<- ob$nws; ws<- ob$ws

  if(is.null(qtl)&&!is.null(eps)){
    stop("\a  eps: should among qtl...")
  }else if(is.null(qtl)&&is.null(eps)){
    dists<- NULL
    P<- matrix(1,nrow=ny,ncol=1)
    G<- numeric(0)
    ngs<- numeric(0)
    gs<- numeric(0)
    b<- numeric(0)
    out<- mymtcmim(y=y,W=xx,nws=nws,ws=ws,a=a,sigma=sigma,P=P,G=G,
      ngs=ngs,gs=gs,b=b,init=init,iter=iter,tol=tol)
  }else{
    nQ<- dim(dists)[1] #number of QTLs
    A<- Array(nQ)
    dv<- div(mpos,len=len)

    gg<- fgg(qtl,eps,nQ)
    ob<- fgs(gg); ngs<- ob$ngs; gs<- ob$gs
    P<- fP(A,mpos,mdat,dists,pp)
    G<- fG(A,eps)

    out<- mymtcmim(y=y,W=xx,nws=nws,ws=ws,a=a,sigma=sigma,P=P,G=G,
      ngs=ngs,gs=gs,b=b,init=init,iter=iter,tol=tol)
    lik2<- -Inf
    lik1<- -Inf
    lik<- out$loglik
    la<- Inf

    k<- 0
    mx<- Inf
    lik0<- lik
    dvv<- dv
    while(mx>tol){
      k<- k%%nQ+1;
      dvv<- fs(dv,win,dists,k,range)
      ndvv<- nrow(dvv)
      if(!is.null(ndvv)&&ndvv>1){
        for(i in 1:ndvv){
          dist0<- dists
          dist0$ch[k]<- dvv$ch[i]
          dist0$mid[k]<- dvv$mid[i]
          dist0$d[k]<- dvv$dist[i]
          o<- fo(dists,dist0)
          eps0<- feps(eps,o)
          qtl0<- fqtl(qtl,o)
          gg0<- fgg(qtl0,eps0,nQ)
          dvv0<- fgs(gg0); ngs0<- dvv0$ngs; gs0<- dvv0$gs
          P0<- fP(A,mpos,mdat,dist0,pp)
          ch0<- colSums(P0)>0
          P0<- as.matrix(P0[,ch0])
          G0<- fG(A,eps0); ncG0<- ncol(G0); G0<- matrix(G0[ch0,],ncol=ncG0)
          out0<- mymtcmim(y=y,W=xx,nws=nws,ws=ws,a=out$a,sigma=out$sigma,P=P0,
            G=G0,ngs=ngs0,gs=gs0,b=out$b,init=init,iter=iter,tol=tol)
          if(out0$loglik>lik0){
            lik0<- out0$loglik
            out<- out0
            qtl<- qtl0
            eps<- eps0
            dists<- dist0
          }
        }
        if(k==nQ){
          lik2 = lik1;
          lik1 = lik;
          lik = out$loglik
          if(lik==lik1) break
          la1 = la
          la = lik1+(lik-lik1)/(1-(lik-lik1)/(lik1-lik2))
          if(abs(la)==Inf){
            mx<- Inf
          }else{
            mx = abs(la-la1)
          }
          iter<- iter-1
#          cat(lik); cat('\n')
        }
        if(iter<0){
          warning("\a  mtcmim: convergence failed...")
          break
        }
      }else break
    }
    o<- forder(dists); o<- order(o,decreasing=F)
    dists=dists[o,] #increasingly ordered
  }

  o<- list(loglik=out$loglik,a=out$a,b=out$b,sigma=out$sigma,
        qtl=qtl,eps=eps,dists=dists)
  class(o)<- "mtcmim"
  o
}

###################################
# part II: mtcmim model selection #
# --- general case                #
###################################

# create qtl, dists and b after one QTL added to qtl[[j]]
# NOTE: add to the first place
ffa1<- function(j,obj,ch,mid,d){
  obj0<- obj
  obj0$eps<- NULL
  qtls<- unlist(obj0$qtl)
  nn<- length(qtls)
  np<- nrow(obj$sigma)
  if(nn<1){
    obj0$qtl<- vector("list",np)
    obj0$qtl[[j]]<- 1
    obj0$dists<- data.frame(mid=mid,ch=ch,d=d)
    obj0$b<- 0
  }else{
    obj0$dists<- obj0$dists[order(obj0$dists$ch,obj0$dists$d),]
    obj0$dists<- data.frame(mid=c(mid,obj0$dists$mid),
                             ch=c(ch,obj0$dists$ch),
                              d=c(d,obj0$dists$d))
    idx<- forder(obj0$dists)
    ii<- nn-length(unlist(obj0$qtl[j:np]))
    if(ii>0){
      obj0$b<- c(obj0$b[1:ii],0,obj0$b[-c(1:ii)])
    }else obj0$b<- c(0,obj0$b)
    for(i in 1:np){
      if(i!=j){
        obj0$qtl[[i]]<- idx[-1][obj0$qtl[[i]]]
      }else{
        obj0$qtl[[i]]<- c(idx[1],idx[-1][obj0$qtl[[i]]])
      }
    }
  }

  obj0
}

# create qtl, dists and b after qtl[[i]][j] dropped 
ffd1<- function(i,j,obj){
  obj0<- obj
  obj0$eps<- NULL
  qtls<- unlist(obj0$qtl)
  nn<- length(qtls)

  nn.i<- length(obj0$qtl[[i]])
  if(nn<2){
    obj0$qtl<- NULL
    obj0$dists<- NULL
  }else if(nn.i>0){
    if(nn.i<j || j<1) stop("wrong i or j.")
    obj0$b<- obj0$b[1:nn]
    ii<- j
    if(i>1) ii<- ii+length(unlist(obj0$qtl[1:(i-1)]))
    obj0$b<- obj0$b[-ii]
    obj0$qtl[[i]]<- obj0$qtl[[i]][-j]

    xx<- unique(unlist(obj0$qtl))
      xx<- sort(xx,decreasing=F)
    oo<- order(xx,decreasing=F)
    xxx<- 1:max(xx)
    xxx[xx]<- oo
    for(jj in 1:length(obj0$qtl))
      if(length(obj0$qtl[[jj]])>0)
        obj0$qtl[[jj]]<- xxx[obj0$qtl[[jj]]]
    obj0$dists<- obj0$dists[order(obj0$dists$ch,obj0$dists$d),]
    obj0$dists<- obj0$dists[xx,]
  }

  obj0
}

# rearrange qtl and b so that qtl[[i]] in increasing order
ffo<- function(obj){
  ob<- obj
  np<- length(ob$qtl)
  if(np>0){
    ii<- NULL
    nb<- 0
    for(i in 1:np){
      x<- ob$qtl[[i]]; nx<- length(x)
      if(nx>0){
        o<- order(x,decreasing=F)
        ob$qtl[[i]]<- x[o]
        ii0<- nb+1:nx
        ii<- c(ii,ii0[o])
        nb<- nb+nx
      }
    }
    if(!is.null(ii)) ob$b<- ob$b[ii]
  }
  ob
}

# add one QTL to the model
mtcmim.Add1<- function(object,y,x,xid,mpos,mdat,pp=1,len=1,
  iter=10000,tol=1e-12,ext=FALSE){
# object: object is an object from mtcmim or alike
# xid: xid[[j]] defines which columns of x to be covariates for y_j
  np<- ncol(y)
  bos<- vector("list",np); bliks<- rep(-Inf,np)
  dvv<- div(mpos,len=len)
  if(ext){
     dvv<- dvv[match(unique(dvv$ch),dvv$ch),]
     win<- Inf
  }else win<- 0
  for(j in 1:np){
    loglik<- -Inf
    o0<- NULL
    for(n in 1:nrow(dvv)){
      ob<- ffa1(j,object,ch=dvv$ch[n],mid=dvv$mid[n],d=dvv$d[n])
      o0Tmp<- mtcmim(y,mpos,mdat,x=x,xid=xid,dists=ob$dists,a=ob$a,b=ob$b,sigma=ob$sigma,
        qtl=ob$qtl,eps=ob$eps,win=win,range=-1,pp=pp,len=len,init=1,iter=iter,tol=tol)
      if(o0Tmp$loglik>loglik){
         loglik<- o0Tmp$loglik
         o0<- o0Tmp
      }
    }
    bos[[j]]<- o0
    bliks[j]<- o0$loglik
  }
  tt<- order(bliks,decreasing=TRUE); tt<- tt[1]
  o<- ffo(bos[[tt]])

  add<- FALSE
  if(o$loglik>object$loglik+1e-8) add<- TRUE
  o$add<- add

  o
}

#drop one QTL from the model
mtcmim.Drop1<- function(object,y,x,xid,mpos,mdat,pp=1,len=1,
  iter=10000,tol=1e-12,ext=FALSE){
# object: object from ffg0 or alike
  np<- ncol(y)
  drop<- TRUE
  nn<- 0; for(jj in 1:np) nn<- nn+length(object$qtl[[jj]])
  if(nn<1){
#    cat("\a  no terms to drop...\n")
    drop<- FALSE
    o<- object
  }

  if(ext){
     win<- Inf
  }else win<- 0
  if(drop){
    lik0<- -Inf
    o1<- NULL
    for(i in 1:np){
      nn0<- length(object$qtl[[i]])
      if(nn0>0) for(j in 1:nn0){
        ob<- ffd1(i,j,object)
        o0<- mtcmim(y,mpos,mdat,x=x,xid=xid,dists=ob$dists,a=ob$a,b=ob$b,sigma=ob$sigma,
          qtl=ob$qtl,eps=ob$eps,win=win,range=-1,pp=pp,len=len,init=1,iter=iter,tol=tol)
        if(o0$loglik>lik0+1e-8){
          lik0<- o0$loglik
          o1<- o0
        }
      }
    }
    o<- ffo(o1)
  }
  o$drop<- drop

  o
}

mtcmim.Step<- function(object,y,x,xid,mpos,mdat,cv=0,
  direction=c("both","backward","forward"), pp=1,len=1,
  iter=10000,tol=1e-12,ext=FALSE){
  direction<- match.arg(direction)
  
  o<- object
  if(direction=="both"){
    yes<- TRUE
    if(missing(cv)){
      stop("\a\n  cv: should be positive but missing...\n\n")
    }
    while(yes){
      yes<- FALSE
      od<- mtcmim.Drop1(o,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,pp=pp,len=len,
        iter=iter,tol=tol,ext=ext)
      if(od$drop && 2*(o$loglik-od$loglik)<cv){
        o<- od
        yes<- TRUE
#        cat("-")
      }else{
        oa<- mtcmim.Add1(o,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,pp=pp,len=len,
          iter=iter,tol=tol,ext=ext)
        if(oa$add && 2*(oa$loglik-o$loglik)>cv){
          o<- oa
          yes<- TRUE
#          cat("+")
        }else{
          od<- mtcmim.Drop1(oa,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,pp=pp,len=len,
            iter=iter,tol=tol,ext=ext)
          if(od$drop && od$loglik>o$loglik+1e-8){
             o<- od
             yes<- TRUE
#             cat("=")
          }
        }
      }
    }
  }else if(direction=="backward"){
    if(missing(cv)){
      cv<- Inf
#      cat("\a\n  cv: replaced by Inf...\n\n")
    }
    yes<- TRUE
    while(yes){
      yes<- FALSE
      od<- mtcmim.Drop1(o,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,pp=pp,len=len,
        iter=iter,tol=tol,ext=ext)
      if(od$drop && 2*(o$loglik-od$loglik)<cv){
        o<- od
        yes<- TRUE
#        cat("-")
      }
    }
  }else if(direction=="forward"){
    if(missing(cv)){
      cv<- 0
#      cat("\a\n  cv: replaced by 0...\n\n")
    }
    yes<- TRUE
    while(yes){
      yes<- FALSE
      oa<- mtcmim.Add1(o,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,pp=pp,len=len,iter=iter,tol=tol,ext=ext)
      if(oa$add && 2*(oa$loglik-o$loglik)>cv){
        o<- oa
        yes<- TRUE
#        cat("+")
      }
    }
  }
#  cat("\n")

  o<- list(loglik=o$loglik,a=o$a,b=o$b,sigma=o$sigma,qtl=o$qtl,eps=o$eps,dists=o$dists)
  class(o)<- "mtcmim"
  o
}

####################################
# part III: mtcmim model selection #
# --- special case (zeng 2000)     #
####################################

# create qtl, dists and b after adding one QTL
# NOTE: add to the first place
ffa1All<- function(obj,ch,mid,d){
  obj0<- obj
  obj0$eps<- NULL
  np<- nrow(obj$sigma)
  nn<- length(obj0$qtl[[1]])
  if(nn<1){
    obj0$qtl<- vector("list",np)
    for(i in 1:np) obj0$qtl[[i]]<- 1
    obj0$dists<- data.frame(mid=mid,ch=ch,d=d)
    obj0$b<- rep(0,np)
  }else{
    obj0$dists<- obj0$dists[order(obj0$dists$ch,obj0$dists$d),]
    obj0$dists<- data.frame(mid=c(mid,obj0$dists$mid),
                             ch=c(ch,obj0$dists$ch),
                              d=c(d,obj0$dists$d))
    idx<- forder(obj0$dists)
    obj0$b<- rep(0,(nn+1)*np)
    for(i in 1:np){
      obj0$qtl[[i]]<- c(idx[1],idx[-1][obj0$qtl[[i]]]) #c(1,obj0$qtl[[i]]+1)
    }
  }

  obj0
}

# create qtl, dists and b after dropping the j-th QTL 
ffd1All<- function(j,obj){
  obj0<- obj
  obj0$eps<- NULL
  np<- length(obj0$qtl)
  nn<- length(obj0$qtl[[1]])
  if(nn>1){
    for(k in 1:np) obj0$qtl[[k]]<- 1:(nn-1)
    o<- forder(obj0$dists); oo<- order(o,decreasing=F)
    obj0$dists<- obj0$dists[oo[-j],]
    obj0$b<- rep(0,(nn-1)*np)
  }else{
    obj0$qtl<- NULL
    obj0$dists<- NULL
    obj0$b<- numeric(0)
  }

  obj0
}

# rearrange qtl and b so that qtl[[i]] in increasing order
ffoAll<- function(obj){
  ob<- obj
  np<- length(ob$qtl)
  if(np>0){
    qts<- ob$qtl[[1]]; nn<- length(qts)
    if(nn>0) for(j in 1:np) ob$qtl[[j]]<- 1:nn
  }
  ob
}

# test if qtl[[i]] == qtl[[i]] 
ffequalAll<- function(obj){
  np<- length(obj$qtl)
  if(np>1) for(j in 2:np){
    if(!setequal(obj$qtl[[j]],obj$qtl[[1]])) stop("qtl[[i]] != qtl[[j]]...")
  }
}

#add one QTL to the model
mtcmim.Add1All<- function(object,y,x,xid,mpos,mdat,pp=1,len=1,
  iter=10000,tol=1e-12,ext=FALSE){
# object: object is an object from mtcmim or alike
# xid: xid[[j]] defines which columns of x to be covariates for y_j
  ffequalAll(object)
  np<- ncol(y)
  bos<- vector("list",np); bliks<- rep(-Inf,np)
  dvv<- div(mpos,len=len)
  if(ext){
     dvv<- dvv[match(unique(dvv$ch),dvv$ch),]
     win<- Inf
  }else win<- 0

  loglik<- -Inf
  o0<- NULL
  for(n in 1:nrow(dvv)){
    ob<- ffa1All(object,ch=dvv$ch[n],mid=dvv$mid[n],d=dvv$d[n])
    o0Tmp<- mtcmim(y,mpos,mdat,x=x,xid=xid,dists=ob$dists,a=ob$a,b=ob$b,sigma=ob$sigma,
      qtl=ob$qtl,eps=ob$eps,win=win,range=-1,pp=pp,len=len,init=1,iter=iter,tol=tol)
    if(o0Tmp$loglik>loglik){
      loglik<- o0Tmp$loglik
      o0<- o0Tmp
    }
  }
  o<- ffoAll(o0)

  add<- FALSE
  if(o$loglik>object$loglik+1e-8) add<- TRUE
  o$add<- add

  o
}

# drop one QTL from the model
mtcmim.Drop1All<- function(object,y,x,xid,mpos,mdat,pp=1,len=1,
  iter=10000,tol=1e-12,ext=FALSE){
# object: object from ffg0 or alike
  ffequalAll(object)

  np<- ncol(y)
  drop<- TRUE
  nn<- length(object$qtl[[1]])
  if(nn<1){
#    cat("\a  no terms to drop...\n")
    drop<- FALSE
    o<- object
  }

  if(ext){
     win<- Inf
  }else win<- 0
  if(drop){
    lik0<- -Inf
    o1<- NULL
    for(j in 1:nn){
      ob<- ffd1All(j,object)
      o0<- mtcmim(y,mpos,mdat,x=x,xid=xid,dists=ob$dists,a=ob$a,b=ob$b,sigma=ob$sigma,
        qtl=ob$qtl,eps=ob$eps,win=win,range=-1,pp=pp,len=len,init=1,iter=iter,tol=tol)
      if(o0$loglik>lik0+1e-8){
        lik0<- o0$loglik
        o1<- o0
      }
    }
    o<- ffoAll(o1)
  }
  o$drop<- drop

  o
}

mtcmim.StepAll<- function(object,y,x,xid,mpos,mdat,cv=0,
  direction=c("both","backward","forward"), pp=1,len=1,
  iter=10000,tol=1e-12,ext=FALSE){
  ffequalAll(object)
  direction<- match.arg(direction)
  
  o<- object
  if(direction=="both"){
    yes<- TRUE
    if(missing(cv)){
      stop("\a\n  cv: should be positive but missing...\n\n")
    }
    while(yes){
      yes<- FALSE
      od<- mtcmim.Drop1All(o,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,pp=pp,len=len,
        iter=iter,tol=tol,ext=ext)
      if(od$drop && 2*(o$loglik-od$loglik)<cv){
        o<- od
        yes<- TRUE
#        cat("-")
      }else{
        oa<- mtcmim.Add1All(o,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,pp=pp,len=len,
          iter=iter,tol=tol,ext=ext)
        if(oa$add && 2*(oa$loglik-o$loglik)>cv){
          o<- oa
          yes<- TRUE
#          cat("+")
        }else{
          od<- mtcmim.Drop1All(oa,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,pp=pp,
            len=len,iter=iter,tol=tol,ext=ext)
          if(od$drop && od$loglik>o$loglik+1e-8){
             o<- od
             yes<- TRUE
#             cat("=")
          }
        }
      }
    }
  }else if(direction=="backward"){
    if(missing(cv)){
      cv<- Inf
#      cat("\a\n  cv: replaced by Inf...\n\n")
    }
    yes<- TRUE
    while(yes){
      yes<- FALSE
      od<- mtcmim.Drop1All(o,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,pp=pp,len=len,
        iter=iter,tol=tol,ext=ext)
      if(od$drop && 2*(o$loglik-od$loglik)<cv){
        o<- od
        yes<- TRUE
#        cat("-")
      }
    }
  }else if(direction=="forward"){
    if(missing(cv)){
      cv<- 0
#      cat("\a\n  cv: replaced by 0...\n\n")
    }
    yes<- TRUE
    while(yes){
      yes<- FALSE
      oa<- mtcmim.Add1All(o,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,pp=pp,len=len,
        iter=iter,tol=tol,ext=ext)
      if(oa$add && 2*(oa$loglik-o$loglik)>cv){
        o<- oa
        yes<- TRUE
#        cat("+")
      }
    }
  }
#  cat("\n")

  o<- list(loglik=o$loglik,a=o$a,b=o$b,sigma=o$sigma,qtl=o$qtl,eps=o$eps,dists=o$dists)
  class(o)<- "mtcmim"
  o
}

# unified step function
mtcmimAdd1.default<- function(object,y,x,xid,mpos,mdat,pp=1,len=1,type=1,
  iter=10000,tol=1e-12,ext=FALSE){
  if(type==1){
    mtcmim.Add1(object=object,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,pp=pp,
      len=len,iter=iter,tol=tol,ext=ext)
  }else if(type==2){
    qtl<- unlist(object$qtl)
       qtl<- sort(unique(qtl))
    for(n in 1:length(object$qtl)) object$qtl[[n]]<- qtl
    mtcmim.Add1All(object=object,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,pp=pp,
      len=len,iter=iter,tol=tol,ext=ext)
  }else stop("type 1 or 2 only.")
}

mtcmimDrop1.default<- function(object,y,x,xid,mpos,mdat,pp=1,len=1,type=1,
  iter=10000,tol=1e-12,ext=FALSE){
  if(type==1){
    mtcmim.Drop1(object=object,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,pp=pp,
      len=len,iter=iter,tol=tol,ext=ext)
  }else if(type==2){
    qtl<- unlist(object$qtl)
       qtl<- sort(unique(qtl))
    for(n in 1:length(object$qtl)) object$qtl[[n]]<- qtl
    mtcmim.Drop1All(object=object,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,pp=pp,
      len=len,iter=iter,tol=tol,ext=ext)
  }else stop("type 1 or 2 only.")
}

mtcmimStep.default<- function(object,y,x,xid,mpos,mdat,cv=0,
  direction=c("both","backward","forward"),
  pp=1,len=1,type=1,iter=10000,tol=1e-12,ext=FALSE){
  if(type==1){
    mtcmim.Step(object=object,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,cv=cv,
      direction=direction,pp=pp,len=len,iter=iter,tol=tol,ext=ext)
  }else if(type==2){
    qtl<- unlist(object$qtl)
       qtl<- sort(unique(qtl))
    for(n in 1:length(object$qtl)) object$qtl[[n]]<- qtl
    mtcmim.StepAll(object=object,y=y,x=x,xid=xid,mpos=mpos,mdat=mdat,cv=cv,
      direction=direction,pp=pp,len=len,iter=iter,tol=tol,ext=ext)
  }else stop("type 1 or 2 only.")
}

mtcmimAdd1<- 
   function(object,
            y,
            x,
            xid,
            mpos,
            mdat,
            pp=1,
            len=1,
            type=1,
            iter=10000,
            tol=1e-12,
            ext=FALSE)
{
   UseMethod("mtcmimAdd1")
}

mtcmimDrop1<- 
   function(object,
            y,
            x,
            xid,
            mpos,
            mdat,
            pp=1,
            len=1,
            type=1,
            iter=10000,
            tol=1e-12,
            ext=FALSE)
{
   UseMethod("mtcmimDrop1")
}

mtcmimStep<- 
   function(object,
            y,
            x,
            xid,
            mpos,
            mdat,
            cv=0,
            direction=c("both","backward","forward"),
            pp=1,
            len=1,
            type=1,
            iter=10000,
            tol=1e-12,
            ext=FALSE)
{
   UseMethod("mtcmimStep")
}

################################################################################
# the end #
###########

