
############################
#  miscellaneous functions #
################################################################################

#####################
# generic functions #
#####################
haldane<-  function(c1){
# c to d (M)
  -0.5*log(1-2*c1)
}

haldane.inv<- function(d){
# d (M) to c
  (1-exp(-2*d))/2
}

recomb.rate.ac<- function(r.ab,r.bc){
# calculate recombination rate ("a-b-c")
# no interference assumed
  r.ab+r.bc-2*r.ab*r.bc
}

recomb.rate.bc<- function(r.ab,r.ac){
  (r.ac-r.ab)/(1-2*r.ab)
}

#pp: population, 1=BC, 2=RIL-selfing, 3=RIL-brother-sister-mating
getp<- function(m1,m2,d,pp=1){
# P(m2|m1)
# d: cM, dist between m1 and m2
  c<- haldane.inv(d/100)
  if(pp==1){
    cc<- c
  }else if(pp==2){
    cc<- 2*c/(1+2*c)
  }else if(pp==3){
    cc<- 4*c/(1+6*c)
  }else stop("\a  Only BC-1 or RIL (selfing-2, brother-sister mating-3) considered...")
  if(m1==m2) cc<- 1-cc
  cc
}

getprob<- function(m1,m2,r1,r){
# cond. prob. of QTL(QQ) given flanking markers
# r:m1-m2 and r1:m1-qtl
  r2<- recomb.rate.bc(r1,r)
  (m1-r1)*(m2-r2)/((m1-1)*m2+m1*(m2-1)+1-r)
}

mtoi<- function(m,nmark){
# m:vector of marker positions
# nmark:vector of numbers of intervals on each chromosome
# return interval positions
  if(is.null(nmark))return(NULL)
  I<- NULL
  for(k in m){
    if(is.element(k,cumsum(nmark))){
      stop("\a  wrong information!")
    }else {
      ind.tmp<- (cumsum(nmark)<k)
      I<- c(I,sum(nmark[ind.tmp]-1)+k-sum(nmark[ind.tmp]))
    }
  }
  I
}

itom<- function(i,nmark){
  if(is.null(nmark))return(NULL)
  I<- NULL
  for(k in i){
    if(!is.na(k)){
      ind.tmp<- (cumsum(nmark-1)<k)
      I<- c(I,sum(nmark[ind.tmp])+k-sum(nmark[ind.tmp]-1))
    }
  }
  I
}

div<- function(mpos,len=1){
# len: steplength (cM)
  mid<- NULL
  ch<- NULL
  dist<- NULL

  N<- dim(mpos)[1]
  for(n in 1:(N-1)){
    tmp1<- mpos$dist[n]
    tmp2<- mpos$dist[n+1]
    dff<- tmp2-tmp1
    tmp<- tmp1
    if(dff>0){
      len0<- max(1,ceiling(dff/len))
      len0<- dff/len0
      while(tmp+len0/10 < tmp2){
        mid<- c(mid,mpos$id[n])
        ch<- c(ch,mpos$ch[n])
        dist<- c(dist,tmp)
        tmp<- tmp+len0
      }
    }else{
      mid<- c(mid,mpos$id[n])
      ch<- c(ch,mpos$ch[n])
      dist<- c(dist,mpos$dist[n])
    }
  }
  mid<- c(mid,mpos$id[N])
  ch<- c(ch,mpos$ch[N])
  dist<- c(dist,mpos$dist[N])

  data.frame(mid,ch,dist)
}

#########################
# create mpos and dists #
#########################
gv2mpos<- function(gmap,v){
# v: g$v where g is a SUR object
  mpos<- gmap
  if(!is.null(mpos$ch) && !is.numeric(mpos$ch))
    stop("gmap: chromosome IDs should be integers.")
  mpos<- gmap[order(mpos$ch,mpos$dist),]
  chrs<- unique(mpos$ch)
  mid<- NULL
  for(n in 1:length(chrs))
    mid<- c(mid,1:sum(mpos$ch==chrs[n]))
  mpos<- data.frame(id=1:nrow(mpos),ch=mpos$ch,m=mid,dist=mpos$dist)

  idx<- NULL
  for(j in 1:length(v)) idx<- c(idx,v[[j]])
    idx<- sort(unique(idx))
  dists<- data.frame(ch=mpos$ch[idx],mid=mpos$id[idx],d=mpos$dist[idx])

  list(mpos=mpos,dists=dists)
}

###################
# xid: covariates #
###################
xid1ch<- function(mpos,v,k){
# v: g$v where g is a SUR object
  np<- length(v)
  mid<- mpos$id[mpos$ch==k]
  xid<- vector("list",np)
  for(j in 1:np){
    xid[[j]]<- sort(setdiff(v[[j]],mid))
  }
  xid
}

#############################
# extract covariate effects #
#############################
#a: object$a
xeff<- function(a,xid){
  eff<- vector("list",length(xid))
  ii<- 0
  for(k in 1:length(xid)){
    x<- xid[[k]]; nx<- length(x)+1 # intercept not included in x
    eff[[k]]<- a[(ii+1):(ii+nx)]
    ii<- ii+nx
  }

  eff
}

#######################
# extract QTL effects #
#######################
#object: mtcmim object
qeff<- function(object){
  eff<- vector("list",length(object$qtl))
  ii<- 0
  for(k in 1:length(object$qtl)){
    x<- object$qtl[[k]]; nx<- length(x)
    if(nx>0) eff[[k]]<- object$b[(ii+1):(ii+nx)]
    ii<- ii+nx
  }

  list(main=eff,eps=object$b[-c(1:ii)])
}

################################################################################
# the end #
###########

