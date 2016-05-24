
#####################################
#  functions for variable selection #
################################################################################

# test statistics #
varT2<- function(x,group,equalVar=T){
# Hotelling's T^2
  x<- as.matrix(x)
  grps<- unique(group)
  if(length(grps)!=2) stop("two and only two allowed.")
  group<- (group==grps[1])

  x1<- as.matrix(x[group,])
  n1<- nrow(x1)
  mu1<- colMeans(x1)
  s1<- var(x1)

  x2<- as.matrix(x[!group,])
  n2<- nrow(x2)
  mu2<- colMeans(x2)
  s2<- var(x2)

  if(equalVar){
    t(mu1-mu2)%*%solve((1/n1+1/n2)*((n1-1)*s1+(n2-1)*s2)/(n1+n2-2))%*%(mu1-mu2)
  }else{
    t(mu1-mu2)%*%solve(s1/n1+s2/n2)%*%(mu1-mu2)
  }
}

# permutation
# Note: covariances must equal
perm<- function(x,seed){
  if(!missing(seed)) set.seed(seed)

  x<- as.matrix(x)
  n<- nrow(x)
  x<- as.matrix(x[sample(1:n,n,replace=F),])
}

# variable selection #
varAdd1<- function(x,group,vin=NULL,scope=1:ncol(x),k=0){
# x: data matrix or frame
# group: class indicator (1 or not 1) of the data
# vin: variables already selected
# scope: variables from which to choose one to add
# k: entry
  nr<- dim(x)[1]
  nc<- dim(x)[2]
  if(length(setdiff(vin,1:nc))>0) stop("wrong vin...")
  v<- nr-2
  p<- length(vin)
  vl<- setdiff(scope,vin); if(length(vl)==0) return(NULL)
  f0<- -Inf
  if(is.null(vin)||length(vin)==0){
    t1<- 0
  }else{
    t1<- varT2(x[,vin],group)
  }
  for(i in vl){
    vin0<- sort(union(vin,i))
    t2<- varT2(x[,vin0],group)
    f1<- (v-p)*(t2-t1)/(v+t1)
    if(f1>f0){
      f0<- f1
      vad<- i
    }
  }
  if(f0 >= k+10^(-5)){
    return(vad)
  }else{
    return(NULL)
  }
}

varDrop1<- function(x,group,vin=1:ncol(x),k=0){
# x: data matrix or frame
# group: class indicator (1 or not 1) of the data
# vin: variables already selected
# k: stay
  nr<- dim(x)[1]
  nc<- dim(x)[2]
  if(length(setdiff(vin,1:nc))>0) stop("wrong vin...")
  v<- nr-2
  p<- length(vin)-1
  if(length(vin)==0||is.null(vin)) return(NULL)
  f0<- Inf
  t2<- varT2(x[,vin],group)
  for(i in vin){
    vin0<- sort(setdiff(vin,i))
    if(length(vin0)==0){
      t1<- 0
    }else{
      t1<- varT2(x[,vin0],group)
    }
    f1<- (v-p)*(t2-t1)/(v+t1)
    if(f1<f0){
      f0<- f1
      vdrp<- i
    }
  }
  if(f0 < k-10^(-5)){
    return(vdrp)
  }else{
    return(NULL)
  }
}

varStep<- function(x,group,scope,k,kf=k/2,
  direction=c("both","forward","backward")){
# x: data matrix or frame
# group: class indicator (1 or not 1) of the data
# scope: variables to add or drop
# kf: entry/stay during stepwise forward
# k: entry/stay during stepwise backward
  if(missing(direction)) dir<- "both" else dir<- match.arg(direction)
  nr<- nrow(x)
  nc<- ncol(x)
  if(dir=="backward") vin<- scope else vin<- NULL
  if(dir=="both"){
    while(1){
      vout<- sort(setdiff(scope,vin))
      vad<- varAdd1(x,group,vin,vout,kf)
      if(is.null(vad)){
        break
      }else{
        vin<- sort(union(vin,vad))
      }
      vdrp<- varDrop1(x,group,vin,kf)
      vin<- sort(setdiff(vin,vdrp))
    }
    if(is.null(vin)||length(vin)==0) return(NULL)
    while(1){
      vdrp<- varDrop1(x,group,vin,k)
      if(is.null(vdrp)){
        break
      }else{
        vin<- sort(setdiff(vin,vdrp))
      }
      vout<- sort(setdiff(scope,vin))
      vad<- varAdd1(x,group,vin,vout,k)
      vin<- sort(union(vin,vad))
    }
  }else if(dir=="forward"){
    while(1){
      vout<- sort(setdiff(scope,vin))
      vad<- varAdd1(x,group,vin,vout,kf)
      if(is.null(vad)){
        break
      }else{
        vin<- sort(union(vin,vad))
      }
      vdrp<- varDrop1(x,group,vin,kf)
      vin<- sort(setdiff(vin,vdrp))
    }
    if(is.null(vin)||length(vin)==0) return(NULL)
  }else if(dir=="backward"){
    while(1){
      vdrp<- varDrop1(x,group,vin,k)
      if(is.null(vdrp)){
        break
      }else{
        vin<- sort(setdiff(vin,vdrp))
      }
      vout<- sort(setdiff(scope,vin))
      vad<- varAdd1(x,group,vin,vout,k)
      vin<- sort(union(vin,vad))
    }
  }else stop("wrong direction...")
  vin
}

varGroup<- function(x,z,zscope=NULL,k=qf(0.95,1,nrow(x)-2),kf=k/2,
  method=c("pool","best"),direction=c("both","forward","backward")){
# x: data matrix/frame
# z: variables used for grouping (1 or not)
# zscope: which columns of z to be used for grouping
# k: entry/stay in backward stepwise
# kf: entry/stay in forward stepwise
  if(length(setdiff(zscope,1:ncol(z)))>0) stop("wrond scope...")
  if(is.null(zscope)) zscope<- 1:ncol(z)
  if(missing(direction)) dir<- "both" else dir<- match.arg(direction)
  if(missing(method)){
    method<- "pool"
  }else{
    method<- match.arg(method)
  }
  vl<- 1:ncol(x)
  g<- list()
  idx<- 0
  if(method=="pool"){
    while(!is.null(vl)||length(vl)>0){
      vin<- NULL
      for(j in zscope){
        vin0<- varStep(x,z[,j]==1,vl,k,kf,direction=dir)
        vin<- union(vin,vin0)
      }
      if(length(vin)>0){
        idx<- idx+1
        g$tmp<- sort(vin)
        names(g)[idx]<- paste("group",idx,sep="")
        vl<- sort(setdiff(vl,vin))
      }else break
    }
  }else if(method=="best"){
    while(!is.null(vl)||length(vl)>0){
      vin<- NULL
      for(j in zscope){
        vin0<- varStep(x,z[,j]==1,vl,k,kf,direction=dir)
        if(length(vin0)>length(vin)){
          vin<- vin0
        }
      }
      if(length(vin)>0){
        idx<- idx+1
        g$tmp<- vin
        names(g)[idx]<- paste("group",idx,sep="")
        vl<- sort(setdiff(vl,vin))
      }else break
    }
  }else stop("wrong method...")
  if(length(vl)>0){
    g$tmp<- vl
    names(g)[idx+1]<- "remainder"
  }
  g
}

varSelect<- function(x,group,scope,nv,direction=c("backward","forward")){
# x: data matrix or frame
# group: class indicator (1 or not 1) of the data
# scope: variables to add or drop
  if(missing(scope)) scope<- 1:ncol(x)
  dir<- match.arg(direction)
  if(nv >= length(scope)){
    vin<- scope
    return(vin)
  }
  nr<- nrow(x)
  nc<- ncol(x)

  if(dir=="forward"){
     vin<- NULL
     vin0<- vin
     while(1){
       vout<- sort(setdiff(scope,vin))

       if(length(vin)>nv){
         vdrp<- varDrop1(x,group,vin,Inf)
         vin<- sort(setdiff(vin,vdrp))
       }else{
         vin0<- vin
         vad<- varAdd1(x,group,vin,vout,0)
         vin<- sort(union(vin,vad))
       }

       if(setequal(vin0,vin)) break
     }
   }else if(dir=="backward"){
     vin<- scope
     vin0<- vin
     while(1){
       vout<- sort(setdiff(scope,vin))

       if(length(vin)>=nv){
         vin0<- vin
         vdrp<- varDrop1(x,group,vin,Inf)
         vin<- sort(setdiff(vin,vdrp))
       }else{
         vad<- varAdd1(x,group,vin,vout,0)
         vin<- sort(union(vin,vad))
       }

       if(setequal(vin0,vin)) break
     }
   }else stop("direction specification wrong.")

  vin
}

################################################################################
# the end #
###########

