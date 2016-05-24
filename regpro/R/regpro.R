additive2<-function(x,y,arg,h=1,kernel="gauss",M=2)
{
d<-length(arg)
n<-length(y)

if (kernel=="gauss") ker<-function(t){ return( exp(-t^2/2) ) }
if (kernel=="uniform") ker<-function(t){ return( (abs(t) <= 1) ) }
if (kernel=="bart") ker<-function(t){ return( (1-t) ) }

hatc<-mean(y)
residual<-matrix(y-hatc,n,1)

for (m in 1:M){
   for (j in 1:d){
      colu<-matrix(x[,j],n,1)
      estim<-matrix(0,n,d)
      if (j==1) rest<-seq(2,d) else if (j==d) rest<-seq(1:(d-1)) 
      else rest<-c(seq(1:(j-1)),seq((j+1):d))
      for (l in rest){
          for (nn in 1:n){
             curarg<-x[nn,l] 
             w<-kernesti.weights(curarg,colu,h=h,kernel=kernel)
             estim[nn,l]<-t(w)%*%residual
          }
       }
       residual<-y-hatc-matrix(rowSums(estim),n,1)
   }
}

valuevec<-matrix(0,d,1)
for (i in 1:d){
    curx<-matrix(x[,i],n,1)
    w<-kernesti.weights(arg[i],curx,h=h,kernel=kernel)
    valuevec[i]<-t(w)%*%residual
}

return(hatc+valuevec)
}


additive3<-function(x,y,arg,h=1,kernel="gauss",M=2)
{
d<-length(arg)
n<-length(y)

hatc<-mean(y)
estim<-matrix(0,n,d)
for (m in 1:M){
   itestim<-estim
   for (j in 1:d){
      colu<-matrix(x[,j],n,1)
      for (nn in 1:n){
          curarg<-x[nn,j]
          jestim<-itestim
          jestim[,j]<-0 
          ycur<-y-hatc-matrix(rowSums(jestim),n,1)
          w<-kernesti.weights(curarg,colu,h=h,kernel=kernel)
          itestim[nn,j]<-t(w)%*%ycur
      }
      itestim[,j]<-itestim[,j]-mean(itestim[,j])      
   }
   estim<-itestim
}

valuevec<-matrix(0,d,1)
for (j in 1:d){
    curx<-matrix(x[,j],n,1)
    w<-kernesti.weights(arg[j],curx,h=h,kernel=kernel)
    jestim<-estim
    jestim[,j]<-0 
    ycur<-y-hatc-matrix(rowSums(jestim),n,1)
    valuevec[j]<-t(w)%*%ycur
}

return(estim)
#return(hatc+valuevec)
}



additive.old<-function(x,y,arg,h=1,kernel="gauss",M=2)
{
d<-length(arg)
n<-length(y)

if (kernel=="gauss") ker<-function(t){ return( exp(-t^2/2) ) }
if (kernel=="uniform") ker<-function(t){ return( (abs(t) <= 1) ) }
if (kernel=="bart") ker<-function(t){ return( (1-t) ) }

G<-matrix(0,n,d)    # estimators g_j evaluated at x_i^j 
hatc<-mean(y)
residual<-matrix(y-hatc,n,1)

for (m in 1:M){
   for (j in 1:d){
      colu<-x[,j]
      pairdiffe<-matrix(colu,n,n,byrow=FALSE)-matrix(colu,n,n,byrow=TRUE)
      Wj<-ker(pairdiffe)
      Wj<-Wj/colSums(Wj)
      G[,j]<-t(Wj)%*%residual
      residual<-y-hatc-matrix(rowSums(G),n,1)
   }
}

argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
W<-ker((x-argu)/h)/h # W<-matrix(0,n,d) kernel weights 
W<-W/colSums(W)
valuevec<-t(W)%*%residual

return(valuevec)
}




additive<-function(x,y,arg=NULL,eval=NULL,
h=1,kernel="gauss",M=2,vect=FALSE)
{
d<-dim(x)[2]
n<-length(y)
hatc<-mean(y)

if (!vect){

if (is.null(eval)){
 eval<-matrix(0,n,d)
 for (m in 1:M){
   for (j in 1:d){
      colu<-matrix(x[,j],n,1)
      jeval<-eval
      jeval[,j]<-0 
      ycur<-y-hatc-matrix(rowSums(jeval),n,1)
      for (nn in 1:n){
          curarg<-x[nn,j]
          w<-kernesti.weights(curarg,colu,h=h,kernel=kernel)
          eval[nn,j]<-t(w)%*%ycur
      }
      eval[,j]<-eval[,j]-mean(eval[,j])      
   }
 }
}

if (is.null(arg)){ 
  value<-NULL
  valuevec<-NULL
}
else{
  valuevec<-matrix(0,d,1)
  for (j in 1:d){
     curx<-matrix(x[,j],n,1)
     w<-kernesti.weights(arg[j],curx,h=h,kernel=kernel)
     jeval<-eval
     jeval[,j]<-0 
     ycur<-y-hatc-matrix(rowSums(jeval),n,1)
     #ycur<-matrix(eval[,j],n,1)
     valuevec[j]<-t(w)%*%ycur
  }
  value<-sum(valuevec)+hatc
}
}

######################################################
if (vect){

if (is.null(eval)){
 eval<-matrix(0,n,d)
 for (m in 1:M){
   for (j in 1:d){
      colu<-matrix(x[,j],n,1)
      jeval<-eval
      jeval[,j]<-0 
      ycur<-y-hatc-matrix(rowSums(jeval),n,1)
      #############################
      xarg<-matrix(x[,j],n,1)
      W<-kernesti.weights(xarg,colu,h=h,kernel=kernel,vect=TRUE) 
      eval[,j]<-t(W)%*%ycur
      #############################
      eval[,j]<-eval[,j]-mean(eval[,j])      
   }
 }
}
if (is.null(arg)){ 
  value<-NULL
  valuevec<-NULL
}
else{
  valuevec<-matrix(0,d,1)
  for (j in 1:d){
     curx<-matrix(x[,j],n,1)
     w<-kernesti.weights(arg[j],curx,h=h,kernel=kernel)
     jeval<-eval
     jeval[,j]<-0 
     ycur<-y-hatc-matrix(rowSums(jeval),n,1)
     #ycur<-matrix(eval[,j],n,1)
     valuevec[j]<-t(w)%*%ycur
  }
  value<-sum(valuevec)+hatc
}
}

return(list(eval=eval,value=value,valvec=valuevec))
}


additive.stage2<-function(arg,x,y,h=1,kernel="gauss",B=2)
{
d<-length(arg)

if (kernel=="gauss") ker<-function(t){ return( exp(-t^2/2) ) }
if (kernel=="uniform") ker<-function(t){ return( (abs(t) <= 1) ) }

ynow<-y
resu<-0
ssr<-matrix(0,d,1)
for (ii in 1:B){
    for (jj in 1:d){
        xjj<-x[,jj]
        xjj<-matrix(xjj,length(xjj),1)
        argu<-matrix(arg[jj],length(xjj),1)
        w<-ker((xjj-argu)/h)/h^d
        w<-w/sum(w)
        yhat<-w%*%ynow
        ssr[jj]<-sum((yhat-ynow)^2)
    }  
    dstar<-which.min(ssr)
    arg.now<-arg[dstar]
    x.now<-x[,dstar]

    w<-kernesti.weights(arg.now,x.now,h=h,kernel=kernel)
    w<-matrix(w,length(w),1)    
    ynow<-matrix(ynow,length(ynow),1)
    notna<-(!is.na(ynow))
    w<-notna*w
    w<-w/sum(w)
    mu<-w*ynow
    neweva<-sum(mu,na.rm=TRUE)

    resu<-resu+neweva    
    residu<-1 
    ynow<-residu        
}

return(resu)
}







additive.stage<-function(x,y=NULL,arg=NULL,residu=NULL,deet=NULL,
h=1,kernel="gauss",M=2,vect=FALSE)
{
d<-dim(x)[2]
n<-dim(x)[1]

if (is.null(residu)){
  residu<-matrix(0,n,M)
  deet<-matrix(0,M,1)
  estim<-matrix(0,n,1)
  eval<-matrix(0,n,d)
  for (m in 1:M){
    residu[,m]<-y-estim
    ssr<-matrix(0,d,1)
    estimat<-matrix(0,n,d)
    for (j in 1:d){
      colu<-matrix(x[,j],n,1)
      ycur<-residu[,m]
      if (!vect){
         for (nn in 1:n){
             curarg<-x[nn,j]
             w<-kernesti.weights(curarg,colu,h=h,kernel=kernel)
             estimat[nn,j]<-t(w)%*%ycur
         }
      }
      else{
         estimat[,j]<-kernesti.regr(colu,colu,ycur,h=h,kernel=kernel,vect=vect)  
      }
      ssr[j]<-sum((ycur-estimat[,j])^2)
    }
    dstar<-which.min(ssr)
    deet[m]<-dstar
    eval[,dstar]<-eval[,dstar]+estimat[,dstar]
    estim<-estim+estimat[,dstar]
  }
}
else eval<-NULL

if (is.null(arg)){ 
  value<-NULL
  valuevec<-NULL
}
else{
  valuevec<-matrix(0,d,1)
  for (m in 1:M){
     ycur<-matrix(residu[,m],n,1)
     xcur<-matrix(x[,deet[m]],n,1)
     w<-kernesti.weights(arg[deet[m]],xcur,h=h,kernel=kernel)
     curvalue<-t(w)%*%ycur
     valuevec[deet[m]]<-valuevec[deet[m]]+curvalue
  }
  value<-sum(valuevec)
}

return(list(eval=eval,residu=residu,deet=deet,value=value,valvec=valuevec))
}



bagg<-function(x,y,arg,B=10,seed=1,method="worpl",propor=0.5,
estimator="greedy",M=2,m=5,splitfreq=1)
{
n<-length(y)
d<-dim(x)[2]
vals<-matrix(0,B,1)
bootstrap<-function(n,seed,method.propor){1}

if (estimator=="greedy")
for (i in 1:B){
   seed<-seed+i
   boot<-bootstrap(n,seed=seed,method=method,propor=propor)
   if (d==1) xsub<-matrix(x[boot],length(boot),1) else xsub<-x[boot,]
   ysub<-matrix(y[boot],length(boot),1)
   vals[i]<-greedy(xsub,ysub,arg,M,m,splitfreq=splitfreq)$val
}
else  # estimator=="linear"
for (i in 1:B){
   seed<-seed+i
   boot<-bootstrap(n,seed=seed,method=method,propor=propor)
   if (d==1) xsub<-matrix(x[boot],length(boot),1) else xsub<-x[boot,]
   ysub<-matrix(y[boot],length(boot),1)
   lin<-linear(xsub,ysub)
   vals[i]<-lin$beta0+lin$beta1%*%arg
}

val<-mean(vals)
return(val)
}










copula.trans<-function(dendat,marginal=rep("gauss",dim(dendat)[2]),remna=TRUE)
{

n<-dim(dendat)[1]
d<-dim(dendat)[2]
copdat<-dendat

for (ii in 1:d){
   if ((marginal[ii]=="gauss")||(marginal[ii]=="uniform")){
      or<-order(dendat[,ii])
      mones<-matrix(0,n,1)
      for (i in 1:n) mones[or[i]]<-i  
      copdat[,ii]<-mones/(n+1)  # copdat[or,ii]<-seq(1:n)/(n+1)
   }
   if (marginal[ii]=="gauss"){ 
      copdat[,ii]<-qnorm(copdat[,ii])
      for (jj in 1:d){
           copdat[(copdat[,jj]==Inf),jj]<-NA
           copdat[(copdat[,jj]==-Inf),jj]<-NA
      }
   }
}

# we remove those rows where there is at least one NA
if (remna){
   for (ii in 1:d){
      indexes<-!is.na(copdat[,ii])
      copdat<-matrix(copdat[indexes,],sum(indexes),d)
   }
}
return(copdat)
}


digit<-function(luku,base){
#Gives representation of luku for system with base
#
#luku is a natural number >=0
#base is d-vector of integers >=2, d>=2, 
#base[d] tarvitaan vain tarkistamaan onko luku rajoissa
#
#Returns d-vector of integers.
#
#example: digit(52,c(10,10)), returns vector (2,5)
#
d<-length(base)
digi<-matrix(0,d,1)
jako<-matrix(0,d,1)
jako[d]<-base[1]
for (i in (d-1):1){
  jako[i]<-base[d-i+1]*jako[i+1]
}
vah<-0
for (i in 1:(d-1)){
  digi[i]<-floor((luku-vah)/jako[i+1]) #if digi[i]>base[i], then ERROR
  vah<-vah+digi[i]*jako[i+1]
}
digi[d]<-luku-vah
# annetaan vastaus kaanteisesti se 2354 annetaan c(4,5,3,2)
# talloin vastaavuus sailyy base:n kanssa 
#apu<-matrix(0,d,1)
#for (i in 1:d){
#  apu[i]<-digi[d-i+1]
#}
apu<-digi[d:1]
return(apu)
}
emp.distribu<-function(arg,dendat)
{
d<-length(arg)
n<-dim(dendat)[1]

if (d>1){
   val<-matrix(0,d,1)
   for (kk in 1:d){
        cateca<-c(dendat[,kk],arg[kk])
        or<-rank(cateca)
        val[kk]<-or[n+1]/(n+2)
   }
}
else{
        cateca<-c(dendat,arg)
        or<-rank(cateca)
        val<-or[n+1]/(n+2) 
}

return(val)
}


emp.quantile<-function(arg,dendat)
{
d<-length(arg)

if (d>1){
   val<-matrix(0,d,1)
   n<-dim(dendat)[1] 
   for (kk in 1:d){
        or<-rank(dendat[,kk])
        uni<-or/(n+1)
        cateca<-c(uni,arg[kk])
        or2<-rank(cateca)
        inter<-or2[length(cateca)]
        if (inter==(n+1)) inter<-n
        so<-sort(dendat[,kk])
        val[kk]<-so[inter]

   }   
}
else{
        n<-length(dendat) 
        or<-rank(dendat)
        uni<-or/(n+1)
        cateca<-c(uni,arg)
        or2<-rank(cateca)
        inter<-or2[length(cateca)]
        if (inter==(n+1)) inter<-n
        val<-sort(dendat)[inter]
}

return(val)
}





greedy<-function(x,y,arg,M,m,splitfreq=1)
{
d<-length(arg)
n<-length(y)
#s<-splitsearch(x,y,arg,splitfreq=splitfreq)

inx<-matrix(0,n*d+1,1)
for (i in 1:n){
    for (j in 1:d){
        inx[1+(i-1)*d+j]=x[i,j]
    }
}
iny<-matrix(0,n+1,1)
iny[2:(n+1)]<-y
inarg<-matrix(0,d+1,1)
inarg[2:(d+1)]<-arg
kg<-1
#kg<-.C("splitSearch",
#           as.double(inx),
#           as.double(iny),
#           as.double(inarg),
#           as.double(splitfreq),
#           as.integer(n),
#           as.integer(d),
#           indeces = integer(n+1),
#           lkm = integer(1)
#)
s<-kg$indeces[2:(kg$lkm+1)]

xnew<-matrix(x[s,],length(s),d)
ynew<-matrix(y[s],length(s),1)
num<-length(s)
lkm<-2
while ((lkm<=M) && (num>m)){
    #s<-splitsearch(xnew,ynew,arg,splitfreq=splitfreq)    

    n<-length(ynew)
    inx<-matrix(0,n*d+1,1)
    for (i in 1:n){
       for (j in 1:d){
          inx[1+(i-1)*d+j]=xnew[i,j]
       }
    }
    iny<-matrix(0,n+1,1)
    iny[2:(n+1)]<-ynew
    inarg<-matrix(0,d+1,1)
    inarg[2:(d+1)]<-arg
    #kg<-.C("splitSearch",
    #           as.double(inx),
    #           as.double(iny),
    #           as.double(inarg),
    #           as.double(splitfreq),
    #           as.integer(n),
    #           as.integer(d),
    #           indeces = integer(n+1),
    #           lkm = integer(1)
    #)
    s<-kg$indeces[2:(kg$lkm+1)]

    num<-length(s)
    if (num>=m){
       xnew<-matrix(xnew[s,],length(s),d)
       ynew<-matrix(ynew[s],length(s),1)
    }
    lkm<-lkm+1
}

val<-mean(ynew)
return(list(val=val,x=xnew,y=ynew))
}


kernesti.der<-function(arg,x,y,h=1,direc=1,kernel="gauss",vect=FALSE)
{
d<-dim(x)[2]

if (d>1){

n<-dim(x)[1]

if (kernel=="gauss"){
   ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
   dker<-function(xx){ 
         return( -(2*pi)^(-d/2)*xx[,direc]*exp(-rowSums(xx^2)/2) ) }
}

argu<-matrix(arg,n,d,byrow=TRUE)
we<-ker((argu-x)/h)/h^d
w<-we/sum(we)
u<-dker((argu-x)/h)/h^(d+1)
q<-1/sum(we)*(u-w*sum(u))  
value<-q%*%y 
return(value)
}

if (d==1){

if (kernel=="gauss"){
   ker<-function(xx){ return( exp(-xx^2/2) ) }
   dker<-function(xx){ return( -xx*exp(-xx^2/2) ) }
   dker2<-function(t){ return( -t*(2*pi)^(-1/2)*exp(-t^2/2) ) }
}

if (!vect){
 w<-ker((arg-x)/h)/h^1
 we<-w/sum(w)
 u<-dker((arg-x)/h)/h^(1+1)
 q<-1/sum(w)*(u-we*sum(u))  
 value<-sum(y*q)  #y%*%q
 return(value)
}
if (vect){
  n<-length(x)
  x<-matrix(x,length(x),1)
  arg<-matrix(x,length(arg),1) 
  xu<-matrix(x,n,n)
  argu<-matrix(arg,n,n)
  w<-ker((argu-xu)/h)/h^1
  we<-t(t(w)/colSums(w))
  u<-dker((argu-xu)/h)/h^(1+1)

  
  q<-1/sum(w)*(u-we*sum(u))  


  y<-matrix(y,1,n)
  value<-y%*%q
  return(value)
}

}

}



kernesti.quantile<-function(arg,x,y,h=1,p=0.5,kernel="gauss")
{
n<-length(y)
lkm<-length(p)
quan<-matrix(0,lkm,1)

w<-kernesti.weights(arg,x,h=h,kernel=kernel)
or<-order(y)
weet<-w[or]
i<-1
zum<-0
for (ii in 1:lkm){
    pcur<-p[ii]
    while ( (i<=n) && (zum<pcur) ){
        zum<-zum+weet[i]
        i<-i+1
    }
    if (i>n) quan[ii]<-max(y) else quan[ii]<-y[or[i]]
}

return(quan)
}




kernesti.regr<-function(arg,x,y,h=1,kernel="gauss",g=NULL,gernel="gauss",
vect=FALSE)
{
w<-kernesti.weights(arg,x,h=h,kernel=kernel,vect=vect,g=g,gernel=gernel)
y<-matrix(y,length(y),1)

if (!vect){
  notna<-(!is.na(y))
  w<-matrix(w,length(w),1)    
  w<-notna*w
  w<-w/sum(w)
  mu<-w*y
  est<-sum(mu,na.rm=TRUE)
}
else{
  est<-t(w)%*%y  #t(y)%*%w
}

return(est)
}



kernesti.weights<-function(arg,x,h=1,kernel="gauss",g=NULL,gernel="gauss",
vect=FALSE)
{
if (!vect) d<-length(arg) else d<-dim(x)[2]

if (d>1){

n<-dim(x)[1]

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ ans<-(rowSums(xx^2) <= 1) 
                      return( ans ) }
if (kernel=="exponential") ker<-function(xx){ return( exp(-rowSums(xx)) ) }

argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
w<-ker((x-argu)/h)/h^d

if ((n==1)&&(is.na(w[1]))) w[1]<-1
w[is.na(w)]<-0   # make NA:s to zero

if (sum(w)==0) weights<-rep(1,n)/n else weights<-w/sum(w)

if (!is.null(g)){

   if (gernel=="bart") 
   ger<-function(xx){ return( (1-rowSums(xx^2))*(rowSums(xx^2)<= 1) ) }
   if (gernel=="gauss") 
   ger<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
   if (gernel=="uniform") 
   ger<-function(xx){ ans<-(rowSums(xx^2)<= 1) 
                      return( ans ) }

   argui<-matrix(seq(n,1,-1),n,1)
   w<-ker((x-argu)/h)/h^d*ger((n-argui)/g)/g
   weights<-w/sum(w)
}
}
else{  # d==1  #########################################

n<-length(x)

if (kernel=="gauss") ker<-function(xx){ return( exp(-xx^2/2) ) }
if (kernel=="uniform") ker<-function(xx){ return( (abs(xx) <= 1) ) }
if (kernel=="uniform01") ker<-function(xx){ return( (abs(2*xx-1) <= 1) ) }
if (kernel=="exp") ker<-function(xx){ return( (xx>=0)*exp(-xx) ) }

if (!vect){
  x<-matrix(x,length(x),1)
  argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
  w<-ker((argu-x)/h)/h^d
  if ((n==1)&&(is.na(w[1]))) w[1]<-1
  w[is.na(w)]<-0   # make NA:s to zero
  if (sum(w)==0) weights<-rep(1,n)/n else weights<-w/sum(w)
}

if (vect){
  n<-length(x)
  x<-matrix(x,length(x),1)
  arg<-matrix(x,length(arg),1) 
  xu<-matrix(x,n,n)
  argu<-matrix(arg,n,n)
  argu<-t(argu)
  w<-ker((argu-xu)/h)/h
  w[is.na(w)]<-0   # make NA:s to zero
  #weights<-w%*%diag(1/colSums(w))
  weights<-t(t(w)/colSums(w))
}

if (!is.null(g)){

   n<-length(x)
   if (gernel=="bart") 
   ger<-function(xx){ return( (1-rowSums(xx^2))*(rowSums(xx^2)<= 1) ) }
   if (gernel=="gauss") 
   ger<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
   if (gernel=="uniform") 
   ger<-function(xx){ ans<-(rowSums(xx^2)<= 1) 
                      return( ans ) }
   if (gernel=="exp") ger<-function(xx){ return( exp(-rowSums(xx))*(xx>=0) ) }

   argui<-matrix(seq(1,n,1),n,1)   # matrix(seq(n,1,-1),n,1)
   w<-ker((x-arg)/h)/h^d*ger((n-argui)/g)/g
   weights<-w/sum(w)
}

}

return(weights=weights)
}



linear.quan<-function(x,y,p=0.5)
{
y<-matrix(y,length(y),1)
n<-dim(x)[1]
d<-dim(x)[2]

rho<-function(t){
    t*(p-(t<0))
}

fn<-function(b) { 
    b2<-matrix(b[2:(d+1)],d,1)
    gx<-b[1]+x%*%b2
    ro<-rho(y-gx)
    return(sum(ro)/n)
}

li<-linear(x,y)
par<-c(li$beta0,li$beta1)    # initial value
par.lower<-rep(-1,d)
par.upper<-rep(1,d)

op.method<-"L-BFGS-B"
op<-optim(par=par,fn=fn,method=op.method)
theta<-op$par

beta0<-theta[1]
beta1<-theta[2:(d+1)]

return(list(beta0=beta0,beta1=beta1))
}


linear<-function(x,y,eleg=TRUE,lambda=0)
{
y<-matrix(y,length(y),1)
n<-dim(x)[1]
d<-dim(x)[2]

if (d==1){
  beta1<-cov(x,y)/var(x)
  beta0<-mean(y)-beta1*mean(x)
  return(list(beta0=beta0,beta1=beta1))
}
else{
  if (lambda==0){
  if (eleg){
  xnew<-matrix(0,n,d+1)
  xnew[,1]<-1
  xnew[,2:(d+1)]<-x
  alku<-t(xnew)%*%xnew
  invalku<-solve(alku,diag(rep(1,d+1)))
  beta<-invalku%*%t(xnew)%*%y
  return(list(beta0=beta[1],beta1=beta[2:(d+1)]))
  }
  else{
  kov<-cov(x)
  inkov<-solve(kov,diag(rep(1,d))) 
  beta<-inkov%*%cov(x,y)
  beta0<-mean(y)-t(beta)%*%t(t(colMeans(x)))
  return(list(beta0=beta0,beta1=beta))
  }
  }
  if (lambda>0){
  Y<-y-mean(y)
  #X<-(x-colMeans(x))/sqrt(colMeans(x^2)-colMeans(x)^2) 
  X<-x-colMeans(x)
  X<-X/sqrt(colMeans(X^2))
  XtX<-t(X)%*%X+lambda*diag(rep(1,d))
  invXtX<-solve(XtX,diag(rep(1,d)))
  beta<-invXtX%*%t(X)%*%Y
  return(list(beta0=mean(y),beta1=beta))
  }
}
}


loclin<-function(arg,x,y,h=1,kernel="gauss",type=0)
{
d<-length(arg)
n<-length(y)

if (d>1){

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( (2*pi)^(-d/2)*exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ ans<-(rowSums(xx^2) <= 1) 
                      return( ans ) }

argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
w<-ker((x-argu)/h)/h^d
weights<-w/sum(w)

X<-cbind(matrix(1,n,1),x-argu)
W<-diag(weights)
A<-t(X)%*%W%*%X     
invA<-solve(A,diag(rep(1,d+1))) 
B<-t(X)%*%W%*%y
esti<-invA%*%B
est<-esti[type+1]

}
else{  # d==1  #########################################

if (kernel=="gauss") ker<-function(xx){ return( exp(-xx^2/2) ) }
if (kernel=="uniform") ker<-function(xx){ return( (abs(xx) <= 1) ) }

x<-matrix(x,length(x),1)
w<-ker((x-arg)/h)/h^d   
weights<-w/sum(w)

X<-cbind(matrix(1,n,1),x-arg)
W<-diag(c(weights))
A<-t(X)%*%W%*%X     
invA<-solve(A,diag(rep(1,d+1))) 
B<-t(X)%*%W%*%y
esti<-invA%*%B
est<-esti[type+1] 

other<-FALSE
if (other){
w<-ker((arg-x)/h); p<-w/sum(w)
barx<-sum(p*x); bary<-sum(p*y)
q<-p*(1-((x-barx)*(barx-arg))/sum(p*(x-barx)^2))


s1<-sum(w*(x-arg))
s2<-sum(w*(x-arg)^2)
q<-w*(s2-(x-arg)*s1)/sum(w*(s2-(x-arg)*s1))
}

}

return(est)
}

ma<-function(x,h=1,kernel="exp",k=length(x))
{
if (kernel=="exp") 
   ker<-function(xx){ return( exp(-xx) ) }
if (kernel=="bart") 
   ker<-function(xx){ return( 1-xx^2 ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( exp(-xx^2/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ ans<-(xx^2 <= 1) 
                      return( ans ) }

t<-seq(0,k-1,1)
w<-ker(t/h)
w<-w/sum(w)
w<-w[length(w):1]
ans<-sum(w*x)

return(ans)
}


pcf.additive<-function(x,y,h,N,kernel="gauss",support=NULL,
M=2,eval=NULL,direc=NULL)
{
d<-length(N)
n<-length(y) 
hatc<-mean(y)

if (is.null(eval)){
 eval<-matrix(0,n,d)
 for (m in 1:M){
   for (j in 1:d){
      colu<-matrix(x[,j],n,1)
      jeval<-eval
      jeval[,j]<-0 
      ycur<-y-hatc-matrix(rowSums(jeval),n,1)
      for (nn in 1:n){
          curarg<-x[nn,j]
          w<-kernesti.weights(curarg,colu,h=h,kernel=kernel)
          eval[nn,j]<-t(w)%*%ycur
      }
      eval[,j]<-eval[,j]-mean(eval[,j])      
   }
 }
}

if (is.null(direc)){

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ return( (rowSums(xx^2) <= 1) ) } 

recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

support<-matrix(0,2*d,1)
for (i in 1:d){
    support[2*i-1]<-min(x[,i])
    support[2*i]<-max(x[,i])
}
lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2

     valli<-additive(x,y,arg,h=h,M=M,eval=eval)$value 
     value[i]<-valli
     index[i,]<-inde
}

down<-index-1
high<-index

pcf<-list(
value=value,index=index,
down=down,high=high,  
support=support,N=N)

}

if (!is.null(direc)){  

y<-matrix(y,1,length(y))
x<-matrix(x,length(x),1)

if (kernel=="gauss") ker<-function(xx){ return( exp(-xx^2/2) ) }
if (kernel=="uniform") ker<-function(xx){ return( (abs(xx) <= 1) ) }

index<-seq(1:N)
len<-length(index)
down<-matrix(0,len,1)
high<-matrix(0,len,1)
down[,1]<-index-1
high[,1]<-index

value<-matrix(0,N,1)
if (is.null(support)){
   support<-matrix(0,2,1)
   support[1]<-min(x)
   support[2]<-max(x)
}
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

for (i in 1:N){
     inde<-i
     argu<-lowsuppo+step*inde-step/2
     arg<-rep(argu,d)
     valli<-additive(x,y,arg,h=h,M=M,eval=eval)$valvec[direc] 
     value[i]<-valli
}

pcf<-list(
value=value,
down=down,high=high,
support=support,N=N)

}

return(pcf)
}
pcf.kernesti.der<-function(x,y,h,N,kernel="gauss",support=NULL,direc=1,
method="ratio")
{
d<-length(N)
n<-length(y)

if (d>1){

if (kernel=="gauss"){
   ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
   dker<-function(xx){ 
         return( -(2*pi)^(-d/2)*xx[,direc]*exp(-rowSums(xx^2)/2) ) }
}
if (kernel=="gauss") radi<-2*h else radi<-h

recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

support<-matrix(0,2*d,1)
for (i in 1:d){
    support[2*i-1]<-min(x[,i])
    support[2*i]<-max(x[,i])
}
lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2
     argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
     neigh<-(rowSums((argu-x)^2) <= radi^2)
     if (sum(neigh)>=2){     # if there are obs in the neigborhood

       xred<-x[neigh,]
       yred<-y[neigh]
       argu<-matrix(arg,dim(xred)[1],d,byrow=TRUE)

       we<-ker((argu-xred)/h)/h^d
  
       w<-we/sum(we)
       u<-dker((argu-xred)/h)/h^(d+1)
       q<-1/sum(we)*(u-w*sum(u))  
  
       valli<-q%*%yred 
     }
     else valli<-mean(y)

     value[i]<-valli
     index[i,]<-inde
}

down<-index-1
high<-index

pcf<-list(
value=value,index=index,
down=down,high=high,  
support=support,N=N)

}
else{  # d==1  #########################################

y<-matrix(y,1,length(y))
x<-matrix(x,length(x),1)

if (kernel=="gauss"){
     ker<-function(xx){ return( exp(-xx^2/2) ) }
     dker<-function(xx){ return( -xx*exp(-xx^2/2) ) }
     dker2<-function(t){ return( -t*(2*pi)^(-1/2)*exp(-t^2/2) ) }
}

index<-seq(1:N)
len<-length(index)
down<-matrix(0,len,1)
high<-matrix(0,len,1)
down[,1]<-index-1
high[,1]<-index

value<-matrix(0,N,1)
if (is.null(support)){
   support<-matrix(0,2,1)
   support[1]<-min(x)
   support[2]<-max(x)
}
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

for (i in 1:N){
     inde<-i
     argu<-lowsuppo+step*inde-step/2
     w<-ker((argu-x)/h)/h^1
     we<-w/sum(w)

     u<-dker((argu-x)/h)/h^(1+1)
     q<-1/sum(w)*(u-we*sum(u))  

     value[i]<-y%*%q

     if (method!="ratio"){
        xs<-sort(x)
        ys<-sort(y)
        dife<-matrix(xs[2:n]-xs[1:(n-1)],n-1,1)
        ydife<-matrix(ys[2:n],1,n-1)
        q<-dife*dker2((argu-xs[2:n])/h)/h^(1+1)
        #q<-dker2((argu-x)/h)/h^(1+1)/n
        #value[i]<-y%*%q
        value[i]<-ydife%*%q
     }
  

}

pcf<-list(
value=value,
down=down,high=high,
support=support,N=N)

}

return(pcf)
}
pcf.kernesti.marg<-function(x,y,h,N,kernel="gauss",coordi=1)
{
#center=rep(0,dim(x)[2]),direc=c(1,rep(0,dim(x)[2]-1)),radius=1)

n<-dim(x)[1]
d<-dim(x)[2]

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ return( (rowSums(xx^2) <= 1) ) } 

index<-seq(1:N)
len<-length(index)
down<-matrix(0,len,1)
high<-matrix(0,len,1)
down[,1]<-index-1
high[,1]<-index

value<-matrix(0,N,1)
support<-matrix(0,2,1)
support[1]<-min(x[,coordi])
support[2]<-max(x[,coordi])
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

for (i in 1:N){
    arg1d<-lowsuppo+step*i  #-step/2
    q<-matrix(0,1,n)
    for (ii in 1:n){
        weet<-matrix(0,n,1)
        for (j in 1:n){
            arg<-x[j,]
            arg[coordi]<-arg1d
            arg<-matrix(arg,d,1)
            w<-kernesti.weights(arg,x,h=h,kernel=kernel)
            weet[j]<-w[ii]
        }
        q[ii]<-mean(weet)
    }    
    value[i]<-q%*%y 
}

pcf<-list(
value=value,
down=down,high=high,
support=support,N=N)

return(pcf)
}




pcf.kernesti<-function(x,y,h,N,kernel="gauss",support=NULL)
{
d<-length(N)

if (d>1){

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ return( (rowSums(xx^2) <= 1) ) } 

if (kernel=="gauss") radi<-2*h else radi<-h

recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

support<-matrix(0,2*d,1)
for (i in 1:d){
    support[2*i-1]<-min(x[,i])
    support[2*i]<-max(x[,i])
}
lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2
     argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
     neigh<-(rowSums((argu-x)^2) <= radi^2)
     if (sum(neigh)>=2){     # if there are obs in the neigborhood

       xred<-x[neigh,]
       yred<-y[neigh]
       argu<-matrix(arg,dim(xred)[1],d,byrow=TRUE)

       w<-ker((xred-argu)/h)/h^d
       w<-w/sum(w)
       valli<-w%*%yred 
     }
     else valli<-mean(y)

     value[i]<-valli
     index[i,]<-inde
}

down<-index-1
high<-index

pcf<-list(
value=value,index=index,
down=down,high=high,  
support=support,N=N)

}
else{  # d==1  #########################################

y<-matrix(y,1,length(y))
x<-matrix(x,length(x),1)

if (kernel=="gauss") ker<-function(xx){ return( exp(-xx^2/2) ) }
if (kernel=="uniform") ker<-function(xx){ return( (abs(xx) <= 1) ) }

index<-seq(1:N)
len<-length(index)
down<-matrix(0,len,1)
high<-matrix(0,len,1)
down[,1]<-index-1
high[,1]<-index

value<-matrix(0,N,1)
if (is.null(support)){
   support<-matrix(0,2,1)
   support[1]<-min(x)
   support[2]<-max(x)
}
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

for (i in 1:N){
     inde<-i
     argu<-lowsuppo+step*inde-step/2
     w<-ker((x-argu)/h)/h^d
     w<-w/sum(w)
     value[i]<-y%*%w
}

pcf<-list(
value=value,
down=down,high=high,
support=support,N=N)

}

return(pcf)
}
pcf.kernesti.slice<-function(x,y,h,N,kernel="gauss",
coordi=1,p=0.5,
center=NULL,direc=NULL,radius=NULL)
{
# the slice goes through vector "vecci"
#center=rep(0,dim(x)[2]),direc=c(1,rep(0,dim(x)[2]-1)),radius=1)

n<-dim(x)[1]
d<-dim(x)[2]

if (is.null(center)){
    sl<-slice.vec(x,coordi=coordi,p=p)
    center<-sl$center
    direc<-sl$direc
    radius<-sl$radius
}

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ return( (rowSums(xx^2) <= 1) ) } 

index<-seq(1:N)
len<-length(index)
down<-matrix(0,len,1)
high<-matrix(0,len,1)
down[,1]<-index-1
high[,1]<-index

value<-matrix(0,N,1)
support<-matrix(0,2,1)
support[1]<--radius
support[2]<-radius
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

for (i in 1:N){
     arg1d<-lowsuppo+step*i  #-step/2
     argDd<-center+arg1d*direc
     arg<-matrix(argDd,n,d,byrow=TRUE)
     w<-ker((x-arg)/h)/h^d
     w<-w/sum(w)
     value[i]<-w%*%y 
}

pcf<-list(
value=value,
down=down,high=high,
support=support,N=N)


return(pcf)
}




pcf.kern.quan<-function(x,y,h,N,p=0.5,kernel="gauss",support=NULL)
{
d<-length(N)

if (d>1){

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ return( (rowSums(xx^2) <= 1) ) } 

if (kernel=="gauss") radi<-2*h else radi<-h

recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

support<-matrix(0,2*d,1)
for (i in 1:d){
    support[2*i-1]<-min(x[,i])
    support[2*i]<-max(x[,i])
}
lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2
     argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
     neigh<-(rowSums((argu-x)^2) <= radi^2)
     if (sum(neigh)>=2){     # if there are obs in the neigborhood

       xred<-x[neigh,]
       yred<-y[neigh]
       #argu<-matrix(arg,dim(xred)[1],d,byrow=TRUE)
       
       valli<-kernesti.quantile(arg,xred,yred,h=h,p=p)
     }
     else valli<-mean(y)

     value[i]<-valli
     index[i,]<-inde
}

down<-index-1
high<-index

pcf<-list(
value=value,index=index,
down=down,high=high,  
support=support,N=N)

}
else{  # d==1  #########################################

y<-matrix(y,1,length(y))
x<-matrix(x,length(x),1)

if (kernel=="gauss") ker<-function(xx){ return( exp(-xx^2/2) ) }
if (kernel=="uniform") ker<-function(xx){ return( (abs(xx) <= 1) ) }

index<-seq(1:N)
len<-length(index)
down<-matrix(0,len,1)
high<-matrix(0,len,1)
down[,1]<-index-1
high[,1]<-index

value<-matrix(0,N,1)
if (is.null(support)){
   support<-matrix(0,2,1)
   support[1]<-min(x)
   support[2]<-max(x)
}
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

for (i in 1:N){
     inde<-i
     argu<-lowsuppo+step*inde-step/2
     value[i]<-kernesti.quantile(argu,x,y,h=h,p=p)
}

pcf<-list(
value=value,
down=down,high=high,
support=support,N=N)

}

return(pcf)
}
pcf.loclin<-function(x,y,h,N,type=0,kernel="gauss",support=NULL,
alt=FALSE,alt2=FALSE)
{
# type=0 ,jos regfunc, type=1, jos 1. muuttujan osit.deriv,
d<-length(N)
n<-length(y)

if (d>1){

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( (2*pi)^(-d/2)*exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ return( (rowSums(xx^2) <= 1) ) } 

if (kernel=="gauss") radi<-2*h else radi<-h

recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

if (is.null(support)){
  support<-matrix(0,2*d,1)
  for (i in 1:d){
     support[2*i-1]<-min(x[,i])
     support[2*i]<-max(x[,i])
  }
}
lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2

     argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
     w<-ker((x-argu)/h)/h^d
     painot<-w/sum(w)

     X<-cbind(matrix(1,n,1),x-argu)
     W<-diag(painot)
     A<-t(X)%*%W%*%X     
     invA<-solve(A,diag(rep(1,d+1))) 
     B<-t(X)%*%W%*%y
     valli<-invA%*%B

     value[i]<-valli[1+type]
     index[i,]<-inde
}

down<-index-1
high<-index

pcf<-list(
value=value,index=index,
down=down,high=high,  
support=support,N=N)

}
else{  # d==1  #########################################

x<-matrix(x,length(x),1)

if (kernel=="gauss") ker<-function(xx){ return( (2*pi)^(-1/2)*exp(-xx^2/2) ) }
if (kernel=="uniform") ker<-function(xx){ return( (abs(xx) <= 1) ) }
if (kernel=="bart") ker<-function(xx){ (1-xx^2)*(xx^2<=1) }

index<-seq(1:N)
len<-length(index)
down<-matrix(0,len,1)
high<-matrix(0,len,1)
down[,1]<-index-1
high[,1]<-index

value<-matrix(0,N,1)
if (is.null(support)){
   support<-matrix(0,2,1)
   support[1]<-min(x)
   support[2]<-max(x)
}
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

for (i in 1:N){
     inde<-i
     argu<-lowsuppo+step*inde-step/2
     w<-ker((x-argu)/h)/h
     painot<-w/sum(w)
     #if (!is.null(painot)) value[i]<-t(painot)%*%w else value[i]<-mean(w)

     X<-cbind(matrix(1,n,1),x-argu)
     W<-diag(c(painot))   # huom matriisi muutetaan jonoksi
     A<-t(X)%*%W%*%X     
     invA<-solve(A,diag(rep(1,d+1))) 
     B<-t(X)%*%W%*%y
     valli<-invA%*%B

     value[i]<-valli[1+type]
     index[i]<-inde     

     if (alt2==TRUE){
       s2<-c(t(painot)%*%x^2)
       s1<-c(t(painot)%*%x)
        if (type==0){
           q0<-painot*(s2-s1*x)/(s2-s1^2)
           q1<-(painot*x-s1)/(s2-s1^2)
           value[i]<-t(q0)%*%y+t(q1)%*%y*argu
        }
        else{
           q<-(painot*x-s1)/(s2-s1^2)
           value[i]<-t(q)%*%y
        }
         index[i]<-inde   
     }

    if (alt==TRUE){
       s2<-c(t(painot)%*%(x-argu)^2)
       s1<-c(t(painot)%*%(x-argu))
        if (type==0){
           q<-painot*(s2-s1*(x-argu))/(s2-s1^2)
           value[i]<-t(q)%*%y
        }
        else{
           mx<-c(t(painot)%*%x)
           q<-(painot*(x-mx))/(s2-s1^2)
           value[i]<-t(q)%*%y
        }
         index[i]<-inde   
     }

}

pcf<-list(
value=value,
down=down,high=high,
support=support,N=N)

}

return(pcf)
}


pcf.make<-function(func,N,support=NULL)
{
d<-length(N)

if (d>1){

recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

if (func=="phi"){

  phi<-function(x){ return( (2*pi)^(-d/2)*exp(-sum(x^2)/2) ) }
 
  if (is.null(support)){
       support<-matrix(0,2*d,1)
       for (i in 1:d){
           support[2*i-1]<--3
           support[2*i]<-3
       }
  }
  lowsuppo<-matrix(0,d,1)
  for (i in 1:d) lowsuppo[i]<-support[2*i-1]
  step<-matrix(0,d,1)
  for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]


  for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2   # arg<-matrix(arg,d,1)
 
     valli<-phi(arg)    

     value[i]<-valli
     index[i,]<-inde
  }

}


down<-index-1
high<-index

pcf<-list(
value=value,index=index,
down=down,high=high,  
support=support,N=N)

}

return(pcf)
}

pcf.single.index<-function(x,y,h,N,kernel="gauss",support=NULL,
method="poid",argd=colMeans(x),type="si")
{
d<-length(N)

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ return( (rowSums(xx^2) <= 1) ) } 

recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

if (is.null(support)){
support<-matrix(0,2*d,1)
for (i in 1:d){
    support[2*i-1]<-min(x[,i])
    support[2*i]<-max(x[,i])
}
}

lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

if (type=="si"){

theta<-single.index(x,y,h=h,method=method,argd=argd,kernel=kernel)
xcur<-x%*%theta
for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2
     arg<-matrix(arg,d,1)
     acur<-sum(arg*theta)
     valli<-kernesti.regr(acur,xcur,y,h=h,kernel=kernel)
    
     value[i]<-valli
     index[i,]<-inde
}

}
else{

for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2
     arg<-matrix(arg,d,1)
     theta<-single.index(x,y,h=h,method="poid",argd=arg,kernel=kernel)
     acur<-sum(arg*theta) 
     xcur<-x%*%theta
     valli<-kernesti.regr(acur,xcur,y,h=h,kernel=kernel)
     value[i]<-valli
     index[i,]<-inde
}

}

down<-index-1
high<-index

pcf<-list(
value=value,index=index,
down=down,high=high,  
support=support,N=N)

return(pcf)
}


pp.regression<-function(x,y=NULL,arg=NULL,residu=NULL,teet=NULL,
h=1,kernel="gauss",M=2,method="poid",argd=NULL,vect=FALSE,seed=1)
{
set.seed(seed)
obsit<-ceiling(runif(M)*n)

d<-dim(x)[2]
n<-dim(x)[1]

if (is.null(residu)){
  residu<-matrix(0,n,M)
  teet<-matrix(0,M,d)
  eval<-matrix(0,n,1)
  for (m in 1:M){
    residu[,m]<-y-eval
    ycur<-y-eval

    obs<-obsit[m]
    argd<-x[obs,]

    theta<-single.index(x,ycur,h=h,method=method,argd=argd)
    teet[m,]<-theta
    xcur<-x%*%theta
    if (!vect){
       estimat<-matrix(0,n,1)
       for (nn in 1:n){
            #arg<-x[nn,]
            #arg<-matrix(arg,d,1)
            #acur<-sum(arg*theta)
            acur<-xcur[nn]
            est<-kernesti.regr(acur,xcur,ycur,h=h,kernel=kernel)
            estimat[nn]<-est
       }
    }
    else{
        estimat<-kernesti.regr(xcur,xcur,ycur,h=h,kernel=kernel,vect=vect)       
    }
    eval<-eval+estimat
  }
}
else eval<-NULL

if (is.null(arg)){ 
  value<-NULL
}
else{
  value<-0
  for (m in 1:M){
     ycur<-matrix(residu[,m],n,1)
     tcur<-matrix(teet[m,],d,1)
     xcur<-x%*%tcur
     arg<-matrix(arg,d,1)
     carg<-sum(arg*tcur)
     w<-kernesti.weights(carg,xcur,h=h,kernel=kernel)
     curvalue<-t(w)%*%ycur
     value<-value+curvalue
  }
}

return(list(eval=eval,residu=residu,teet=teet,value=value))
}



quantil.emp<-function(y,p)
{
n<-length(y)
or<-order(y)
m<-ceiling(p*n)
quant<-y[or[m]]

return(quant)
}


single.index.aved<-function(x,y,h=1,kernel="gauss",argd=NULL,take=length(y),seed=1)
{
d<-dim(x)[2]
n<-length(y)

if (take==n) ota<-seq(1,n) else{
   set.seed(seed)
   ota<-ceiling(n*runif(take))
}

if (is.null(argd)){
  grad<-matrix(0,n,d)
  for (ii in 1:d){
     for (jj in ota){
        arg<-x[jj,]
        grad[jj,ii]<-kernesti.der(arg,x,y,h=h,direc=ii,kernel=kernel,vect=FALSE)
     }
  }
  theta<-colMeans(grad)
}
else{
  grad<-matrix(0,d,1)
  for (ii in 1:d){
        grad[ii]<-kernesti.der(argd,x,y,h=h,direc=ii,kernel=kernel,vect=FALSE)
  }
  theta<-grad
}

theta<-theta/sqrt(sum(theta^2))

return(theta)
}

single.index.gene<-function(x,y,h,kernel="gauss")
{
n<-dim(x)[1]
d<-dim(x)[2]
y<-matrix(y,n,1)

if (kernel=="bart") 
   ker<-function(t){ return( (1-t) ) }
if (kernel=="gauss") 
   ker<-function(t){ return( exp(-t/2) ) }
if (kernel=="uniform") 
   ker<-function(t){ return( (t <= 1) ) } 

fn<-function(b) { 
       z<-x%*%b           # z is a column vector of new explanatory variables

       A<-z%*%t(z)
       B<-matrix(diag(A),n,n)
       C<-B-2*A+t(B)           # C is the symmetric n*n matrix of mutual 
                               # squared distances among the elements of z
       D<-ker(C/h)/h           # D is the n*n-matrix of weights; row i
                               # of D is a vector of weights associated to 
                               # argument z_i

       W<-D/colSums(D)  # rows of W sum to one, i:th row is the normalized
                        # vector of weights associated to argument z_i

       error<-sum((W%*%y-y)^2)
       return(error) 
}

par<-rep(1,d)/d           # initial value
par.lower<-rep(-1,d)
par.upper<-rep(1,d)

op.method<-"L-BFGS-B"
op<-optim(par=par,fn=fn,method=op.method,lower=par.lower,upper=par.upper)
theta<-op$par

#nlin<-list( function(b){ return( sum(b^2) ) } )
#nlin.lower<-1
#nlin.upper<-1
#control<-donlp2.control(silent=TRUE)
#curp<-donlp2(par=par,fn=fn, # par.lower=par.lower, par.upper=par.upper,
#             nlin=nlin, nlin.upper=nlin.upper, nlin.lower=nlin.lower, 
#             control=control)
#theta<-curp$par

theta<-theta/sqrt(sum(theta^2))
return(theta)
}
single.index.itera<-function(x,y,h=1,M=2,kernel="gauss",vect=FALSE)
{
d<-dim(x)[2]
n<-length(y)
theta0<-matrix(1,d,1)
theta0<-theta0/sqrt(sum(theta0^2))

w<-matrix(0,n,1)
haty<-matrix(0,n,1)
for (m in 1:M){

   xcur<-x%*%theta0
   for (i in 1:n) w[i]<-kernesti.der(xcur[i],xcur,y,h=h,vect=vect)
   weights<-w^2
   for (i in 1:n) haty[i]<-kernesti.regr(xcur[i],xcur,y,h=h,vect=vect)
   z<-xcur+(y-haty)/w

   W<-diag(c(weights),nrow=n,ncol=n)
   A<-t(x)%*%W%*%x     
   invA<-solve(A,diag(rep(1,d))) 
   B<-t(x)%*%W%*%z
   theta0<-invA%*%B
   theta0<-theta0/sqrt(sum(theta0^2))
}

return(theta0)
}

single.index<-function(x,y,arg=NULL,h=1,kernel="gauss",
M=2,method="iter",vect=FALSE,argd=arg,take=length(y),seed=1)
{
d<-dim(x)[2]

if (method=="nume") 
theta<-single.index.gene(x,y,h=h,kernel=kernel) 

if (method=="iter") 
theta<-single.index.itera(x,y,h=h,kernel=kernel,M=M,vect=vect) 

if (method=="aved") 
theta<-single.index.aved(x,y,h=h,kernel=kernel,take=take,seed=seed) 

if (method=="poid"){
   if (is.null(argd)) argd<-colMeans(x)
   theta<-single.index.aved(x,y,h=h,kernel=kernel,argd=argd) 
}

if (!is.null(arg)){
  xcur<-x%*%theta
  arg<-matrix(arg,d,1)
  acur<-sum(arg*theta)
  est<-kernesti.regr(acur,xcur,y,h=h,kernel=kernel,vect=vect)
  return(est)
}
else return(theta)
  
}




slice.vec<-function(x,coordi=1,p=0.5)
{
d<-dim(x)[2]
center<-matrix(0,d,1)
for (i in 1:d) center[i]<-quantile(x[,i],p=p)  #center<-colMeans(x)
direc<-rep(0,d)
direc[coordi]<-1
radius<-(max(x[,coordi])-min(x[,coordi]))/2

return(list(center=center,direc=direc,radius=radius))
}



