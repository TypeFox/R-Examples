
isiindata<-function(i,x){
n<-nrow(x)
ind<-rep(NA,len=n)
for (j in 1:n){
 ind[j]<-all(i==x[j,])
}
ind
}
triang=function(Q,P1,P2,P3){
    q1=Q-P1;a1=atan2(q1[2],q1[1])
    q2=Q-P2;a2=atan2(q2[2],q2[1])
    q3=Q-P3;a3=atan2(q3[2],q3[1])
    a=sort(c(a1,a2,a3))
    return(((a[2]-a[1])<=pi) && ((a[3]-a[2])<=pi) && ((a[3]-a[1])>=pi))
}

 ntriang<-function(q,puntos,plot=FALSE){
#print("entra ntriang")
#   if (dim(x)[2]!=2) stop("The dimension of the data points is not a 2-column matrix")
#   if (dim(x)[1]<3) stop("The dimension of the data points must be at least 3 rows.")
#   if (length(q)!=2) stop("El punto q debe tener dos componentes")
n=dim(puntos)[1]
#nq=q-puntos
nq<-sweep(puntos,2,q,"-")
alfa=atan2(nq[,2],nq[,1])
a=sort(alfa)
di<-ci<-numeric(n);di=numeric(n)
#ci=di=numeric(n)
ind.a<-a>0
k<-ifelse(sum(ind.a)>0,min(which(ind.a)),n)
#k=min(min(which(a>0)),n)
#if (is.infinite(k)) k=n
lista=1:(k-1)
for (i in 1:n){
    ci[i]=max(which((a-a[i])<=pi))
    di[i]=min(c(which((a-a[i])>=pi)),n+1)
    }
e1=ci[k]*(ci[k]+1-2*k)+k*(3-k)-2
d1=n*(n-k+1)+0.5*k*(k-1)-0.5*n*(n+1)
e=(lista+n-1)*ci[lista]
f=(ci[lista]+1)*ci[lista]+(di[lista]-1)*(di[lista]-2*lista-2)
sume=sum(e)+0.5*n*e1
sumf=sum(f/2)+n*d1
ttriang=n*(n-1)*(n-2)/6 #n sobre 3
return(list("triang"=sume-sumf,"ttriang"=ttriang))
}

mdepth.SD=function(x,xx=NULL,scale=FALSE){
   if (dim(x)[2]!=2) stop("The dimension of the data points is not a 2-column matrix")
   if (dim(x)[1]<3) stop("The dimension of the data points must be at least 3 rows.")
   if (is.data.frame(x)) {
      x<-as.matrix(x)
    }
    if (is.null(xx)) {
      xx<-x
      al<-TRUE
      }
   else {
    al<-FALSE
    xx<-as.matrix(xx)
    }
   nam<-colnames(xx)

       if ( is.null(rownames(x)))  rownames(x)<-1:nrow(x)
       nms<-rownames(x)
       m0<-nrow(x)
       xx<-na.omit(xx)
       x<-na.omit(x)
       nas<-na.action(x)
       nullans<-!is.null(nas) 
          
   n=nrow(x)
   nn=nrow(xx)
   rownames(xx)<-NULL#1:nn
   d=ncol(xx)
   dep=numeric(n)
 if (al){
# print("entra al l")
    ttriang=nn*(nn-1)*(nn-2)/6    
#    print(al);    print("all(x==xx)")
    fac<-(nn-1)*(nn-2)/2
    for (i in 1:n){
#        dep[i]=(ntriang(x[i,],x[-i,],TRUE)$triang+(n-1)*(n-2)/2)/ttriang
        dep[i]<-ntriang(x[i,],xx[-i,])$triang      
#        cat("i ",i);        print(dep[i])
   }
   dep<-(dep+fac)/ttriang
   }
  else{
    ttriang=(nn+1)*(nn)*(nn-1)/6
    fac<-(nn)*(nn-1)/2
    for (i in 1:n){
          dep[i]=ntriang(x[i,],xx)$triang
        }      
   dep<-(dep+fac)/ttriang 
   }
   if (scale){ 
    dep<-dep/max(dep)
    } 
 if  (nullans){
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-dep
        dep<-ans1      
        }
   names(dep)<-nms         

 return(invisible(list(dep = dep)))
}

