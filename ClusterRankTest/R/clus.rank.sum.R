clus.rank.sum <-
function(Cluster,X,grp=NULL,Y=NULL,test=c("DS","DD","SDS")) {

if(as.integer(is.null(grp))==1 && as.integer(is.null(Y))==1) {stop("exactly one of 'grp' and 'Y' must be null")}


if(test=="DS" && as.integer(is.null(grp))==0){
#####calculate quantity 2 (using the pooled estimate of F)
      n<-length(X)
      F.hat<-numeric(n)
      for (i in 1:n){
           F.hat[i]<-(sum(X<=X[i])+sum(X<X[i]))/(2*n)
                    }
#####calculate quantity 1 (using ECD-F for each cluster)
#### M is No. of clusters, n is No. of observations
     M<-length(unique(Cluster)) 
     n.i<-table(Cluster)
     F.prop<-numeric(n)
     for(ii in 1:n){
               F.j<-numeric(M)
               for (i in 1:M){
                    F.j[i]<-(sum(X[Cluster==i]<X[ii])+0.5*sum(X[Cluster==i]==X[ii]))/(n.i[i])
                             } 
               F.prop[ii]<-sum(F.j[-Cluster[ii]])
                   }   

###########calculate S=E(W*|X,g)
a<-numeric(M)
b<-1+F.prop
for (i in 1:M){
      a[i]<-sum((grp[Cluster==i]*b[Cluster==i])/(n.i[i]))
                  }    
c<-1/(M+1)
S<-c*sum(a)
########note: for m groups maybe can use grp[Cluster==i&grp=m]

#########Calculate E(S)=E(W*)
n.i1<-table(Cluster[grp==1])
d<-n.i1/n.i
E.S<-(1/2)*sum(d)

#######Calculate estimate of variance of S
W.hat<-numeric(M)        #####first calculate W.hat for each cluster
a<-n.i1/n.i
for (i in 1:M){
    b<-1/(n.i[i]*(M+1))
    c<-(grp[Cluster==i])*(M-1)
    d<-sum(a[-i])
    W.hat[i]<-b*sum((c-d)*F.hat[Cluster==i])
              }
a<-n.i1/n.i
E.W<-(M/(2*(M+1)))*(a-sum(a)/M)    ##second, calculate E(W)

var.s<-sum((W.hat-E.W)^2) #calculate var(s)
statistic <-(S-E.S)/sqrt(var.s)   #calculate the test statistic
p.value<-2*pnorm(abs(statistic),lower.tail=F)
}


if(test=="DD" && as.integer(is.null(grp))==0){

################# new rank sum test ########################

rn<-function(dv){
ik=dv[1]
x=dv[2]
ds1=data[data[,3]==1,]
vs1=(kh==2)*(ds1[,2]<x)+(kh==1)*(ds1[,2]<=x)
sl1=aggregate(vs1,list(ds1[,1]),mean)[,2]

ds2=data[data[,3]==0,]
vs2=(kh==2)*(ds2[,2]<x)+(kh==1)*(ds2[,2]<=x)
sl2=aggregate(vs2,list(ds2[,1]),mean)[,2]

fg=(sl1+sl2)/2
fg[ik]=0
return(fg)
}

rst<-function(il){
	#only for variance
	ly=sum(mat[-which(dw[,1]==il),-il])
#ly=apply(mat[-which(dw[,1]==il),-il],1,sum)
	return(ly)
}




data= cbind(Cluster,X, grp)
m=length(unique(data[,1]))
dw=data[(data[,3]==1),]
ns=(dw[,1])
nv=as.vector(table(ns)[match(ns,names(table(ns)))])

kh=1
mat=t(apply(cbind(dw[,1:2]),1,rn))/nv
 vf1=apply(cbind(seq(1,m)),1,rst) # variance part check
sFs1=sum(mat) #-estimate part

kh=2
mat=t(apply(cbind(dw[,1:2]),1,rn))/nv
vf2=apply(cbind(seq(1,m)),1,rst) #check
sFs2=sum(mat)

v1=((sFs1+sFs2)/4)+(m/2) #estimate --matches original
vd= ((vf1+vf2)/4)+(m-1)/2 # 
   #v1 is fine
   #check vf1,vf2 and vd
   h=1
statistic <- v1
E.T<- 0.25*m*(m+1)
test=(m/m^h)*v1-((m-1)/(m-1)^h)*vd
v.test=var(test)
v_hat=(((m^h)^2)/(m-1))*v.test
v.hat=ifelse(v_hat==0,0.00000001,v_hat)
Z<- (statistic-E.T)/sqrt(v.hat)
p.value<- 2*pnorm(abs(Z), lower.tail=FALSE)
}


if(test=="SDS" && as.integer(is.null(Y))==0){
Xij <- X-Y
ni <- as.vector(table(Cluster))
g <- length(ni)
n <- sum(ni)

 


cni <- cumsum(ni)
cni <- c(0,cni)

Fi <- function(x,i) { Xi <- Xij[(cni[i]+1):(cni[i+1])];
                      (sum(abs(Xi)<=x)+sum(abs(Xi)<x))/(2*ni[i])}
Ftot <- function(x) { st <- 0;
                      for (i in 1:g) st <- st + Fi(x,i);
                      return(st)}
Fcom <- function(x) { st <- 0;
                      for (i in 1:g) st <- st + Fi(x,i)*ni[i];
                      return(st/n)}

# SIGNED RANK TEST STATISTIC
TS <- VTS <- 0
for (i in 1:g) {

Xi <- Xij[(cni[i]+1):(cni[i+1])]
first <- (sum(Xi>0)-sum(Xi<0))/length(Xi)
second <- 0
third <- 0
for (x in Xi) { second <- second + sign(x)*(Ftot(abs(x))-Fi(abs(x),i));
                third <- third + sign(x)*Fcom(abs(x))}

TS <- TS + first+second/length(Xi)
VTS <- VTS + (first+ (g-1)*third/length(Xi))^2
                }

statistic=TS/sqrt(VTS) ##
p.value= 2*(1-pnorm(abs(statistic))) ##
}

structure(list(p.value=p.value,TestStat=statistic),class="Cluster.Test")
}
