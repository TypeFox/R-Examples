dpglasso <-
function(Sigma, X=NULL,invX=NULL,rho,outer.Maxiter=100,obj.seq=FALSE,outer.tol=10^-5) {

if (is.null(Sigma)) { print("Sigma is required as input"); return(); } else { p=nrow(Sigma) }

check1<-(nrow(Sigma)==ncol(Sigma))&&(Sigma==t(Sigma));
if (check1==0) { print("Sigma sould be square and symmetric"); return(); } 

if(is.null(X))  X= diag( 1/(rep(rho,p) + diag(Sigma)) ) 

if(is.null(invX)) { U.mat = diag(rep(rho,p)); invX<- Sigma + U.mat } else {U.mat <- invX - Sigma; }

p1<-nrow(X); p11<-ncol(X); p2<-nrow(invX); p22<-ncol(invX);

if (p1!=p11) { print("check dimensions of X"); return(); }
if (p2!=p22) { print("check dimensions of invX"); return(); }


check.dim<- (p1==p2)&&(p1==p); 
if (check.dim==0) { print("check dimensions of X,invX,Sigma"); return(); }

ids<-rep(1:p,outer.Maxiter)
time.counter.QP<-array(0,dim=c(length(ids),3));
rel.err<-rep(0,outer.Maxiter);
obj.vals<-rep(0,outer.Maxiter);
sparse.nos<-rep(0,outer.Maxiter);

ii<-0; kk<-0;
tol=10;

vec.diagsold<-X;

for (cd.iter in ids) {
##print(ii)
ii<-ii+1;
t<-proc.time()
obj<-box_qp_f(X[-cd.iter,-cd.iter],u=U.mat[-cd.iter,cd.iter],b=Sigma[-cd.iter,cd.iter],rho,
              Maxiter=1000,tol=10^-7); 
t<- (proc.time() - t);

time.counter.QP[ii,]<- as.numeric(c(t[1],t[2],t[3]));
theta.hat<--obj$grad_vec*0.5/(Sigma[cd.iter,cd.iter]+ rho);
theta.hat[abs(obj$u)< rho*.99999]=0; theta.hat[abs(theta.hat)< 10^-5]=0;

X[cd.iter,-cd.iter]<-theta.hat;
X[-cd.iter,cd.iter]<-X[cd.iter,-cd.iter];
X[cd.iter,cd.iter]=sum((obj$u+Sigma[cd.iter,-cd.iter])*X[cd.iter,-cd.iter]);
X[cd.iter,cd.iter]<-(1 - X[cd.iter,cd.iter])/(Sigma[cd.iter,cd.iter]+ rho)

U.mat[-cd.iter,cd.iter]<-as.numeric(obj$u,ncol=1); U.mat[cd.iter,-cd.iter]<-U.mat[-cd.iter,cd.iter];
U.mat[cd.iter,cd.iter]=rho;


if (cd.iter==p) {
  kk<-kk+1;
#  vec.diagsnew<-diag(X); 
vec.diagsnew<-X;
  tol = max(abs(vec.diagsnew-vec.diagsold))/max(abs(vec.diagsold));
  rel.err[kk]<-tol;
  vec.diagsold<-vec.diagsnew;

sparse.nos[kk]<-sum(abs(X)<=10^-9);
  
  if (obj.seq==TRUE) obj.vals[kk]<- - sum(log(abs(diag(qr(X)$qr))))  + sum(Sigma*X) + rho*sum(abs(X))

  if ((tol<outer.tol)&(kk>1))  break
                }

                    }

time.counter.QP=time.counter.QP[1:ii,];
time.counter.QP=colSums(time.counter.QP);
if (obj.seq==TRUE) {
return(list(time.counter.QP=time.counter.QP,rel.err=rel.err[1:kk],obj.vals=obj.vals[1:kk],
            X=X,invX= Sigma + U.mat,sparse.nos=sparse.nos[1:kk]))
                  } else {
return(list(time.counter.QP=time.counter.QP,rel.err=rel.err[1:kk],
            X=X,invX= Sigma + U.mat,sparse.nos=sparse.nos[1:kk]))
}

                                                                                }
