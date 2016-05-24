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
