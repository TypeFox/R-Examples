pcf.func<-function(func, N,
sig=rep(1,length(N)), support=NULL, theta=NULL, 
g=1, M=NULL, p=NULL, mul=3, t=NULL, 
marginal="normal", r=0,
mu=NULL, xi=NULL, Omega=NULL, alpha=NULL, df=NULL, 
a=0.5, b=0.5, distr=FALSE, std=1, lowest=0) # contrast="loglik")   
{
# t<-rep(1,length(N))

d<-length(N)

if (d>1){

  if (marginal=="unif") support<-c(0,sig[1],0,sig[2])

  recnum<-prod(N)
  value<-matrix(0,recnum,1)
  index<-matrix(0,recnum,d)

  # new ############################################

  if (func=="mixt"){ 

     if (is.null(support)){
       support<-matrix(0,2*d,1)
       for (i in 1:d){
           support[2*i-1]<-min(M[,i]-mul*sig[,i])
           support[2*i]<-max(M[,i]+mul*sig[,i])
       }
     }
     lowsuppo<-matrix(0,d,1)
     for (i in 1:d) lowsuppo[i]<-support[2*i-1]
     step<-matrix(0,d,1)
     for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]
     mixnum<-length(p)

     numpositive<-0
     for (i in 1:recnum){
        inde<-digit(i-1,N)+1
        point<-lowsuppo+step*inde-step/2
 
        if (!is.null(theta)){
           rotmat<-matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
           point<-rotmat%*%point
        }

        valli<-0
        for (mi in 1:mixnum){
            evapoint<-(point-M[mi,])/sig[mi,]
            valli<-valli+p[mi]*evanor(evapoint)/prod(sig[mi,])
        }
        if (valli>lowest){
           numpositive<-numpositive+1
           value[numpositive]<-valli
           index[numpositive,]<-inde
        }
     }
     value<-value[1:numpositive]
     index<-index[1:numpositive,]
     down<-index-1
     high<-index
  }


  else if (func=="student"){ 
     lowsuppo<-matrix(0,d,1)
     for (i in 1:d) lowsuppo[i]<-support[2*i-1]
     step<-matrix(0,d,1)
     for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

     numpositive<-0
     for (i in 1:recnum){
        inde<-digit(i-1,N)+1
        x<-lowsuppo+step*inde-step/2

        #valli<-eva.student(x,t,marginal,sig,r,df)

        margx<-matrix(0,d,1)
        u<-matrix(0,d,1)

        if (marginal=="unif"){
           for (j in 1:d){
             u[j]<-x[j]/sig[j]  #+1/2
             margx[j]<-1/sig[j]
           }
        }
        if ((marginal=="normal")||(marginal=="gauss")){
           for (j in 1:d){
             u[j]<-pnorm(x[j]/sig[j])
             margx[j]<-evanor(x[j]/sig[j])/sig[j]
           }
        }
        if (marginal=="student"){
          for (j in 1:d){
             u[j]<-pt(x[j]/sig[j],df=t[j])
             margx[j]<-dt(x[j]/sig[j],df=t[j])/sig[j]
          }
        }
        
        x1<-qt(u[1],df=df)
        x2<-qt(u[2],df=df)

        d<-2
        vakio<-gamma((df+d)/2)*gamma(df/2)/gamma((df+1)/2)^2
        nelio<-(x1^2+x2^2-2*r*x1*x2)/(1-r^2)
        prod<-(1+x1^2/df)^((1+df)/2)*(1+x2^2/df)^((1+df)/2)
        copuval<-vakio*(1-r^2)^(-1/2)*prod*(1+nelio/df)^(-(df+d)/2)

        valli<-copuval*margx[1]*margx[2]

        ###############################################

        if (valli>0){
           numpositive<-numpositive+1
           value[numpositive]<-valli
           index[numpositive,]<-inde
        }
     }
     value<-value[1:numpositive]
     index<-index[1:numpositive,]
     down<-index-1
     high<-index
  }

  else if (func=="gauss"){ 
     lowsuppo<-matrix(0,d,1)
     for (i in 1:d) lowsuppo[i]<-support[2*i-1]
     step<-matrix(0,d,1)
     for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

     numpositive<-0
     for (i in 1:recnum){
        inde<-digit(i-1,N)+1
        x<-lowsuppo+step*inde-step/2

        #valli<-eva.copula(x,type="gauss",marginal=marginal,sig=sig,r=r,t=t)

        margx<-matrix(0,d,1)
        u<-matrix(0,d,1)

        if (marginal=="unif"){
           for (j in 1:d){
             u[j]<-x[j]/sig[j]  #+1/2
             margx[j]<-1/sig[j]
           }
        }
        if ((marginal=="normal")||(marginal=="gauss")){
           for (j in 1:d){
             u[j]<-pnorm(x[j]/sig[j])
             margx[j]<-evanor(x[j]/sig[j])/sig[j]
           }
        }
        if (marginal=="student"){
          for (j in 1:d){
             u[j]<-pt(x[j]/sig[j],df=t[j])
             margx[j]<-dt(x[j]/sig[j],df=t[j])/sig[j]
          }
        }
        
        x1<-qnorm(u[1],sd=1)
        x2<-qnorm(u[2],sd=1)

        nelio<-(x1^2+x2^2-2*r*x1*x2)/(1-r^2)
        copuval<-(1-r^2)^(-1/2)*exp(-nelio/2)/exp(-(x1^2+x2^2)/2)

        valli<-copuval*margx[1]*margx[2]

        ########################################

        if (valli>0){
           numpositive<-numpositive+1
           value[numpositive]<-valli
           index[numpositive,]<-inde
        }
     }
     value<-value[1:numpositive]
     index<-index[1:numpositive,]
     down<-index-1
     high<-index
  }


else{

# old #########################################################

if (is.null(support)){

   if (func=="epan"){
      if (is.null(sig)) sig<-c(1,1)
      support<-matrix(0,2*d,1)
      for (i in 1:d){
          support[2*i-1]<--sig[i]
          support[2*i]<-sig[i]
      }
   }

}

if ((marginal=="unif")) support<-c(0,sig[1],0,sig[2])
# && (is.null(support))) 
#support<-c(-sig[1]/2,sig[1]/2,-sig[2]/2,sig[2]/2)


lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

numpositive<-0
for (i in 1:recnum){
    inde<-digit(i-1,N)+1
    #if ((inde[1]==0) && (inde[2]==N[2])) inde<-c(0,0)
    point<-lowsuppo+step*inde-step/2

    if (!is.null(theta)){
         rotmat<-matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
         point<-rotmat%*%point
    }

    if (func=="prod") valli<-eva.prod(point,marginal,g)
    if (func=="skewgauss") valli<-eva.skewgauss(point,mu,sig,alpha)
    #if (func=="dmsn") valli<-dmsn(point,xi,Omega,alpha)
    if (func=="gumbel") valli<-eva.copula(point,
        type="gumbel",marginal=marginal,sig=sig,r=r,t=t,g=g)
    if (func=="frank") valli<-eva.copula(point,
        type="frank",marginal=marginal,sig=sig,t=t,g=g)
    if (func=="plackett") valli<-eva.plackett(point,t,marginal,sig)
    if (func=="clayton2") valli<-eva.clayton(point,t,marginal,sig,df)
    if (func=="clayton") valli<-eva.copula(point,
        type="clayton",marginal=marginal,sig=sig,r=r,t=t,g=g)
    if (func=="cop6") valli<-eva.cop6(point,t,marginal,sig)
    if (func=="epan") valli<-epan(point)
    if (func=="normal") 
        valli<-eva.gauss(point,t=t,marginal=marginal,sig=sig,r=r)   
    if (func=="hat") valli<-eva.hat(point,a=a,b=b)

    if (valli>0){
       numpositive<-numpositive+1
       value[numpositive]<-valli
       index[numpositive,]<-inde
    }
}

value<-value[1:numpositive]
index<-index[1:numpositive,]
down<-index-1
high<-index

}


pcf<-list(
value=value,index=index,
down=down,high=high,  #step=delta,
support=support,N=N)

  #pcf<-eval.func.dD(func,N,
  #sig=sig,support=support,theta=theta,g=g,
  #M=M,p=p,mul=mul,
  #t=t,marginal=marginal,r=r, 
  #mu=mu,xi=xi,Omega=Omega,alpha=alpha,df=df,a=a,b=b)

}

else{  # (d==1){ ######################################################

  pcf<-eval.func.1D(func,N,
  support=support,g=g,std=std,distr=distr,
  M=M,sig=sig,p=p,
  a=a,b=b,d=2)

}


return(pcf)

}

