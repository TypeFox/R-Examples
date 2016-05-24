kmc.el<-function(delta,omega,S){
  n<-length(S)
    sum(llog(omega[delta==1],1/n^2/100))+sum(llog(S[delta==0],1/n^2/100))
  omega[delta==1] -> val1
  S[delta==0]     -> val2
  eps=1e-7
    sum(llog(val1,.001/n^2)[val1>eps])+sum(llog(val2,.001/n^2)[val2>eps])
}

kmc.clean <- function(kmc.time,delta){
  #TASK:
  #1 sort T
  #2 the first is uncen!
  n=length(kmc.time);
  tmp <- sort(kmc.time,index.return=TRUE);
  kmc.time=kmc.time[tmp$ix];
  delta=delta[tmp$ix];
  FirstUnCenLocation<-which(delta==1)[1];
  if (FirstUnCenLocation==n) {stop('Only one uncensored point.');}
  if (FirstUnCenLocation!=1){
    delta=delta[FirstUnCenLocation:n];
    kmc.time=kmc.time[FirstUnCenLocation:n];
  }
  delta[length(kmc.time)]=1;
  #U=kmc_find0loc(delta);
#  cat("len: ",U);
#if (U==0) stop("Not enough event points!");
#  return (list(kmc.time=kmc.time[1:U],delta=delta[1:U]));
 return (list(kmc.time=kmc.time,delta=delta));
}




omega.lambda<-cmpfun(function(kmc.time,delta,lambda,g,gt.mat){
  #iter
  p=length(g);# the number of constraint
  n=length(kmc.time);
  uncen.loc<- which(delta==1);
  cen.loc<- which(delta==0);
  delta[n]=1
  #########################################

  u.omega<-numeric(n);
  u.omega[1]<-1/(n-sum(lambda*gt.mat[,1]));
  for (k in 2:n){
    if (delta[k]==1){
      S <- 1-cumsum(u.omega);#need to update every kmc.time(add in one entry each kmc.time)
      SCenLoc<-cen.loc[cen.loc%in%(1:(k-1))];
      S.cen=0;
      if (length(SCenLoc)!=0){S.cen=sum(1/S[SCenLoc])}
      u.omega[k]=1/(n-sum(lambda*gt.mat[,k])-S.cen)
      #cat(':::',sum(omega),'\n')
    }
  }
  return(list(S=S,omega=u.omega,gt=gt.mat));
}
)

kmc.data <- cmpfun(function(kmc.time,delta,lambda,g, gt.mat,using.C=F){
  #my omega contains 0, it is length of n
  #need S omega
  if(using.C) {
      #tmp<-lambdaoo(kmc.time,delta,lambda,gt.mat);
      return( kmcdata_rcpp(kmctime=kmc.time,delta=delta,lambda=lambda,gtmat=gt.mat));
  }else{
      tmp<-omega.lambda(kmc.time,delta,lambda,g, gt.mat)
  }
  b= delta * tmp$omega;
#check.constriant=apply(t(b*t(tmp$gt)),1,sum);
  check.constriant=rowSums(t(b*t(tmp$gt)));
  gama=1/(tmp$S);
  return (list(omega=tmp$omega,gamma=gama,S=tmp$S,chk=check.constriant));
})

omega.lambda12<-cmpfun(function(kmc.time,delta,lambda,g, gt.mat){
  #iter
  cumsumx<-function(x) apply(x,1,cumsum)
  p=length(g);# the number of constraint
  n=length(kmc.time);
  uncen.loc<- which(delta==1);
  cen.loc<- which(delta==0);
  delta[n]=1
  #########################################
  #gt.mat<-matrix(0,p,n);
  #for(i in 1:p){for (j in 1:n) { gt.mat[i,j]=g[[i]](kmc.time[j])}}
  #for(i in 1:p) gt.mat[i,]=g[[i]](kmc.time);
  u.omega<-numeric(n);
  udev.omega<-matrix(0,p,n);
  u.omega[1]<-1/(n-sum(lambda*gt.mat[,1]));
  udev.omega[,1]<-u.omega[1]^2*gt.mat[,1];
  for (k in 2:n){
    if (delta[k]==1){
      S <- 1-cumsum(u.omega);#need to update every kmc.time(add in one entry each kmc.time)
      SCenLoc<-cen.loc[cen.loc%in%(1:(k-1))];
      S.cen=0;
      S.cen2=0;
      if (length(SCenLoc)!=0){
        S.cen=sum(1/S[SCenLoc])
        S.cen2=sum( 1/(S[SCenLoc]^2)*cumsumx(matrix(udev.omega[,1:(k-1)],ncol=k-1))[SCenLoc] );
      }
      u.omega[k]=1/(n-sum(lambda*gt.mat[,k])-S.cen)
      udev.omega[k]=u.omega[k]^2*(gt.mat[,k]+S.cen2)
      #cat(':::',sum(omega),'\n')
    }
  }
  return(list(S=S,omega=u.omega,gt=gt.mat,omega.dev=sum(delta*gt.mat*udev.omega)));
}
)

kmc.data12 <- function(kmc.time,delta,lambda,g, gt.mat){
  #my omega contains 0, it is length of n
  #need S omega
  tmp<-omega.lambda12(kmc.time,delta,lambda,g, gt.mat)
  b= delta * tmp$omega;
  #check.constriant=apply(t(b*t(tmp$gt)),1,sum);
  check.constriant=rowSums(t(b*t(tmp$gt)));
  gama=1/(tmp$S);
  return (list(omega=tmp$omega,gamma=tmp$gamma,S=tmp$S,chk=check.constriant,domega=tmp$omega.dev));
}



kmc.solve<-function(x,d,g,em.boost=T,using.num=T,using.Fortran=T,using.C=F,tmp.tag=T,rtol=1E-9,control=list(nr.it=20,nr.c=1,em.it=3),...){
    #n=length(x)
    ###### checking PHASE 1         ######
    if (length(unique(d))!=1 ){if (!setequal(unique(d),c(0,1))) stop("Status must be 0/1");}else{if (d[1]!=1) stop("Status must be 0/1");}# check d = 0/1 or all 1's
    if (sum(d)<length(g)) stop("Number of observation MUST be greater than numbers of constraints");
    if ('nr.it'%in%names(control)){ nr.it=control$nr.it; if (nr.it<10) nr.it=10}else{nr.it=20} # NR iteration
    if ('nr.c'%in%names(control)){nr.c=control$nr.c;if (nr.c>1) {warning("In N-R iteration, C should be between 0 and 1");nr.c=1}  }else{nr.c=1} #NR scaler
    if ('em.it'%in%names(control)){ em.it=control$em.it; if (em.it>10) em.it=10}else{em.it=3} # EM iteration
    
    ###### end of checking PHASE 1  ######
    kmc.clean(kmc.time=x,delta=d)->re
    kmc.time=re$kmc.time;
    delta=re$delta; ## use global now, mod latter
    p=length(g)
    if (tmp.tag)delta[1:p]=1## if two constraints
    n=length(delta)
    gt.mat<-matrix(0,p,n);
    for(i in 1:p) gt.mat[i,]=g[[i]](kmc.time);
    
    ######TODO: checking PHASE 2        ######
    ##    Is the feasible region = NULL?
    
    kmc.comb<-function(x){
      #OLD kmc.data(kmc.time,delta,lambda=x,g, gt.mat= gt.mat,using.C=using.C)-> re;
            kmc.data(kmc.time,delta,lambda=x,g, gt.mat= gt.mat,using.C=using.C)-> re;
        re$chk
    }
    
    kmc.comb12<-Vectorize(function(x){
        kmc.data12(kmc.time,delta,lambda=x,g, gt.mat= gt.mat)-> re;
      list(x=re$chk,dev=re$domega);
    })
    
	kmc.comb123<-function(x){
        kmc_routine4(lambda=x,delta=delta,gtmat=gt.mat)->re;
        return(re);
    }
	
    multiroot.nr<-function(f_,xinit,it=nr.it,C=nr.c,trace=FALSE,tol=1E-9){
      if (C*tol>1) C=ceiling(1/tol/10);
      re=xinit;
      if (trace) cat('\nx\tf(x)\tdf(x)')
      for (i in 1:it){
        f_(re)->tmp
        if ( abs(tmp[[1]])< tol ) break;
        re=re-tmp[[1]]/tmp[[2]]*C
        if (trace){
          cat('\n',re,'\t',tmp[[1]],'\t',tmp[[2]]);
        }
      }
      if (i==it){ cat('\nMay not converge.\n')}else cat('\nConverged!\n')
      re
    }
    
  if (em.boost){
    if (length(g)==1){
        u.lambda1<-function(re=el.cen.EM.kmc(x=kmc.time,d=delta,fun=g[[1]],mu=0,maxit=em.it,debug.kmc=F)){
          (n-1/re$prob[1])/g[[1]](re$times[1])
        }
        init.lam=u.lambda1()
      }else{
        if(length(g)==2){
          u.lambda2<-function(re=el.cen.EM2.kmc(x=kmc.time,d=delta,fun=function(x){cbind(g[[1]](x),g[[2]](x))},mu=c(0,0),maxit=5,debug.kmc=F)){
            del.loc=which(delta==1)[1:2];
            tmp=c(0,0);
            if (del.loc[2]!=2) tmp[2]=sum(as.numeric(delta[1:(del.loc[2]-1)]==0)/( rep(1-re$prob[1],2) ) )
            UD<-cbind(g[[1]](re$times[1:2]),g[[2]](re$times[1:2]))
            uu.lambda=as.vector(
              solve(UD)%*%(n-1/re$prob[del.loc]-tmp)
              )
              # debug oupur lambda: print(uu.lambda)
            uu.lambda
          }
          init.lam=u.lambda2()
        }else{
          u.lambda3<-function(){return(0);}
          init.lam=u.lambda3()
        }
      }
    
  }else{init.lam=rep(0,length(g))}
  
	if (using.num || ( length(g)!=1) ){
	  multiroot(kmc.comb123,start=init.lam,ctol=rtol,useFortran =using.Fortran)$root -> lambda
	  
	}else{
	  multiroot.nr(f_=kmc.comb12,xinit=init.lam,it=15,C=1,FALSE,tol=rtol) -> lambda;
	}	
    if(em.boost & (length(g)==1) ){
		loglik.null<- WKM(kmc.time,delta)$logel
	}else{
        omega.lambda(kmc.time=kmc.time,delta=delta,lambda=0,g=g,gt.mat=gt.mat)->re0; ## set lambda=0, it compute KM-est
        loglik.null<- kmc.el(delta,re0$omega,re0$S)
 	}
    result<-tryCatch(
        result<-omega.lambda(kmc.time,delta,lambda,g,gt.mat=gt.mat)
        ,
        error=function(cond) {
            message(cond)
            return(list(S=NA,omega=NA,gt=NA))
        }
	)
 if (!is.na(result$S[1])){
			loglik.ha<-kmc.el(delta,result$omega,result$S)
			re.tmp <- list(loglik.null=loglik.null,loglik.h0=loglik.ha,"-2llr"=-2*(loglik.ha-loglik.null),g=g,time=x,status=d,phat=result$omega,pvalue=1-qchisq(-2*(loglik.ha-loglik.null),df=length(g)),lambda=lambda);
			if (re.tmp[["-2llr"]]>100) warning('\nThe results may be not feasible!\n');
		}else{
			re.tmp <- list(loglik.null=loglik.null,loglik.h0=NA,"-2llr"=NA,g=g,time=x,status=d,phat=NA,pvalue=NA,df=NA,lambda=NA);
		}
    class(re.tmp)<-"kmcS3";
	return(re.tmp);
}

kmc.bjtest<-function(
  y, d, x, beta,init.st="naive"
){
  
  n <- length(y)
  x <- as.matrix(x)
  xdim <- dim(x)
  if (xdim[1] != n) 
    stop("check dim of x")
  if (length(beta) != xdim[2]) 
    stop("check dim of x and beta")
  e <- y - as.vector(x %*% beta)
  ordere <- order(e, -d)
  esort <- e[ordere]
  dsort <- d[ordere]
  xsort <- as.matrix(x[ordere, ])
  dsort[length(dsort)] <- 1
  temp0 <- WKM(esort, dsort, zc = 1:n)
  pKM <- temp0$jump
  temp <- redistF(y = esort, d = dsort, Fdist = pKM)
  weight <- temp$weight/n
  A <- matrix(0, ncol = xdim[2], nrow = n)
  for (i in 1:n) if (dsort[i] == 1) {
    A[i, ] <- t(as.matrix(weight[1:i, i])) %*% xsort[1:i,]
    A[i, ] <- A[i, ]/pKM[i]
  }
  
  gt.matrix = t(A*esort);
  delta = dsort;
  kmc.time = esort;
  
  kmc.comb123<-function(x){
    kmc_routine4(lambda=x,delta=delta,gtmat=gt.matrix)->re;
    return(re);
  }
  
  u.lambda2<-function(re=el.cen.EM2.kmc(x=kmc.time,d=delta,fun= function(t, q) {
        t * q
    },mu=c(0,0),maxit=5,debug.kmc=F,q=A)){
    del.loc=which(delta==1)[1:2];
    tmp=c(0,0);
    if (del.loc[2]!=2) tmp[2]=sum(as.numeric(delta[1:(del.loc[2]-1)]==0)/( rep(1-re$prob[1],2) ) )
    UD<-cbind(gt.matrix[1:2,1],gt.matrix[1:2,2])
    uu.lambda=as.vector(
      solve(UD)%*%(n-1/re$prob[del.loc]-tmp)
    )
    # debug oupur lambda: print(uu.lambda)
    uu.lambda
  }
  
  if (init.st=="naive"){ init.lam=c(0,0)}else{init.lam=u.lambda2()}
  cat("init.lam:\t",init.lam,'\tLAM')
  multiroot(kmc.comb123,start=init.lam,useFortran = T,rtol = 1e-9, atol = 1e-9, ctol = 1e-9)$root -> lambda
  cat(lambda,'\n')
  omega.lambda(kmc.time=esort,delta=delta,lambda=lambda,g=NULL,gt.mat=gt.matrix) -> result; ## set lambda=0, it compute KM-est  
  temp2<-kmc.el(delta,result$omega,result$S)
  pnew <- result$omega
  logel1 <- temp0$logel
  logel2 <- temp2
  list(prob = pnew, logel = logel1, logel2 = logel2, `-2LLR` = 2 * 
         (logel1 - logel2))
}


plotkmc2D <-function(resultkmc,flist=list(f1=function(x){x},f2=function(x){x^2}),range0=c(0.2,3,20)){
	tmp.df=length(resultkmc$g);
	xx<-resultkmc[["-2llr"]];
	xl<-seq(0,max(6,xx+2),0.01)
	plot(xl,dchisq(xl,df=tmp.df),type='l',main='Kaplan-Meier Estimator with Constraint',xlab="X",ylab="Probabilty");
	points(xx,dchisq(xx,df=tmp.df),col='red',lty=2,type='h');
	
	if (tmp.df==2){
		#X11();
		theta0=-do.call(c,lapply(resultkmc$g,function(x)(x(0))));
		x.grid=seq((theta0[1]-range0[1]),(theta0[1]+range0[1]),length.out=range0[3]);
		y.grid=seq((theta0[2]-range0[2]),(theta0[2]+range0[2]),length.out=range0[3]);
		tmp.z=matrix(0,range0[3],range0[3]);
		for (ii in 1:range0[3]){
			for (jj in 1:range0[3]){
				tmpg=list(f1=function(xuu){flist[[1]](xuu)-tmp1},f2=function(xuu){flist[[2]](xuu)-tmp2});
				tmp1=x.grid[ii];
				tmp2=y.grid[jj];
				tmp.z[ii,jj]=kmc.solve(resultkmc$time,resultkmc$status,tmpg)[[2]]
			}
		}
		contour(x.grid,y.grid,tmp.z)
        points(theta0[1],theta0[2],main='CI',col="red");
	}
	
	par(mfrow=c(1,1))
	return (list(X=x.grid,Y=y.grid,Z=tmp.z));
}

