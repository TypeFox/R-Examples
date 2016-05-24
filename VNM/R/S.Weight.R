S.Weight<-function(X,P,lambda,delta,epsilon_w=10^-6)
{
	np <- length(P)
	if(np==2) T <- c(1,-P[2],P[2]*log(P[1]),0)
      if(np==3) T <- c(P[1],-P[3],P[3]*log(P[2]),0)
      if(np==4) T <- c(P[2],-P[4],P[4]*log(P[3]),P[1])

      D <- search_weight(X,T,epsilon_w,delta,np,M_weight,lambda) 
      W <- D[2,1:ncol(D)-1]
      X <- D[1,]
      inv <- ginv(upinfor(W,T,X,np))
      p <- length(W)
      f1 <- rep(0,p)
      f2 <- matrix(c(rep(f1,p)),nrow=p,ncol=p,byrow=F)

      for (i in 1:p) 
           f1[i] <- lambda[1]*d1(T,X[i],X[ncol(D)],inv,np)/np-lambda[2]*d2(T,X[i],X[ncol(D)],inv,np)-(1-lambda[1]-lambda[2])*d3(T,X[i],X[ncol(D)],inv,delta,np)
      
      for(i in 1:p) 
           for(j in 1:p)
                 f2[i,j] <- lambda[1]*dd1(T,X[i],X[j],X[ncol(D)],inv,np)/np-lambda[2]*dd2(T,X[i],X[j],X[ncol(D)],inv,np)-(1-lambda[1]-lambda[2])*dd3(T,X[i],X[j],X[ncol(D)],inv,delta,np)
     
      return(new("SW",Opt.W=D,First.C=f1,Second.C=diag(f2)))
}
      

