Deff<-function(weight,P,dose,LB,UB,r=10,grid=0.01,epsilon=.001,epsilon_w=10^-6)
{
	np <- length(P)
      if(np==2) {T <- c(1,-P[2],P[2]*log(P[1]),0); X<-runif(2,LB,UB)}
      if(np==3) {T <- c(P[1],-P[3],P[3]*log(P[2]),0); X<-runif(3,LB,UB)}
      if(np==4) {T <- c(P[2],-P[4],P[4]*log(P[3]),P[1]); X <- c(LB,LB+(UB-LB)/3,LB+2*(UB-LB)/3,UB)}

      LB <- round(LB,2)
      UB <- round(UB,2)

      k <- length(T)

      W <- rep(1/length(X),length(X)-1)

      M <- upinfor(W,T,X,np)
      n <- 1
      while(n<r){
           x <- seq(LB,UB,grid)
           n1 <- length(x)
           ds <- rep(0,n1)
           inv <- ginv(M)
           for (i in 1:n1) 
                 ds[i] <- t(f234(T,x[i],np))%*%inv%*%f234(T,x[i],np)
           newX <- x[which(ds==max(ds))]
           newX <- round(newX[1],2)
           an <- 1/(n+1)
           p <- abs(max(ds)-np)
           newM <- (1-an)*M+an*f234(T,newX,np)%*%t(f234(T,newX,np))
           M <- newM
           X <- c(X,newX)
           n <- n+1
      }

      X <- unique(X[(length(X)-np):length(X)])
      X <- sort(X,decreasing=F)
      if (length(X)==1) X <- c(X,runif(1,LB,UB))

      it <- 1
      while(p>epsilon) {
           x <- seq(LB,UB,grid)
           n1 <- length(x)
           ds <- rep(0,n1)
           D <- search_weight(X,T,epsilon_w,dt=1,np,D_weight)
           X <- D[1,]
           W <- D[2,1:length(X)-1]
           inv <- ginv(upinfor(W,T,X,np))
           for (i in 1:n1) 
                ds[i] <- t(f234(T,x[i],np))%*%inv%*%f234(T,x[i],np)
           newX <- x[which(ds==max(ds))]
           newX <- round(newX[1],2)
           X <- c(X,newX)
           X <- sort(X,decreasing=F)
           X <- unique(X)
           newp <- abs(max(ds)-np)
           if(abs(newp-p)<.0000000001) newp <- 10^-20
           if(it>20) newp <- 10^-20
           p <- newp
           it <- it+1
      }

      X <- D[1,]
      W <- D[2,1:length(X)-1]
      M <- upinfor(W,T,X,np)
      x <- seq(LB,UB,grid)
      ds <- rep(0,length(x))
      for (i in 1:length(x)) 
           ds[i] <- t(f234(T,x[i],np))%*%ginv(M)%*%f234(T,x[i],np)-np

      weight <- weight[1:length(dose)-1]
      eff <- (det(upinfor(weight,T,dose,np))/det(upinfor(W,T,X,np)))^(1/np)
    
      R=new("PAR",fid="Deff",LB=LB,UB=UB,grid=grid,ds=ds)
      return(new("OPT",Par=R,Opt=D,Eff=eff))
}



