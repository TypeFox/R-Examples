ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
     dnx <- dimnames(X)
     if (is.null(dnx)) 
         dnx <- vector("list", 2)
     s <- svd(X)
     nz <- s$d > tol * s$d[1]
     structure(if (any(nz)) 
         s$v[, nz] %*% (t(s$u[, nz])/s$d[nz])
         else X, dimnames = dnx[2:1])
}

c1_weight <- function(W, T, X, d, dt = 1, order) {
     inv <- ginv(upinfor(W, T, X, order))
     f1 <- c1_weight_1(W,T,X,inv,order)
     f2 <- c1_weight_2(W,T,X,inv,order)
     W - d * (f1 %*% ginv(f2))
}

c_weight <-function(W,T,X,d,dt,order) {
     inv<-ginv(upinfor(W,T,X,order))
     f1<-c_weight_1(W,T,X,inv,dt,order)
     f2<-c_weight_2(W,T,X,inv,dt,order)
     W-d*(f1%*%ginv(f2))
}

D_weight <- function(W,T,X,d,dt=1,order) {
     inv<-ginv(upinfor(W,T,X,order))
     f1<-D_weight_1(W,T,X,inv,order)
     f2<-D_weight_2(W,T,X,inv,order)
     W-d*(f1%*%ginv(f2))
}

M_weight <- function(W,T,X,d,dt,order,lambda) {
     inv<-ginv(upinfor(W,T,X,order))
     f1<-M_weight_1(W,T,X,inv,dt,order,lambda)
     f2<-M_weight_2(W,T,X,inv,dt,order,lambda)
     W-d*(f1%*%ginv(f2))
}

search_weight<-function(X,T,epsilon_w,dt,order,f,...) {
     diff <- 10
     W <- rep(1/length(X),length(X)-1)

     while(diff>epsilon_w) {
           d <- 0.2
           NW <- f(W,T,X,d,dt,order,...)
           minW <- min(min(NW),1-sum(NW))
           while(minW<0 & d>.0001) {
                d <-d/2
                NW <- f(W,T,X,d,dt,order,...)
                minW <- min(min(NW),1-sum(NW))
           }
           NW <- c(NW,1-sum(NW))
           n <- length(NW)
           minW <- min(NW)
           diff <- max(abs(W-NW[1:n-1]))
           if (abs(minW)<.0000001||minW<0) {
                for(i in 1:n) 
                     if (NW[i]==minW) NW[i] <- 0          
           }
           D <- rbind(X,NW)
           for (i in 1:n) 
                if (D[2,i]==0) D[,i] <- NA
           X <- na.omit(D[1,])
           W <- na.omit(D[2,])
           W <- W[1:length(X)-1]
     }
     W <- c(W,1-sum(W))
     D <- rbind(X,W)
     D
}


