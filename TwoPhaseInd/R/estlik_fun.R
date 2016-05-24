

###################################################################
##### This is a collection of functions written for compute   #####
##### point and variance estimate in estimated likelihood     #####
##### method by Pepe 91, in a two-stage study with biased     #####
##### sampling, July 2006                                     #####
###################################################################


scoreMat <- function(X1,X,y1,x1,y,pzw,beta,Nyx){
   fit1  <- drop(exp(X1 %*% beta)/(1+ exp(X1 %*% beta)))
   pfit1 <- ifelse(y1==1,fit1,1-fit1)
   pyxz<- pfit1*pzw
   pyx <- rep(0,length(y1))
   pyx[x1==0 & y1==0] <- sum(pyxz[x1==0 & y1==0])
   pyx[x1==1 & y1==0] <- sum(pyxz[x1==1 & y1==0])
   pyx[x1==0 & y1==1] <- sum(pyxz[x1==0 & y1==1])
   pyx[x1==1 & y1==1] <- sum(pyxz[x1==1 & y1==1])  
   S  <- (y1-fit1)*X1
   B  <- S*pfit1*pzw
   V  <- pfit1*(1-pfit1)
   der.B <- t(S[y1==1 &x1==1,])%*% ((V*X1*pzw/pyx)[y1==1&x1==1,]) - t((X1*V*pfit1*pzw/pyx)[y1==1&x1==1,]) %*% (X1[y1==1&x1==1,])
   der.A <- apply((X1*V*pzw)[y1==1 & x1==1,], 2, sum)
   B1    <- apply(B[y1==1 & x1==1,],2,sum)
   inf11  <- der.B -outer(der.A,B1/((unique(pyx[y1==1 & x1==1]))^2))

   der.B <- t(S[y1==1 &x1==0,])%*% ((V*X1*pzw/pyx)[y1==1&x1==0,]) - t((X1*V*pfit1*pzw/pyx)[y1==1&x1==0,]) %*% (X1[y1==1&x1==0,])
   der.A <- apply((X1*V*pzw)[y1==1 & x1==0,], 2, sum)
   B1    <- apply(B[y1==1 & x1==0,],2,sum)
   inf10  <- der.B -outer(der.A,B1/((unique(pyx[y1==1 & x1==0]))^2))

   der.B <- (-t(S[y1==0 & x1==1,])%*% ((V*X1*pzw/pyx)[y1==0 & x1==1,]) - t((X1*V*pfit1*pzw/pyx)[y1==0 & x1==1,]) %*% (X1[y1==0 & x1==1,])) 
   der.A <- (-apply((X1*V*pzw)[y1==0 & x1==1,], 2, sum))
   B0    <- apply(B[y1==0 & x1==1,],2,sum)
   inf01 <- der.B -outer(der.A,B0/((unique(pyx[y1==0 & x1==1]))^2))

   der.B <- (-t(S[y1==0 & x1==0,])%*% ((V*X1*pzw/pyx)[y1==0 & x1==0,]) - t((X1*V*pfit1*pzw/pyx)[y1==0 & x1==0,]) %*% (X1[y1==0 & x1==0,])) 
   der.A <- (-apply((X1*V*pzw)[y1==0 & x1==0,], 2, sum))
   B0    <- apply(B[y1==0 & x1==0,],2,sum)
   inf00 <- der.B -outer(der.A,B0/((unique(pyx[y1==0 & x1==0]))^2))

   fit  <- drop(exp(X %*% beta)/(1+ exp(X %*% beta)))
   v    <- fit*(1-fit)
   vinf <- t(X)%*% (v*X)
   cheese <- (-inf11*unique(Nyx[y1==1 & x1==1]) - inf10*unique(Nyx[y1==1 & x1==0]) -inf01*unique(Nyx[y1==0 & x1==1]) - inf00*unique(Nyx[y1==0 & x1==0])  + vinf)
   
   score  <- apply(rbind((y-fit)*X,Nyx*B/pyx),2,sum)   
   return(list(cheese=cheese,score=score))
}


##### The extra variance drawn from estimated density when using independence assumption #####


VarDensity <- function(beta,X1,y1,x1,pzw,Nyx,z.from.y,n1,n0,n, N1, N0, N){
   fit1  <- drop(exp(X1 %*% beta)/(1+ exp(X1 %*% beta)))
   pfit1 <- ifelse(y1==1,fit1,1-fit1)
   D1 <- ifelse(y1==1,pfit1*(1-pfit1),-pfit1*(1-pfit1))*X1
   DZW1 <- D1*pzw
   D2 <- matrix(0,length(y1),ncol(X1))
   D2[x1==0 & y1==0,] <- matrix(rep(apply(DZW1[x1==0 & y1==0,],2,sum),n),nrow=n,byrow=T)
   D2[x1==1 & y1==0,] <- matrix(rep(apply(DZW1[x1==1 & y1==0,],2,sum),n),nrow=n,byrow=T)
   D2[x1==0 & y1==1,] <- matrix(rep(apply(DZW1[x1==0 & y1==1,],2,sum),n),nrow=n,byrow=T)
   D2[x1==1 & y1==1,] <- matrix(rep(apply(DZW1[x1==1 & y1==1,],2,sum),n),nrow=n,byrow=T)
   pyx <- rep(0,length(y1))
   pyxz<- pfit1*pzw
   pyx[x1==0 & y1==0] <- sum(pyxz[x1==0 & y1==0])
   pyx[x1==1 & y1==0] <- sum(pyxz[x1==1 & y1==0])
   pyx[x1==0 & y1==1] <- sum(pyxz[x1==0 & y1==1])
   pyx[x1==1 & y1==1] <- sum(pyxz[x1==1 & y1==1])
   weight <- ifelse(z.from.y==1,1/n1,1/n0)   
   W <- (D1/pyx - D2/(pyx)^2 *pfit1)*Nyx * weight
   Wbar <- matrix(0,n,ncol(X1))
   num <- rep(1:n,4)
   for (i in 1:n) Wbar[i,] <- apply(W[num==i,],2,sum)         
   from.y <- z.from.y[1:n]
   I <- var(Wbar[from.y==1,])*(N1/N)^2*n1 + var(Wbar[from.y==0,])*(N0/N)^2*n0
   I
}

