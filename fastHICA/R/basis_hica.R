basis_hica <-
function(X,maxlev=dim(X)[2]-1,dim.subset=512){
 Xin=X
 n <- dim(X)[1]
 p <- dim(X)[2]
 if (maxlev<1 || maxlev>(p-1)){
  stop(paste("the maximum level must to be between 1 and ",p-1,sep=""))
 }
 Atot <- diag(rep(1,p))
 Aparz <- Atot
 basis <- vector("list",maxlev)
 sim_hica <- similarity_hica(X,dim.subset)
 stat_value <- sim_hica$similarity_matrix
 sam <- sim_hica$subset
 ind=NULL
 for (h in 1:maxlev){
  print(paste("Performing the step ",h," out of ",maxlev,sep=""))
  index1 <- stat_value[which.max(stat_value[,1]),2]
  index2 <- stat_value[which.max(stat_value[,1]),3]
  if (max(stat_value[,1])<0.99999){
   p_sdev <- princomp(cbind(X[,index1],X[,index2]))$sdev
   ica <- fastICA(cbind(X[,index1],X[,index2]),2,w.init=cbind(c(p_sdev[1],0),c(0,p_sdev[2])))
   A <- t(ica$A)
   val1 <- sum(A[,1]^2)
   val2 <- sum(A[,2]^2)
  }
  if (max(stat_value[,1])>=0.99999){
   pca <- princomp(cbind(X[,index1],X[,index2]))
   A <- pca$loa
   val1 <- 50
   val2 <- 10
  }
  A <- cbind(A[,1]/sqrt(sum(A[,1]^2)),A[,2]/sqrt(sum(A[,2]^2)))
  Aparz[c(index1,index2),c(index1,index2)] <- A
  Atot <- Atot%*%Aparz
  basis[[h]] <- Atot
  if(val1<val2){
   ind <- rbind(ind,c(index2,index1,max(stat_value[,1])))
   X[,index1] <- X[,index1]-mean(X[,index1])
   X[,index2] <- X[,index2]-mean(X[,index2])
   X <- X%*%solve(t(Aparz))
   Aparz <- diag(rep(1,dim(Aparz)[2]))
   del <- NULL
   for (cc in 1:dim(stat_value)[1]){
    if (stat_value[cc,2]==index1 || stat_value[cc,2]==index2 || stat_value[cc,3]==index1 || stat_value[cc,3]==index2){
     del <- c(del,cc)
    }
   }
   stat_value <- stat_value[-del,]
   for (f in 1:dim(X)[2]){
    if (f != index2 && sum(f==ind[,2])==0){
     aaa <- dcor(X[sam,f],X[sam,index2])
     stat_value <- rbind(stat_value,c(aaa,f,index2))
    }
   }
  }
  if(val1>val2){
   ind <- rbind(ind,c(index1,index2,max(stat_value[,1])))
   X[,index1] <- X[,index1]-mean(X[,index1])
   X[,index2] <- X[,index2]-mean(X[,index2])
   X <- X%*%solve(t(Aparz))
   Aparz <- diag(rep(1,dim(Aparz)[2]))
   del <- NULL
   for (cc in 1:dim(stat_value)[1]){
    if (stat_value[cc,2]==index1 || stat_value[cc,2]==index2 || stat_value[cc,3]==index1 || stat_value[cc,3]==index2){
     del <- c(del,cc)
    }
   }
   stat_value <- stat_value[-del,]
   for (f in 1:dim(X)[2]){
    if (f != index1 && sum(f==ind[,2])==0){
     aaa <- dcor(X[sam,f],X[sam,index1])
     stat_value <- rbind(stat_value,c(aaa,f,index1))
    }
   }
  }
 }
 list(X=Xin, basis=basis, aggregation=ind)
}
