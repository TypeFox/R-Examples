B2R <-
function(B,n) {

      T=nrow(B)    
      tau1=ncol(B)
  
    if (nargs()<2)
       { side=tau1
         } else {
        side=min(n,tau1)}
  
      nT=floor((side-1)/T)+1
   
      Bex=matrix(0,(T*nT),tau1)
      Bex[1:T,]=B
   
    if (nT>1) 
       { for (n in 2:nT) { Bex[((n-1)*T+1):(n*T),]=B }  }
  
      diags1<-list(Bex[1:side,1])    
      R = Matrix::bandSparse(length(Bex[1:side,1]), k = 0, diag = diags1)
      R=as.matrix(R)
   
    if (side>1)
       {  for (tau in 1:(side-1))
         { diags2<-list(Bex[1:(side-tau),(tau+1)])
           MR=Matrix::bandSparse(nrow(R), k = tau, symmetric = TRUE, diagonals = diags2)
           MR=as.matrix(MR)
           R=R+MR      
         }
       }
     result = list(R=R) 
     class(result) = "B2R"
     result
}

