RSA <- function(Data,N_subset,Stepscale,Total_Iteration,Warm)
{  
  DataCol = ncol(Data);
  DataNum = nrow(Data);
  beta = rep(0,DataCol-2);
  phi = 0;
  sigmasq = 0;
  tausq = 0;
  RetVec2 = .C("RSAc",
               as.numeric(Data),
               as.integer(DataCol),
               as.integer(DataNum),
               as.integer(N_subset),
               as.integer(Stepscale),
               as.integer(Total_Iteration),
               as.integer(Warm),             
               as.numeric(beta),
               as.numeric(phi),
               as.numeric(sigmasq),
               as.numeric(tausq)
               )
  beta = RetVec2[[8]];
  phi = RetVec2[[9]];
  sigmasq = RetVec2[[10]];
  tausq = RetVec2[[11]];
  Z = NULL
  z = list(beta = beta,phi=phi,sigmasq=sigmasq,tausq=tausq)
  #return(c(beta0,beta1,phi,sigmasq,tausq));
  z
}
