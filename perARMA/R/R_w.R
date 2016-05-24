R_w <-
function(theta,phi,nlen){

       T=nrow(theta)
       q=ncol(theta) 
  
    if (length(phi)>1)                         
       {T=nrow(phi)
       p=ncol(phi)   
       } else {
        p=0}                             
                                          
     if (p==0)  {R_w_ma<-R_w_ma(theta,1,nlen)
                 R=R_w_ma$R} 
                                           
    del=theta[,1]
    del2=del*del                               
    q=q-1                                       
    lag=max(p,q)                                
    side=(lag+1)*T                              
 
       d=matrix(1,side,1)
       A = Matrix::bandSparse(side, k = 0, diag = list(d))
       A= as.matrix(A)

    for (tau in 0:lag)
    {for (s in 1:T)                         
       { row=s+tau*T
         for (j in 1:p)                                                  
            {i1=s-j
             i2=s-tau
             if (i1 < i2)                     
                {itemp=i2
                i2=i1
                i1=itemp}
            
            good1=matlab::mod(i1-1,T)+1       
            nT=good1-i1                     
            good2=i2+nT                       
            idiff=good1-good2
            col=good1+T*idiff
            A[row,col]=-phi[s,j]}
        }
      }
        
 
    B_x_xi<-matrix(0,T,lag+1) 
    B_x_xi[,1]=theta[,1]                           
    Rold=theta[,1]  

        if (q>1)
          {   for (tau in 1:q) 
          {d=matrix(1,T*(tau+1),1)
           B = Matrix::bandSparse(T*(tau+1), k = 0, diag = list(d))
           B= as.matrix(B)
    
       for (s in 1:T)                             
        { row=s+tau*T
           for (j in 1:p)                      
            {i1=s-j
             i2=s-tau
              if (i1 >= i2)                    
                {good1=matlab::mod(i1-1,T)+1    
                nT=good1-i1                     
                good2=i2+nT                     
                idiff=good1-good2
                col=good1+T*idiff
                B[row,col]=phi[s,j]            
                }
             }
         }
  

      if (tau > (ncol(theta)-1) )                    
        {Rnew=B%*%c(Rold,matrix(0,T,1))    
         }  else  {
        Rnew=B%*%c(Rold,theta[,tau+1])
        }
    
        B_x_xi[,tau+1]=Rnew[(tau*T+1):((tau+1)*T)]
        Rold=Rnew                                                
        }
}

       Yvec<-matrix(0,side,1)
            minlen=q+1
       for (t in 1:T) {Yvec[t,1]=theta[t,1:minlen]%*%B_x_xi[t,1:minlen]}

       for (tau in 1:q)                                 
          { for (s in 1:T)                               
          {row=s+tau*T
           rhs=0
            for (j in 1:minlen)                    
            {i1=s-tau
             i2=s-j+1                            
              if (i1 >= i2)                       
                {good1=matlab::mod(i1-1,T)+1      
                 nT=good1-i1                      
                 good2=i2+nT                     
                 idiff=good1-good2
                 col=good1+T*idiff
                 rhs=rhs+theta[s,j]%*%B_x_xi[good1,idiff+1]}
           }
        Yvec[row]=rhs
      } 
    }

     gamma=qr.solve(A)%*%Yvec

      B<-matrix(0,T,lag+1)
      B[,1]=gamma[1:T]
        for (tau in 1:lag) {B[,tau+1]=gamma[(tau*T+1):((tau+1)*T)]}
                                                
        diags<-list(B[,1])
        R_1 = Matrix::bandSparse(length(B[,1]), k = 0, diag = diags)
        R_1=as.matrix(R_1)

                                                
      for (tau in 1:lag)
        { diags<-list(B[1:(T-tau),(tau+1)])
          MR_1=Matrix::bandSparse(nrow(R_1), k = tau, diag = diags)
          MR_1=as.matrix(MR_1)
          NR_1=Matrix::bandSparse(nrow(R_1), k = -tau, diag = diags)
          NR_1=as.matrix(NR_1)
          R_1 = R_1+MR_1
          R_1 = R_1+NR_1
    }
          R_1=R_1[1:lag,1:lag]


         R_11<-matrix(0,lag,lag)
            for (t in 1:lag)                                   
          { for (s in 1:t)                         
           {ss=s+lag                             
            rhs=0                                 
            for (k in 1:(q+1))                      
            { i1=t
              i2=ss-k+1                                      
              if (i1>=i2)                                    
               { good1=matlab::mod(i1-1,T)+1                
                nT=good1-i1                                 
                good2=i2+nT                                 
                idiff=good1-good2
                rhs=rhs+theta[s,k]%*%B_x_xi[good1,idiff+1]  
               }
             }
        R_11[t,s]=rhs
       } 
    }     
              
       R_w_ma<-R_w_ma(theta,lag+1,nlen-lag)
       R_22=R_w_ma$R
       rindex=R_w_ma$rindex

    result = list( R_1= R_1, R_22=R_22, rindex=rindex) 
    class(result) = "R_w"
    result
   
}