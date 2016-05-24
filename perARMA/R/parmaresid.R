parmaresid <-
function(x,stype,del,phi,...)
{
  parmaresid_full <-function(x,stype,del,phi,theta)

 {      
         if (length(phi)>1)  
           {  Tp=nrow(phi)
              p=ncol(phi)
            } else {
              Tp=0
              p=0}

         if (!is.null(theta))  
           {  Tq=nrow(theta)
              q=ncol(theta)
            } else {
              Tq=0
              q=0}

      ZTHRS=10*.Machine$double.eps
      DISPRANK=1

      nx=length(x)
      T=length(del)
     
      m=max(p,q)  

      if (p)
      {A=cbind(matrix(1,T,1),-phi)
       } else {
       A=matrix(1,T,1)}
  
       B=matrix(1,T,1)

       m=max(p,q)
         
       pf<-parmafil(A,B,x)
       w0<-pf$y 

   if (stype)
         {  if (length(theta)>1)
            {  R_w<-R_w(cbind(del,theta),phi,nx)  
             } else {
               R_w<-R_w(del,phi,nx)  
             }        
            RW1=R_w$R_1
            RW2=R_w$R_22
            rindex=R_w$rindex   
           
            Rc=cbind(rindex,RW2)

            r<-matrix(0,nx,nx)
            r[1:p,1:p]=RW1  
     
            for(i in 1:nrow(Rc))
            {r[Rc[i,1]+p,Rc[i,2]+p]=Rc[i,3] }
            r=as.matrix(r)
     
           decomp<-qr(r)   
           R1=qr.R(decomp)              
           dr1=diag(R1) 
             
           igood<-NULL
           for( i in 1:length(dr1))
               {if(abs(dr1[i])>ZTHRS) 
               {igood=c(igood,i) }} 
                ngood=length(igood) 


           rm=r[1:m,1:m]
           e_m<-eigen(rm,only.values = TRUE)
           eig=e_m$values
            eigsum=base::sum(eig[which(eig>0)])           
           if (eigsum<m)
            {igood[1:m]=0}

                 R=chol(r[igood,igood])
                 L=t(R)
                 Linv<-qr.solve(L)
                 e=Linv%*%w0[igood]
     } else {
         
         if (length(theta)>1)   
         {  RW_ma<-R_w_ma(cbind(del,theta),p+1,nx-p)   
         } else {
           RW_ma<-R_w_ma(del,p+1,nx-p)  
         }        
          RW=RW_ma$R
          rindex=RW_ma$rindex
          r<-matrix(0,nx-p,nx-p)

           Rc=cbind(rindex,RW)
     
            for(i in 1:nrow(Rc))
            {r[Rc[i,1],Rc[i,2]]=Rc[i,3] }
            r=as.matrix(r)

          w0=w0[(p+1):nx]       
                                    
          dr=diag(r)
          igood_r<-NULL
          for( i in 1:length(dr))
          {if(abs(dr[i])>2.2204e-015) 
           {igood_r=c(igood_r,i) }}  
                               
            r1=r[igood_r,igood_r]

           decomp<-qr(r1)
           R1=qr.R(decomp)
           dr1=diag(R1)

          igood_r1<-NULL
          for( i in 1:length(dr1))
                   {if(abs(dr1[i])>ZTHRS) 
                   {igood_r1=c(igood_r1,i) }} 
         
         R=r1[igood_r1,igood_r1]

         R=chol(r1[igood_r1,igood_r1])
         L=t(R)
         Linv<-qr.solve(L)
         w0_r=w0[igood_r]
         w0_r1=w0_r[igood_r1]
         e=Linv%*%w0_r1
     }     
          resids=e

         result = list(resids=resids)
         class(result) = "parmaresid"
         result 

}
   L<-modifyList(list(theta=NULL),list(x = x, stype=stype, phi=phi,del=del,...))
    do.call(parmaresid_full,L)

}
