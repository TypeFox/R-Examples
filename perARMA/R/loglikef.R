loglikef <-
function(ab,x,conpars)
{   
      ZTHRS=10*.Machine$double.eps
      DISPRANK=0
      T=as.numeric(conpars[1])
      p=as.numeric(conpars[2])
      q=as.numeric(conpars[3])
      naf=conpars[4]
      nbf=conpars[5]
      del_mask=conpars[6:(6+T-1)]
      iaf=conpars[(6+T):(6+T-1+naf)]
      end=length(conpars)
      ibf=conpars[(6+T-1+naf+1):(end-1)]
      stype=conpars[end]

      nx=length(x)

     if (p>0)
       {  a=matrix(0,T,p)
          a[iaf]=ab[1:naf]
          phiab<-ab2phth(a)
          phi=phiab$phi
           } else {
          phi<-matrix()  }
 
        b=matrix(0,T,q+1)
        end=length(ab)
        
       b[ibf]=ab[(naf+1):end]

       b1a=as.matrix(b[,1])
       b1ab<-ab2phth(b1a)
       b1a=b1ab$phi

       del=b1a*del_mask 
   
     if (q>0)
      {  b2a=as.matrix(b[,2:(q+1)])
         thetaab<-ab2phth(b2a)
         theta=thetaab$phi
           } else {
         theta<-matrix()  }

        m=max(p,q)

     if(p) {
      A=cbind(matrix(1,T,1),-phi)
       }  else {
      A=matrix(1,T,1)
       }
     B=matrix(1,T,1)      
     x=as.matrix(x)   
            
         pf<-parmafil(A,B,x)
         w0=pf$y

  
    if (stype)
         {   if (length(theta)>1)
            {  R_w<-R_w(cbind(del,theta),phi,nx)  
             } else {
               R_w<-R_w(del,phi,nx)  
             } 
                  
            RW1=R_w$R_1
            RW2=R_w$R_22
            rindex=R_w$rindex       

            r<-matrix(0,nx,nx)
            r[1:p,1:p]=RW1  
             

            Rc=cbind(rindex,RW2)

            
          
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
           l=length(diag(L))
           diagL=diag(L)
           logdiag<-matrix(0,1,l)
           for (i in 1:l)  { logdiag[i]=log(diagL[i])}
             logdetL=base::sum(logdiag)  
 
           zres=t(e)%*%e 
           mse=zres/length(e)

           y=.5*length(e)*log(2*pi) +logdetL + .5*(zres)
           

         result = list(y=y)
         class(result) = "loglikef"
         result 
}