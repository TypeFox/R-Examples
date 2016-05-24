parmaf <-
function(x,T,p,q,af,bf,...)
{
  parmaf_full<-function(x,T,p,q,af,bf,a0,b0,stype)
    {  
           if (!is.null(b0))  
           { del_mask=as.numeric(b0[,1]!=0)
            } else {
             del_mask=matrix(1,T,1)
            }

           if ( q != (ncol(bf)-1) ) { stop("number of columns of bf not equal to q+1")}
           if ( (ncol(af))!=p )  { stop("number of columns of af not equal to p")}

             nx=length(x)
             naf=base::sum(af[,]!=0)
             nbf=base::sum(bf[,]!=0)
            
             iaf=which(af[,]!=0)
             ibf=which(bf[,]!=0)

          if (nargs() < 7 | is.null(a0) | is.null(b0) )            
          {  
           if (p)  
           { estimators<-perYW(x,T,p,NaN)
             phi0=estimators$phi
             del0=estimators$del
             phi0=as.matrix(phi0)
             del0=as.matrix(del0)

             phtha<-phth2ab(phi0)
             a0=phtha$a
             a0=as.matrix(a0)

             phthb<-phth2ab(del0)
             b1=phthb$a
             b1=as.matrix(b1)

             a00<-matrix()
             for (i in 1:p)  {a00[seq(1,T)+(i-1)*T]=a0[,i]} 
              } else {
             phi0<-matrix()
             del0=matrix(1,T,1)

             phthb<-phth2ab(del0)
             b1=phthb$a
             b1=as.matrix(b1)
             a00<-matrix()
            }
          }

             b0=cbind(b1,matrix(0,T,q))
           

        if ( length(which(a00!=0)))
        {  

           if (p)
           {  estimators<-perYW(x,T,p,NaN)
              phi0=estimators$phi
              del0=estimators$del
              phi0=as.matrix(phi0)
              del0=as.matrix(del0)

              phtha<-phth2ab(phi0)
              a0=phtha$a
              a0=as.matrix(a0)

              del0=del0*del_mask

              phthb<-phth2ab(del0)
              b1=phthb$a
              b1=as.matrix(b1)

              b0=cbind(b1,matrix(0,T,q))

             a00<-matrix()
             for (i in 1:p) { a00[seq(1,T)+(i-1)*T]=a0[,i]}  
             a00=a00[iaf]

                } else {
              phi0=matrix(0,T,1) 
              del0=matrix(1,T,1) 

               del0=del0*del_mask

              phthb<-phth2ab(del0)
              b1=phthb$a
              b1=as.matrix(b1)

              b0=cbind(b1,matrix(0,T,q))
              }
              b00<-matrix()
              for (i in 1:(q+1)) { b00[seq(1,T)+(i-1)*T]=b0[,i]}
                b00=b00[ibf]

              ab0=c(a00,b00)
           }  else {
             
              b00<-matrix()
              for (i in 1:(q+1)) { b00[seq(1,T)+(i-1)*T]=b0[,i]}
                b00=b00[ibf]

              ab0=c(b00)
             
          }
      

  conpars=c(T,p,q,naf,nbf,del_mask,iaf,ibf,stype)

   fun<-function(ab) 
   { val<-loglikef(ab,x,conpars)
     val=val$y
     val=as.numeric(val)
     val}
   ans<-optim(ab0,fun, gr=NULL, method="BFGS", lower=-Inf, upper=Inf, hessian=FALSE, control=list())

   ab=ans$par
   ab=as.matrix(ab)
   negloglik=ans$val

         
        endab=length(ab)
        if (p)
          { a<-matrix(0,T,p)
            a[iaf]=ab[1:naf]

            ab2p<-ab2phth(a)
            phi=ab2p$phi
            phi=as.matrix(phi)

               } else {
            a<-matrix()
            phi<-matrix() }
 
          b<-matrix(0,T,q+1) 
          b[ibf]=ab[(naf+1):endab]


          b1<-matrix(0,T,1)
          b1=b[,1]

          ab2pb1<-ab2phth(b1)
          b1ab=ab2pb1$phi

          del=b1ab*del_mask


         if (q)
          { b2<-matrix(0,T,q)
            b2=b[,2:(q+1)]
            ab2ptheta<-ab2phth(b2)
            theta=ab2ptheta$phi
            theta=as.matrix(theta)
              } else {
            theta<-matrix()}


      x=as.matrix(x)
      res<-parmaresid(x,stype,del,phi,theta)

      resids=res$resids

      mse=(t(resids)%*%resids)/length(resids)
      mse=as.matrix(mse)
                       
      netpars=naf+nbf
      nval=(nx-p)
     

      aicval=negloglik + 2*netpars
      bicval=negloglik + netpars*log(nval)                          
      fpeval=mse*(1+netpars/nval)/(1-netpars/nval)

      result = list(a=a,b=b,negloglik=negloglik,aicval=aicval,fpeval=fpeval,bicval=bicval,resids=resids)
      class(result) = "parmaf"
      result
     }

    L<-modifyList(list(a0=NULL,b0=NULL, stype=0),list(x = x, T=T, p=p,q=q, af=af, bf=bf,...))
    do.call(parmaf_full,L)

}



