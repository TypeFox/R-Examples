prateD.exact <-
function(x1,x2,e1,e2,delta.grd, n.grd=15, midp=T)
      {
        fit1=poisson.confint(x1,0.9995); 
        fit2=poisson.confint(x2,0.9995);
        L2=(fit2$upper-fit2$lower)/e2;
        L1=(fit1$upper-fit1$lower)/e1;
        
        rev=F
        if(L2<L1)
          {e3=e2
           x3=x2
           fit3=fit2

           e2=e1
           x2=x1                
           fit2=fit1

           e1=e3
           x1=x3
           fit1=fit3

           delta.grd=rev(-delta.grd)

           rev=T}  

        l1=fit1$lower 
        u1=fit1$upper
        N1=qpois(0.999999, u1)

        lambda1.grd=seq(l1, u1, length=n.grd)/e1; 
        pnull1.tot=matrix(0, N1+1, n.grd)
        for(b in 1:n.grd)
           pnull1.tot[,b]=dpois((0:N1),lambda1.grd[b]*e1)

        u2=max(fit2$upper, e2*(u1/e1+max(delta.grd)))
        N2=qpois(0.999999, u2)

      
        dfnull=matrix(0, N1+1,N2+1)
        sdnull=matrix(0, N1+1,N2+1)
        for(i in 0:N1)
           {lambda1=(i+0.5)/e1
            lambda2=(0:N2+0.5)/e2
            dfnull[i+1, ]=(0:N2)/e2-i/e1
            sdnull[i+1, ]=sqrt(lambda1/e1+lambda2/e2)
           }

        pv1=pv2=numeric(0)
        for(theta in delta.grd)
           {lambda1=(x1+0.5)/e1
            lambda2=(x2+0.5)/e2
            t.stat=(x2/e2-x1/e1-theta)/sqrt(lambda1/e1+lambda2/e2)
            tnull=(dfnull-theta)/sdnull

            pvalue1=pvalue2=rep(0,n.grd); error=1e-6
            
            if(u1/e1+theta<=0)
                    { pv1=c(pv1, 0); pv2=c(pv2, 1) }
            if(u1/e1+theta>0)
              {for(b in 1:n.grd)
                  {lambda1=lambda1.grd[b]
                   lambda2=lambda1+theta
                   if(lambda2>=0)
                     {pnull1=pnull1.tot[,b]
                      pnull2=dpois((0:N2), lambda2*e2)
                      n1.adj=N1+1-max(c(1,(1:(N1+1))[cumsum(sort(pnull1))<error]))
                      n2.adj=N2+1-max(c(1,(1:(N2+1))[cumsum(sort(pnull2))<error]))
                      id1=order(pnull1)
                      id2=order(pnull2)
                      id1=(id1[(N1+1):1])[1:n1.adj]
                      id2=(id2[(N2+1):1])[1:n2.adj]
                      pnull=pnull1[id1]%*%t(pnull2[id2])
                      if(midp==T)
                        {pvalue1[b]=sum(pnull[tnull[id1,id2]>t.stat])+sum(pnull[tnull[id1,id2]==t.stat])*0.5
                         pvalue2[b]=sum(pnull[tnull[id1,id2]<t.stat])+sum(pnull[tnull[id1,id2]==t.stat])*0.5
                        }else{
                              pvalue1[b]=sum(pnull[tnull[id1,id2]>=t.stat])
                              pvalue2[b]=sum(pnull[tnull[id1,id2]<=t.stat])
                              }
                      }
    
                
              }
             pv1=c(pv1, max(pvalue1)+(1-0.9995))
             pv2=c(pv2, max(pvalue2)+(1-0.9995))
           }
          }
        
        if(rev==T)
          {pv3=pv2
           pv2=rev(pv1)
           pv1=rev(pv3)
           }
          
        return(list(pv1=pv1, pv2=pv2))
      }
