priskD.exact <-
function(x1,x2,n1,n2,delta.grd, n.grd=15, midp=T)
      {
        fit1=binom.confint(x1,n1,0.9995,method="exact"); 
        fit2=binom.confint(x2,n2,0.9995,method="exact");

        L2=(fit2$upper-fit2$lower);
        L1=(fit1$upper-fit1$lower);
        
        rev=F
        if(L2<L1)
          {n3=n2
           x3=x2
           fit3=fit2

           n2=n1
           x2=x1                
           fit2=fit1

           n1=n3
           x1=x3
           fit1=fit3

           delta.grd=rev(-delta.grd)

           rev=T}  

        l=fit1$lower; 
        u=fit1$upper


        p1.grd=seq(l, u, length=n.grd); 
        pnull1.tot=matrix(0, n1+1, n.grd)
        for(b in 1:n.grd)
           {p1=p1.grd[b]
            pnull1.tot[,b]=dbinom((0:n1),n1,p1)}
        dfnull=matrix(0, n1+1,n2+1)
        sdnull=matrix(0, n1+1,n2+1)
        for(i in 0:n1)
           {p1=(i+0.5)/(n1+1)
            p2=(0:n2+0.5)/(n2+1)
            dfnull[i+1, ]=(0:n2)/n2-i/n1
            sdnull[i+1, ]=sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
           }
        pv1=pv2=numeric(0)
        for(theta in delta.grd)
           {p1=(x1+0.5)/(n1+1)
            p2=(x2+0.5)/(n2+1)
            t=(x2/n2-x1/n1-theta)/sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
            tnull=(dfnull-theta)/sdnull
            pvalue1=pvalue2=rep(0,n.grd); error=1e-6
            if(min(p1.grd)+theta>=1)
                    { pv1=c(pv1, 1)
                      pv2=c(pv2, 0) }
            if(max(p1.grd)+theta<=0)
                    { pv1=c(pv1, 0)
                      pv2=c(pv2, 1) }
            if(min(p1.grd)+theta<1 && max(p1.grd)+theta>0)
            {
            for(b in 1:n.grd)
               {p1=p1.grd[b]
                p2=p1+theta
                if(p2>=0 && p2<=1)
                  {pnull1=pnull1.tot[,b]
                   pnull2=dbinom((0:n2), n2, p2)
                   n1.adj=n1+1-max(c(1,(1:(n1+1))[cumsum(sort(pnull1))<error]))
                   n2.adj=n2+1-max(c(1,(1:(n2+1))[cumsum(sort(pnull2))<error]))
                   id1=order(pnull1)
                   id2=order(pnull2)
                   id1=(id1[(n1+1):1])[1:n1.adj]
                   id2=(id2[(n2+1):1])[1:n2.adj]
                   pnull=pnull1[id1]%*%t(pnull2[id2])
                   if(midp==T)
                      { pvalue1[b]=sum(pnull[tnull[id1,id2]>t])+sum(pnull[tnull[id1,id2]==t])*0.5
                        pvalue2[b]=sum(pnull[tnull[id1,id2]<t])+sum(pnull[tnull[id1,id2]==t])*0.5
                      }else{
                        pvalue1[b]=sum(pnull[tnull[id1,id2]>=t])
                        pvalue2[b]=sum(pnull[tnull[id1,id2]<=t])
                      }
                  }
    
              }
            pv1=c(pv1, max(pvalue1)+(1-0.9995)); pv2=c(pv2, max(pvalue2)+(1-0.9995))
           }
          }
        
       if(rev==T)
          {pv3=pv2
           pv2=rev(pv1)
           pv1=rev(pv3)
           }
          
        return(list(pv1=pv1, pv2=pv2))
      }
