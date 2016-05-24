ci.RiskD <-
function(x1, x2, n1, n2, cov.prob=0.95, BB.grdnum=2000, midp=T)
          {if(x1*x2==0)
            {p1=(x1+0.5)/(n1+1) 
             p2=(x2+0.5)/(n2+1)
             delta=p2-p1
             sigma=sqrt((1-p1)*p1/(n1+1)+(1-p2)*p2/(n2+1))
             }
          if(x1*x2>0)
            {p1=x1/n1 
             p2=x2/n2
             delta=p2-p1
             sigma=sqrt((1-p1)*p1/n1+(1-p2)*p2/n2)
             }

          lower=min(max(delta-5*sigma, -1), -1/n1-1/n2)
          upper=max(min(delta+5*sigma,  1),  1/n1+1/n2)
          delta.grd=sort(unique(c(0, seq(lower, upper, length=BB.grdnum))))
          BB.grdnum=length(delta.grd) 

          fit=priskD.exact(x1, x2, n1, n2, delta.grd, midp=midp)
          pv1=fit$pv1
          pv2=fit$pv2
          for(j in 1:BB.grdnum)
             pv1[BB.grdnum-j+1]=max(pv1[1:(BB.grdnum-j+1)]); pv2[j]=max(pv2[j:BB.grdnum])
          

          lower=range(delta.grd[pv1>(1-cov.prob)/2 & pv2>(1-cov.prob)/2])[1]
          upper=range(delta.grd[pv1>(1-cov.prob)/2 & pv2>(1-cov.prob)/2])[2]
          est=delta.grd[abs(pv2-pv1)==min(abs(pv2-pv1))][1]
          
          id0=(1:BB.grdnum)[delta.grd==0]
          p=min(1, 2*min(pv1[id0], pv2[id0]))


          status=1
          if(min(delta.grd)==lower || max(delta.grd)==upper)
             status=0

          return(list(est=est, lower=lower, upper=upper, status=status, p=p))
          }
