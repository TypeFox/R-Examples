ci.RateD <-
function(x1, x2, e1, e2, cov.prob=0.95, BB.grdnum=2000, midp=T)
          {if(x1*x2==0)
            {lambda1=(x1+0.5)/(e1+1) 
             lambda2=(x2+0.5)/(e2+1)
             delta=lambda2-lambda1
             sigma=sqrt(lambda1/(e1+1)+lambda2/(e2+1))
             }
          if(x1*x2>0)
            {lambda1=x1/e1 
             lambda2=x2/e2
             delta=lambda2-lambda1
             sigma=sqrt(lambda1/e1+lambda2/e2)
             }

          lower=min(delta-5*sigma, -1/e1-1/e2)
          upper=max(delta+5*sigma,  1/e1+1/e2)
          delta.grd=sort(unique(c(0, seq(lower, upper, length=BB.grdnum))))
          BB.grdnum=length(delta.grd) 

          fit=prateD.exact(x1, x2, e1, e2, delta.grd, midp=midp)
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
