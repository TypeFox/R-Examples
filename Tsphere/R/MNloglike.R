MNloglike <-
function(x,M,Sig,Delt,rhor,rhoc,qr=2,qc=2,Sigi=NULL,Delti=NULL)
{
  n = nrow(x)
  p = ncol(x)
  if(length(Sigi)==0){  Sigi = solve(Sig) }
  if(length(Delti)==0){   Delti = solve(Delt) }
  if(qr==2){ tr = sum((rhor*Sigi)^2)} else { tr = sum(abs(rhor*Sigi))}
  if(qc==2){ tc = sum((rhoc*Delti)^2)} else { tc = sum(abs(rhoc*Delti))}
  if(det(Sigi)==0){ t1 = log(1e-300)
                  }else{ t1 = log(det(Sigi))}
  if(sum(is.na(x))==0)
    {
      val = (p/2)*t1 + (n/2)*log(det(Delti)) - (1/2)*sum(diag(Sigi%*%(x-M)%*%Delti%*%t(x-M))) - tr - tc
    }
  else
    {
      xs = x
      xc = x - M
      tlds = 0
      for(j in 1:p)
        {
          oj = !is.na(x[,j])
          sigioj = solve(Sig[oj,oj])
          tlds = tlds + log(det(sigioj))
          xs[oj,j] = chol(sigioj)%*%xc[oj,j]
        }
      tldd = trt = 0
      for(i in 1:n)
        {
          oi = !is.na(x[i,])
          deltioi = solve(Delt[oi,oi])
          tldd = tldd + log(det(deltioi))
          trt = trt + sum(diag(xs[i,oi]%*%t(xs[i,oi])%*%deltioi))
        }
      val = (1/2)*(tlds + tldd - trt) - tr - tc
    }
  return(val)
}

