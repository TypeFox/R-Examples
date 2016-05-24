simdiffT=function(N,a,mv,sv,ter,vp=1,max.iter=19999,eps=1e-15){
      rt=p=matrix(,N,1)

      for(jj in 1:N){
        drift=rnorm(1,mv,sv)
        p[jj]=exp(a*drift)/(1+exp(a*drift))
        M=pi*vp^2/a^2 * (exp(a*drift/(2*vp^2))+exp(-a*drift/(2*vp^2))) * 1/ (drift^2/(2*vp^2)+pi^2*vp^2 / (2*a^2))
        lmb = drift^2/(2*vp^2) +pi^2*vp^2/(2*a^2)
        ou=c()
        rej=0
       while(length(ou)<1){
         v=runif(1); u=runif(1);
         FF=pi^2*vp^4 * 1/(pi^2*vp^4+drift^2*a^2)
         sh1=1
         sh2=0
         sh3=0
         i=0
         while(abs(sh1-sh2)>eps | abs(sh2-sh3)>eps){
          sh1=sh2
          sh2=sh3
          i=i+1
          sh3= sh2 + (2*i+1)*(-1)^i*(1-u)^(FF*(2*i+1)^2)
          }
          eval=1+(1-u)^-FF * sh3
          if(v<=eval) ou=c(ou,1/lmb*abs(log(1-u)))
          else rej=rej+1
          if(rej==max.iter) stop("Rejection algorithm failed. Increase the 'max.iter' argument or try different true parameter values.")
        }
        rt[jj]=ou+ter
      }
      x=(p>matrix(runif(N)))*1
      return(list(rt=rt,x=x))
}
