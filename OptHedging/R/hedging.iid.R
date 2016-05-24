hedging.iid = function(R,T,K,r,put,n,m,minS,maxS)
{

 N=length(R) 


out0 = .C("HedgingIID",
                  as.double(R),
                  as.double(T),
                  as.double(K),
                  as.double(r),
                  as.integer(put),
                  as.integer(n),
                  as.integer(m),
                  as.double(maxS),
                  as.double(minS),
                  as.integer(N),
                  S         = double(m),
                  Cvec      = double(n*m),
                  avec      = double(n*m),
                  rho       = double(1)
                  )


 S        = out0$S;
 rho      = out0$rho;

 C        = matrix(out0$Cvec,nrow=n,ncol=m);
 a        = matrix(out0$avec,nrow=n,ncol=m);
 
 phi1 = (a[1,] - C[1,]*rho)/S;

out.names = c('S', 'C','a', 'rho','phi1')
  out = list(S=S, C=C,a=a,rho=rho,phi1=phi1)
out
 
 
}
