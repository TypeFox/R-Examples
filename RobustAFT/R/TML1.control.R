"TML1.control" <-
function(tlo=0.0001, mxf=50, mxs=50, ntm=50, tls=1e-6, k0=1.5477, h=100) {              
 old <- comval(); dfcomn2(ipsi=4, xk=1.5477)
 Beta0 <- integrate(Chiphi, -10, 10)$value
 dfcomn2(ipsi=old$ipsi, xk=old$xk)
 list(tlo=tlo,mxf=mxf,mxs=mxs,ntm=ntm,tls=tls,k0=k0,h=h,Beta0=Beta0)}

