#source('poly.R')

fm_unred = function( wave, flux, ebv,
              r_v = 3.1,
              avglmc=FALSE, lmc2=FALSE,
              x0=NULL, gamma=NULL, c4=NULL, c3=NULL, c2=NULL, c1=NULL) {

 x = 10000./ wave                # convert to inverse microns
 curve = x*0.
if(!missing(lmc2) ) {
        if(missing(x0)) x0    =  4.626
        if(missing(gamma) )gamma =  1.05
        if(missing(c4) )c4   =  0.42
        if(missing(c3) )c3    =  1.92
        if(missing(c2) )c2    = 1.31
        if(missing(c1) )c1    =  -2.16
 }
 else if(!missing(avglmc) ){
        if(missing(x0) )x0 = 4.596
        if(missing(gamma) )gamma = 0.91
        if(missing(c4) )c4   =  0.64
        if(missing(c3) )c3    =  2.73
        if(missing(c2) )c2    = 1.11
        if(missing(c1) )c1    =  -1.28
  }
  else {
        if(missing(x0) )x0    =  4.596
        if(missing(gamma) )gamma =  0.99
        if(missing(c3) )c3    =  3.23
        if(missing(c4) )c4   =  0.41
        if(missing(c2) )c2    = -0.824 + 4.717/r_v
        if(missing(c1) )c1    =  2.030 - 3.007*c2
 }

 xcutuv = 10000.0/2700.0
 xspluv = 10000.0/c(2700.0,2600.0)
 iuv = (x>=xcutuv);
 iopir = !iuv

 xuv = c(xspluv,x[iuv])
 yuv = c1  + c2*xuv
 yuv = yuv + c3*xuv^2/((xuv^2-x0^2)^2 +(xuv*gamma)^2)
 yuv = yuv + c4*(0.5392*((xuv>5.9)-5.9)^2+0.05644*((xuv>5.9)-5.9)^3)
 yuv = yuv + r_v
 yspluv  = yuv[1:2]                  # save spline points
 curve[iuv] = yuv[-(1:2)]      # remove spline points

 xsplopir = c(0,10000.0/c(26500.0,12200.0,6000.0,5470.0,4670.0,4110.0))
 ysplir   = c(0.0,0.26469,0.82925)*r_v/3.1
 ysplop   = c(polyidl(r_v, c(-4.22809e-01, 1.00270, 2.13572e-04) ),
             polyidl(r_v, c(-5.13540e-02, 1.00216, -7.35778e-05) ),
             polyidl(r_v, c( 7.00127e-01, 1.00184, -3.32598e-05) ),
             polyidl(r_v, c( 1.19456, 1.01707, -5.46959e-03, 7.97809e-04,
                     -4.45636e-05) ) )

 ysplopir = c(ysplir,ysplop)

 #Arnab: Replaced a call to cspline by spline.
 curve[iopir] = spline(c(xsplopir,xspluv),c(ysplopir,yspluv),xout=x[iopir])$y
 # now apply extinction correction to input flux vector
 curve = ebv*curve
 funred = flux * 10.^(0.4*curve)       #derive unreddened flux

  if(avglmc)
     return(list(funred = funred, gamma=gamma, xo=x0,
              c1=c1, c2=c2, c3=c3, c4=c4,extcurve = curve - r_v))
  else
     return(list(funred = funred, extcurve = curve - r_v))
}
