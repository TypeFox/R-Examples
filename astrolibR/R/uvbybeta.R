#source('deredd.R')
#source('poly.R')

uvbybeta = function(xby,xm1,xc1,xhbeta,xn, eby_in, name) {
	
  rm1 = -0.33 
  rc1 = 0.19
  rub = 1.53         # Parameter values
  
  nstar  = length(xby)
  xub = xc1 + 2*(xm1+xby)
  xflag1 = ifelse(xhbeta==0,1,0)
  if(missing(name)) name = 'star'

  do_eby = missing(eby_in)
  te = numeric(nstar)
  mv = te
  delm0 = te
  radius = te

  hbeta.out = 
    hbeta.status =
      delm0.status =
        by0.out=
          m0.out = 
            c0.out = te
  
  warn.out = character(nstar)
  
  if(do_eby)
    eby = te
  else
    eby = rep(eby_in,nstar)

  if(any(!(xn%in%(1:8)))) stop("Group must be in 1:8")
  
  for(i in 1:nstar) {
    by = xby[i]
    m1 = xm1[i]
    c1 = xc1[i]
    hbeta = xhbeta[i]
    n = as.integer(xn[i])
    ub = xub[i]
    flag1 = xflag1[i] 
    flag2 = 0
    warn = ''

    if(n==1) {
      if(do_eby) eby[i] = ( 13.608*by-ub+1.467 ) / (13.608-rub)
      tmp=deredd(eby[i], by, m1, c1, ub)
      by0=tmp$by0
      m0 = tmp$m0
      c0=tmp$c0
      ub0=tmp$ub0
      if (flag1==1) hbeta = 
        polyidl(c0, c(2.61033, 0.132557, 0.161463, -0.027352) )
      g = log10(hbeta - 2.515) - 1.6*log10(c0 +0.322)
      mv[i] = 3.4994 + 7.2026*log10(hbeta - 2.515) -2.3192*g + 2.9375*g^3
      te[i] = 5040/(0.2917*c0 + 0.2)  
      m0zams = polyidl(c0, c(0.07473, 0.109804, -0.139003, 0.0957758) )
      delm0[i] = m0zams - m0
      flag2 = 1
    }
    else if(n==2) {
      if(do_eby ){
        eub = ( 1.5*c1 - ub + 0.035) / (1.5/(rub/rc1)-1)
        eby[i] = eub/rub
      }
      tmp=deredd(eby[i], by, m1, c1, ub)
      by0=tmp$by0
      m0 = tmp$m0
      c0=tmp$c0
      ub0=tmp$ub0
      if(( flag1==1 ) )hbeta = 0.037*c0 + 2.542
    }
    else if(n==3) {
      if(do_eby ){
        eub = (1.36*c1-ub+0.004) / (1.36/(rub/rc1)-1)
        eby[i] = eub/rub
      }
      tmp=deredd(eby[i], by, m1, c1, ub)
      by0=tmp$by0
      m0 = tmp$m0
      c0=tmp$c0
      ub0=tmp$ub0
      if(flag1==1 )hbeta = 0.047*c0 +2.578
    }
    else if(n==4) {
      if(do_eby ){
        eub = ( 1.32*c1 - ub - 0.056) / ( 1.32 / (rub/rc1)-1 )
        eby[i] = eub/rub
      }
      tmp=deredd(eby[i], by, m1, c1, ub)
      by0=tmp$by0
      m0 = tmp$m0
      c0=tmp$c0
      ub0=tmp$ub0
      if(( flag1==1 ) )hbeta = 0.066*c0+2.59
    }
    else if(n==5) {
      if(do_eby ){
        m = m1 - rm1*by
        by0 = 4.2608*m^2 - 0.53921*m - 0.0235
        while(TRUE) {
          bycorr = by0
          m0 = m1 - rm1*(by-bycorr)
          by0 = 14.0881*m0^2 - 3.36225*m0 + 0.175709
          if( abs(bycorr - by0)<0.001) break
        }
        eby[i] = by - by0
      }
      tmp=deredd(eby[i], by, m1, c1, ub)
      by0=tmp$by0
      m0 = tmp$m0
      c0=tmp$c0
      ub0=tmp$ub0
      if(flag1==1 )hbeta = 2.7905 - 0.6105*by + 0.5*m0 + 0.0355*c0
      r = 0.35*(c1-rc1*by) - (hbeta-2.565)
      a0 = by0+ 0.18*(ub0-1.36)
      mv[i] = 1.5 + 6.0*a0 - 17.0*r
      te[i] =  5040. /(0.7536 *a0 +0.5282)
      m0zams = -3.95105*by0^2 + 0.86888*by0 + 0.1598
      delm0[i] = m0zams - m0
    }

    else if(n==6) {
      if(flag1 ){
        warn = ' estimate of hbeta only validif(star is unreddened'
        hbeta = 3.06 - 1.221*by - 0.104*c1
      }
      m1zams = -2.158*hbeta^2 +12.26*hbeta-17.209
      if(( hbeta<=2.74 ) ){
        c1zams = 3.0*hbeta - 7.56
        mvzams = 22.14 - 7*hbeta
      } else if(( ( hbeta>2.74 ) && ( hbeta<=2.82 ) ) ){
        c1zams = 2.0*hbeta - 4.82
        mvzams = 11.16-3*hbeta
      } else {
        c1zams = 2.0*hbeta-4.83
        mvzams =-88.4*hbeta^2+497.2*hbeta-696.41
      }        
      if(do_eby ){
        delm1 = m1zams - m1
        delc1 = c1-c1zams
        if(delm1<0. )
          by0 = 2.946 - hbeta - 0.1*delc1 - 0.25*delm1 else 
        by0 = 2.946 - hbeta - 0.1*delc1
        eby[i] = by - by0
      } 
      tmp=deredd(eby[i], by, m1, c1, ub)
      by0=tmp$by0
      m0 = tmp$m0
      c0=tmp$c0
      ub0=tmp$ub0
      delm0[i] = m1zams - m0
      delc0 = c0 - c1zams
      mv[i] = mvzams -9.0*delc0
      te[i] = 5040 / (0.771453*by0 + 0.546544)
    }
    else if(n==7) {
      if(flag1 ){ 
        byinit = by
        m1init = m1
        for(ii  in ( 1):(1000)) {
          m1by = 2.5*byinit^2 - 1.32*byinit + 0.345
          bycorr = byinit + (m1by-m1init) / 2.0
          if(( abs(bycorr-byinit)<=0.0001 ) ) break
          byinit = bycorr
          m1init = m1by
        }
        hbeta = 1.01425*bycorr^2 - 1.32861*bycorr + 2.96618 
      }
      m1zams = polyidl(hbeta, c( 46.4167, -34.4538, 6.41701) )
      mvzams = polyidl(hbeta, c(324.482, -188.748, 11.0494, 5.48012))
      if(hbeta<=2.65 )
        c1zams = 2*hbeta - 4.91 else 
      c1zams = 11.1555*hbeta^2-56.9164*hbeta+72.879
      if(do_eby ){
        delm1 = m1zams - m1
        delc1 = c1 - c1zams
        dbeta = 2.72 - hbeta
        by0 = 0.222+1.11*dbeta +2.7*dbeta^2-0.05*delc1-(0.1+3.6*dbeta)*delm1
        eby[i] = by - by0
      }
      tmp=deredd(eby[i], by, m1, c1, ub)
      by0=tmp$by0
      m0 = tmp$m0
      c0=tmp$c0
      ub0=tmp$ub0
      delm0[i] = m1zams - m0
      delc0 = c0 - c1zams
      f = 9.0 + 20.0*dbeta
      mv[i] = mvzams - f*delc0
      te[i] = 5040/(0.771453*by0 + 0.546544)
    }
    else if(n==8) {
      if(( flag1==1 ) )flag1 = 2
      if(( by<=0.65 ) )
        eby[i] = (5.8651*by - ub -0.8975) / (5.8651 - rub) 
      else if(( ( by>0.65 ) && ( by<0.79 ) ) ){ 
        eby[i] = (-0.7875*by - c1 +0.6585)/(-0.7875 - rc1)
        by0 = by - eby[i]
        if(( by0<0.65 ) )
          eby[i] = (5.8651*by - ub -0.8975) / (5.8651-rub)
      } else { 
        eby[i] = ( 0.5126*by - c1 - 0.3645 ) / (0.5126-rc1)
        by0 = by - eby[i]
        if(( by0<0.79 ) ) 
          eby[i] = (-0.7875*by - c1 + 0.6585) / (-0.7875-rc1)
        by0  = by - eby[i]
        if(( by0<0.65 ) ) 
          eby[i] = ( 5.8651*by - ub - 0.8975) / (5.8651-rub)
      } 
      
      tmp=deredd(eby[i], by, m1, c1, ub)
      by0=tmp$by0
      m0 = tmp$m0
      c0=tmp$c0
      ub0=tmp$ub0
      m1zams = polyidl( by0, c(7.18436, -49.43695, 122.1875, -122.466, 42.93678)) 
      if (by0<0.65) {
        c1zams = polyidl(by0, c(3.78514, -21.278, 42.7486, -28.7056 ) )
        mvzams =  
          polyidl(by0, c(-59.2095, 432.156, -1101.257, 1272.503, -552.48))
      }
      else if( (by0>=0.65) && (by0<0.79)) {
        c1zams = -0.631821*by0^2+0.116031*by0+0.33657
        mvzams = 1.37632*by0^2 + 4.97911*by0+3.4305
      }
      else {
        c1zams = -0.010028*by0^2 + 0.530426*by0 - 0.37237
        mvzams =  1.18298*by0^2  + 3.92776*by0 + 4.37507
      }
      delm0[i] = m1zams - m0
      delc0 =c0 - c1zams
      if (by0<=0.505) {
        f = 10. - 80.*(by0-0.38)
        te[i] = 10^(-0.416*by0+3.924)
      }
      else {
        f = 0.0
        te[i] = 10^(-0.341*by0+3.869)
      }
      mv[i] = mvzams - f*delc0 + 3.2*delm0[i] - 0.07
    }
    
    
    if((n>=2) && ( n<=4 ) ){
      betaza = polyidl(c0, c(2.62745, 0.228638, -0.099623, 0.277363, -0.160402 ) )
      b = betaza - 2.5
      mvzams =203.704*b^3 - 206.98*b^2 + 77.18*b - 9.563
      dbeta = betaza - hbeta
      dmv = -121.6*dbeta^2 +61.0*dbeta + 0.08
      mv[i] = mvzams - dmv
      te[i] = 5040 / (0.35866*ub0 + 0.27346)
      flag2 = 2
    }
    if( by0<=0.335  )
      fv = -6.759*by0^3 + 3.731*by0^2 - 1.092*by0 + 3.981 
    else
      fv = -0.534*by0 + 3.959

    radius[i] = 10^(2.*(4.236-0.1*mv[i] - fv))

                                        #Arnab: Collect the scalar computed quantities
    hbeta.out[i] = round(hbeta*1000)/1000
    hbeta.status[i] = flag1
    delm0.status[i] = flag2
    by0.out[i] = by0
    m0.out[i] = m0
    c0.out[i] = c0
    warn.out[i] = warn
  }
  
  
  teff = round(te/10)*10

  res = data.frame(
    name=name,
    group=n,
    by= xby,
    m1=xm1,
    c1=xc1,
    hbeta = hbeta.out,
    hbeta.status = hbeta.status,
    by0 = by0.out,
    m0 = m0.out,
    c0= c0.out,
    eby = eby,
    mv = mv,
    radius= radius,
    delm0.status=delm0.status,
    delm0=delm0,
    teff=teff,
    warn=warn.out,stringsAsFactors=FALSE)
  
  return(res)
}
