DO.DYKE <-
function(a=a, x=x, t=t, k=k, T0=T0, NDIM=TRUE)
  {
    if(missing(a)) { a = 6 }
    if(missing(x)) { x = seq(0,20, length=1000) }
    if(missing(T0)) { T0 = 1000 }
    if(missing(k)) { k = 1e-6  }
    if(missing(t)) { t = seq(from=0,to=1000*24*3600, by=100*24*3600)  }
       if(missing(NDIM)) { NDIM = TRUE }

    TMAT = matrix(ncol=length(t), nrow=length(x))
    
    AX = a-x
    flag1 = AX<0
    
   
    AX = a+x
    flag2 = AX<0

    i = 0
    for(tim in t)
      {
        i = i+1
        tim = t[i]
        e1 = erf(abs(a-x)/(2*sqrt(k*tim)))
        e1[flag1] = -e1[flag1]
        
        
        
        e2 = erf(abs(a+x)/(2*sqrt(k*tim)))
        e2[flag2] = -e2[flag2]
    
        
        Temp  =  0.5*(e1 + e2)
        TMAT[,i] = Temp 
      }

    
   ###  matplot(x/a,TMAT, type='l')

    
    BMAT =   rbind( TMAT[ rev(1:length(x)), ], TMAT)

    EX = c( rev(-x/a), x/a)
   
    if(NDIM==TRUE)
      {
        matplot(EX , BMAT , type='l', xlab="Distance, x/a", ylab="Temperature, T/T0" )
      }
    else
      {
        matplot(a*EX , T0*BMAT , type='l', xlab="Distance, m", ylab="Temperature, K" )
      }
    grid(col='darkgray')

    px = rep(0,length(t))
    py = rep(0,length(t))
    
    tx = a*seq(from=0.01, to=.99, length=length(t))
    i =0
    for(tim in t)
      {
        i = i+1
        tim = t[i]
        e1 = erf(abs(a-tx[i])/(2*sqrt(k*tim)))
        e2 = erf(abs(a+tx[i])/(2*sqrt(k*tim)))
         py[i] =  0.5*(e1 + e2)
      }
    
   #   labs = k*t/a^2
    if(NDIM==TRUE)
      {
        labs = format.default(k*t/a^2, digits=3)
        text(tx/a, py, labels=labs)
      }
    else
      {
        tlabs = t/(24*3600)
        labs = format.default(tlabs, digits=3)
        text(tx, T0*py, labels=labs)

      }
  }

