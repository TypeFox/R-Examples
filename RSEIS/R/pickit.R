`pickit` <-
function(ay, deltat=0.008 ,  MED=225, FRWD=8,  BKWD=8,  sbef=1, saft=6, thresh=2, Tthresh2=7, stretch=1000, flo=0.1, fhi=5.0, Kmin=7, dthresh=.01, threshbot=1.01)
  {
    if(missing(deltat)) { deltat=0.008 }
    if(missing(FRWD)) { FRWD=8 }
    if(missing(BKWD)) { BKWD=8 }
    
    if(missing(sbef)) { sbef = 1 }
    if(missing(saft)) { saft = 6 }
    if(missing(thresh)) { thresh=2 }
    if(missing(dthresh)) { dthresh=.01 }
    if(missing(threshbot)) { threshbot=1.01 }

    
    if(missing(Tthresh2)) { Tthresh2= (sbef+saft) }
    if(missing(stretch)) { stretch=1000 }
    if(missing(flo)) { flo = .1 }
    if(missing(fhi)) { fhi=5.0 }
    if(missing(Kmin)) { Kmin = 7 }
    if(missing(MED)) { MED = 1001 }

    if(fhi > 1/(2*deltat)) { fhi = (1/(2*deltat))- 0.05*(1/(2*deltat)) }

    
    ##########   bandpass filter the data, if both flo and fhi are negative, do not filter
    if(flo>0 & fhi>0)
      {
        fy = butfilt(ay , flo , fhi , deltat , "BP" , "BU" )
      }
    else
      {
        fy = ay
      }
    
    LEN1 = FRWD/deltat
    LEN2 = BKWD/deltat

    
    A = STLTcurve(fy, dt=deltat, fwlen = LEN1,  bwlen  = LEN2, stretch=stretch, MED=MED, PLOT=FALSE)

    
###        A = STLTcurve(fy, dt=deltat, fwlen = LEN1,  bwlen  = LEN2, PLOT=TRUE)


    x = 1:length(A$rat)

       argh = jstats(A$rat)
     stthresh = argh$bstats[5]
  
  ###   if(thresh<2*stthresh) { thresh = 2*stthresh }
    if(thresh<stthresh) { thresh = stthresh }
  
    J = Thresh.J(A$rat,thresh)

    if(is.null(J))
      {

        print("Nothing in Thresh.J")
        
        return(NULL)


      }
    
    Kthresh2 = Tthresh2/deltat
    
    Z = J$J[(J$L-J$J)>Tthresh2]

    ##  lower thresh hold until have a minimum number of picks and are above 1.01 S/N
    ##  thresh drops below 1.01 then we are picking noise.
    ##  the min number of picks (7) should be a parameter to adjust to dataset

  if(FALSE)
    {
    while(length(Z)<Kmin & thresh>threshbot )
      {
        thresh = thresh-dthresh
        J = Thresh.J(A$rat,thresh)
        
        Z = J$J[(J$L-J$J)>Kthresh2]
      }
  }

    
    a1 = x[J$J]-sbef/deltat
    a2 = x[J$L]+saft/deltat

###  get rid of overlapping traces
    tt4 = rep(0, length(x))
    for(j in 1:length(a1))
      {
        tt4[x>=a1[j]&x<=a2[j]] = 1
      }
    
    J = Thresh.J(tt4,.5)
    a1 = x[J$J]
    a2 = x[J$L]

    return(list(RAT=A$rat, x=x, ay=ay, fy=fy, deltat=deltat, J=J$J , Z=Z, a1=a1, a2=a2, thresh=thresh, Tthresh2=Tthresh2, Kmin=Kmin) )

    
  }

