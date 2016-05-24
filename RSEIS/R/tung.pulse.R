`tung.pulse` <-
function(r,q, dt)
  {
    ###  plot a small pulse and calculate a characterization
    ### 
    ## chug.pulse(r,q,dt)
    
    jp = hilow(q)
    
    t1 = min(r)
    mid = mean(r)
    dmid = (mid-r[jp$lo])
    
    
    omid = c(max(which(sign(dmid)==1)), min(which(sign(dmid)==-1)))
    mins= sort(jp$lo[c(omid[1], omid[2])])
    
    plot(r,q, type='b')
    title("CLICK TWICE IN WINDOW: either side of pulse: end with middle mouse")
    ##  title(main=paste(sep=' ', i,'of',length(AL))) 
    points(r[jp$hi],q[jp$hi], col=2)
    points(r[jp$lo],q[jp$lo], col=3)
    
    abline(v=r[mins], col=2)
    
    vp = plocator(COL=rgb(1.0,0.5,0.5))
    
    nvp = length(vp$x)

    rvp = sort(vp$x)
    
    if(nvp>1)
      {
        ###   source("/home/lees/Progs/R_stuff/tung.R");
     left =     which.min(abs(r-rvp[nvp-1]) )
     right =    which.min(abs(r-rvp[nvp]) )
    
      }
    mins = c(left, right)
    
    abline(v=r[c(left, right)], col=3)
    
    s = q[r>=r[left]&r<=r[right]]
        
    sum0 = sqrt(sum(s*s)/length(s))
    
    imax = left-1+which.max(q[r>=r[left]&r<=r[right]]^2)

    
    abline(v=r[imax])
   
    ### p1 = q[jp$hi[imax]]-q[jp$lo[c(omid[1], omid[2])]]
    
    ###  tees = t1+r[jp$lo[c(omid[1], omid[2])]]
    ###   t0 = mid+t1

    Ex = r[mins]
    Ey = q[mins]

    ###  find the max point that lies between the two low points


    Cx = r[imax]
    Cy = q[imax]

    
   ###   Cx = r[jp$hi[imax]]
    ###  Cy = q[jp$hi[imax]]

    lines(c(Ex[1],  Ex[2], Cx, Ex[1]), c(Ey[1],  Ey[2], Cy, Ey[1]))

    A1 = c(Ex[1]-Cx, Ey[1]-Cy, 0)
    A2 = c(Ex[2]-Cx, Ey[2]-Cy, 0)
    
    ar2 = 0.5*vlen(xprod(A1,A2))

    tr = r[r>=Ex[1]&r<=Ex[2] ]
    tq = q[r>=Ex[1]&r<=Ex[2] ]

 
     DefInt = integ1(tr, tq)

    ### print(paste(sep=' ',DefInt[1],DefInt[2] , s4, s3))
    #### Ex[1], Ex[2] = left minimum
     #### Ey[1], Ey[2] = right  minimum
     ####    Cx, Cy  = center (max?)
    ####    7:  ar2 = area of triangle
    ####    8:  DefInt[1]  = integral under curve
    ####    9:  DefInt[2]  = integral under curve ( bottom triangle removed)
    ####   10:  sum0   = RMS amplitude
    
    return(c(Ex[1], Ex[2], Ey[1], Ey[2], Cx, Cy, ar2, DefInt[1], DefInt[2], sum0))
  }

