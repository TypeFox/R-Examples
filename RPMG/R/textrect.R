textrect<-function(x, y, lab, textcol="black" ,  col="white", border='black' , off=0.06,  brd=.06 , pos=1, log="" , add=TRUE, ...)
  {
    
    if(missing(col)) col='white'
    if(missing(border)) border='black'
    if(missing(textcol)) textcol='black'
    
    if(missing(brd)) brd = 0.06
    if(missing(off)) off = 0.06
    
    if(missing(pos)) pos=0
    if(missing(add)) add= TRUE
    if(missing(log)) log = ""
    
    ####  to vectorize this, x,y,lab must all have the same length

    n = length(x)
    
    if(length(pos)<n)
      {
       pos =  rep(pos, length=n)
      }

    if(length(col)<n)
      {
       col =  rep(col, length=n)
      }

      if(length(border)<n)
      {
       border =  rep(border, length=n)
      }
       if(length(textcol)<n)
      {
       textcol =  rep(textcol, length=n)
      }
    


    
    u = par("usr")
 
     if(length(grep("x", log))>0 )
       {
         u[1] = 10^u[1]
         u[2] = 10^u[2]
       }
     if(length(grep("y", log))>0 )
       {
         u[3] = 10^u[3]
         u[4] = 10^u[4]
       }
       
    
    f = par("fin")
      
    p = par("pin")

    inchoff = 0.06
    inchdx = (u[2]-u[1])/p[1]
    inchdy = (u[4]-u[3])/p[2]
    
     
    thepos = c(0, seq(from=1, to=4.5, by=.5))

    offang = pi/2 -  seq(from=-pi, to=pi-(pi/4), by=pi/4)
    
    offang*180/pi
    

    cx = cos(offang)
    cy = sin(offang)

  ###  labelpos    0    1    1.5  2   2.5   3   3.5   4  4.5

    offx = off*inchdx
    offy = off*inchdy

    brdx = brd*inchdx
    brdy = brd*inchdy
    
   SHIFTMAT = matrix(
      c(-0.5, -1, -0.5, -1,
        -0.5, -1, -1, -2,
        -1, -2, -1, -2,
        -1, -2, -.5, -1,
        -1, -2, 0, 0,
        -.5, -1, 0, 0,
        0, 0, 0, 0,
        0, 0,-0.5 ,-1, 
        0, 0,-1 , -2), byrow=TRUE, ncol=4)


    rectmatrix = matrix(nrow=n, ncol=4)

    for(J in 1:n)
      {
       
        
        W = strwidth(lab[J], units = "user", ...)
        H = strheight(lab[J], units = "user", ...)


     if(length(grep("x", log))>0 )
       {
        W = 10^W
       }
     if(length(grep("y", log))>0 )
       {
         H =  10^H
       }
       




        
        kpos = which(pos[J]== thepos)

        ipos = thepos[kpos]

        
        if(ipos==0)
          {
            ax =0
            ay =0

            
          }
        else
          {
            j = kpos-1  
            ax = offx*cx[j]
            ay = offy*cy[j]
            
            
          }

        
        
        p1x = x[J]+ax+SHIFTMAT[kpos, 1]*W+SHIFTMAT[kpos, 2]*brdx
        p1y = y[J]+ay+SHIFTMAT[kpos, 3]*H+SHIFTMAT[kpos, 4]*brdy

        p2x =p1x +W+2*brdx
        p2y = p1y+H+2*brdy

        rectmatrix[J, ] = c(p1x, p1y, p2x, p2y)
        

        ##     title(sub=paste(i, theadj[i,1] , theadj[i,2], sep=" ")) 

        
      ##  a1x = (p2x+p1x)/2
      ##  a1y = (p1y+p2y)/2
        
       ## rect(p1x, p1y, p2x, p2y , col=col[J], border=border[J], xpd=TRUE)
     ##   text(a1x, a1y ,labels=lab[J],  col=textcol[J], adj=c(.5, .5), xpd=TRUE, ... )
        ##   points(x+ax,y+ay, pch=3, col="magenta") 
      }


    if(add==TRUE)
      {

        for(J in 1:n)
          {


            p1x =rectmatrix[J, 1]
            p1y = rectmatrix[J, 2]

            p2x = rectmatrix[J, 3]
            p2y = rectmatrix[J, 4]

            
            a1x = (p2x+p1x)/2
            a1y = (p1y+p2y)/2
            
            rect(p1x, p1y, p2x, p2y , col=col[J], border=border[J], xpd=TRUE)
            text(a1x, a1y ,labels=lab[J],  col=textcol[J], adj=c(.5, .5), xpd=TRUE, ... )
          }


      }


    
    invisible(rectmatrix)


    
  }


