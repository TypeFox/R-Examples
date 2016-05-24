`ColorScale` <-
  function(z, loc=list(x=0, y=0), thick=1, len=1, offset=.2, col=rainbow(100),
           border='black', gradcol='black',numbcol='black',  unitscol='black',   units="", SIDE=1, font=1, fontindex =1, cex=1)
  {
##########   creates a horizontal color scale for images
######   at location loc, if loc is not provided - returns locator(2)
    if(missing(units)) { units="" }
    if(missing(SIDE)) {  SIDE=1}
    if(missing(col)) {  col= rainbow(100) }
    if(missing(font)) {  font=1  }
    if(missing(fontindex)) {  fontindex=1  }
    
    if(missing(thick))  thick=.4
    if(missing(len)) len=3
    if(missing(offset))  offset = .2
    if(missing(border))  border  = par("fg")
    ##  if(missing(gradcol))  gradcol  = 'black'
    if(missing(gradcol))  gradcol  = 'black'
    if(missing(numbcol))  numbcol  = par("fg")
    if(missing(unitscol))  unitscol  = par("fg")
    if(missing(cex))  cex=1

    
    
    
#####   loc is the boundary around the image: everything is relative to this boundary
    
    if(missing(loc)) {  loc=locator(2) }


    if(is.null(loc)) {  loc=locator(2) }

    ######   marg is the margin offset in inches
    marg = offset

    typefaces  =  c("serif","sans serif", "script", "gothic english", "serif symbol" , "sans serif symbol")
    fontindeces = c("plain", "italic", "bold", "bold italic", "cyrillic")
    typeface = typefaces[font]
    fontindex = fontindeces[fontindex]
    vfont = c(typeface, fontindex)
###  vfont=vfont


###########   SIDE =1,2,3,4  = bot, left, top, right

    u = par("usr")
    pin = par("pin")
    fin = par("fin")
    mai = par("mai")
    
    dxu = (u[2]-u[1])/fin[1]
    dyu = (u[4]-u[3])/fin[2]


    if(FALSE)
      {
        PxBOT = u[1]-mai[2]*dxu
        PyBOT = u[3]-mai[1]*dyu

        ##
        points(PxBOT, PyBOT, pch=1, col='red', cex=3, xpd=TRUE  )

        PxTOP = u[2]+mai[4]*dxu
        PyTOP = u[4]+mai[3]*dyu

        points(PxTOP, PyTOP, pch=1, col='red', cex=3, xpd=TRUE  )

        rect(PxBOT, PyBOT, PxTOP, PyTOP, border='red', xpd=TRUE)
        rect(u[1], u[3], u[2], u[4], border='green', xpd=TRUE)

      }
#########################

    rngz = range(z, na.rm = TRUE)
    rng = pretty(rngz)
    rng = rng[ rng>=rngz[1] & rng<=rngz[2]  ]
    
    
#######################    
    
    if(SIDE==1)
      {

        
        Y = c( loc$y[1]-(marg+thick)*dyu, loc$y[1]-(marg)*dyu )
        
        X = c(loc$x[1], loc$x[1]+len*dxu)


        i = seq(along = col)
        BX = (X[2]-X[1])/length(i)
        x1 =X[1]+(i-1)*BX
        x2 = x1+BX
        y1 = Y[1]
        y2 =  Y[2]
        rect(x1,y1,x2,y2,  col=col, xpd = TRUE, border=NA)

        if(!is.na(border))  rect(X[1], Y[1], X[2], Y[2], border=border, xpd=TRUE)
        
        alocs = X[1]+(X[2]-X[1])*(rng-rngz[1])/diff(rngz)
        if(!is.na(gradcol))      segments(alocs,rep(Y[2], length=length(alocs)) , alocs, rep(Y[1], length=length(alocs)),
                                          col=gradcol,  xpd=TRUE )

        ## rng[length(rng)] = paste(sep=" ",  rng[length(rng)], units)
        
        if(!is.na(numbcol))   text(alocs,rep(Y[1], length=length(alocs)), labels=rng, pos=1, vfont=vfont, col=numbcol,  cex=cex,  xpd=TRUE)
        ## text(X[2], Y[1], labels=units, adj=c( 0, 1.8 ), vfont=vfont, xpd=TRUE)
        if(!is.na(unitscol))   text(X[2], mean(Y), labels=units, pos=4, vfont=vfont, col=unitscol, cex=cex,  xpd=TRUE)
        
        
        
      }
    if(SIDE==2)
      {
        X = c( loc$x[1]-(marg*dxu+thick*dxu),   loc$x[1]-marg*dxu )
        Y = c( loc$y[2]-len*dyu, loc$y[2] )

        i = seq(along = col)
        BY = (Y[2]-Y[1])/length(i)
        y1 =Y[1]+(i-1)*BY
        y2 = y1+BY
        x1 = X[1]
        x2 =  X[2]
        rect(x1,y1,x2,y2,  col=col, xpd = TRUE, border=NA)

        alocs = Y[1]+(Y[2]-Y[1])*(rng-rngz[1])/diff(rngz)
        locsX1= rep(X[1], length=length(alocs))
        locsX2 = rep(X[2], length=length(alocs))
        if(!is.na(gradcol))    segments( locsX1 ,alocs,locsX2 , alocs, 
                                        col=gradcol,  xpd=TRUE )
        
        if(!is.na(border))  rect(X[1], Y[1], X[2], Y[2], border=border, xpd=TRUE)
        
        if(!is.na(numbcol))   text(locsX1, alocs, labels=rng, pos=2, vfont=vfont, col=numbcol,cex=cex, xpd=TRUE)
        
        if(!is.na(unitscol))   text(mean(X), Y[2], labels=units, adj=c(.5, -1 ), vfont=vfont, col=unitscol, cex=cex,xpd=TRUE)
        
        
      }
    if(SIDE==3)
      {
        Y = c(loc$y[2]+(marg)*dyu, loc$y[2]+(marg+thick)*dyu)
        
        X = c(loc$x[1], loc$x[1]+len*dxu)
        
        rect(X[1], Y[1], X[2], Y[2], border=border, xpd=TRUE)
        
        i = seq(along = col)
        BX = (X[2]-X[1])/length(i)
        x1 =X[1]+(i-1)*BX
        x2 = x1+BX
        y1 = Y[1]
        y2 =  Y[2]
        rect(x1,y1,x2,y2,  col=col, xpd = TRUE, border=NA)

        if(!is.na(border))  rect(X[1], Y[1], X[2], Y[2], border=border, xpd=TRUE)
        
        alocs = X[1]+(X[2]-X[1])*(rng-rngz[1])/diff(rngz)
        if(!is.na(gradcol))   segments(alocs,rep(Y[2], length=length(alocs)) , alocs, rep(Y[1], length=length(alocs)),
                                       col=gradcol,  xpd=TRUE )

        ## rng[length(rng)] = paste(sep=" ",  rng[length(rng)], units)
        
        if(!is.na(numbcol))  text(alocs,rep(Y[2], length=length(alocs)), labels=rng, pos=3, vfont=vfont, col=numbcol, cex=cex, xpd=TRUE)
        
        if(!is.na(unitscol))  text(X[2], mean(Y), labels=units, pos=4, vfont=vfont, col=unitscol, cex=cex, xpd=TRUE)

        
        
      }
    if(SIDE==4)
      {
        X = c(  loc$x[2]+marg*dxu ,   loc$x[2]+(marg*dxu+thick*dxu)   )
        Y = c( loc$y[2]-len*dyu, loc$y[2] )

        i = seq(along = col)
        BY = (Y[2]-Y[1])/length(i)
        y1 =Y[1]+(i-1)*BY
        y2 = y1+BY
        x1 = X[1]
        x2 =  X[2]
        rect(x1,y1,x2,y2,  col=col, xpd = TRUE, border=NA)

        alocs = Y[1]+(Y[2]-Y[1])*(rng-rngz[1])/diff(rngz)
        locsX1= rep(X[1], length=length(alocs))
        locsX2 = rep(X[2], length=length(alocs))
        if(!is.na(gradcol))    segments( locsX1 ,alocs,locsX2 , alocs, 
                                        col=gradcol,  xpd=TRUE )
        
        if(!is.na(border))  rect(X[1], Y[1], X[2], Y[2], border=border, xpd=TRUE)
        
        if(!is.na(numbcol)) text(locsX2, alocs, labels=rng, pos=4, vfont=vfont, col=numbcol, cex=cex, xpd=TRUE)
        
        if(!is.na(unitscol))    text(mean(X), Y[2], labels=units, pos=3, vfont=vfont, col=unitscol, cex=cex, xpd=TRUE)
        
        
        
        
        
        
        
      }
    

    return(list(x=X, y=Y))
    

  }

