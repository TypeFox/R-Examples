plotseis24<-function(JJ, dy=1/18, FIX=24, SCALE=0, FILT=list(ON=FALSE, fl=0.05 , fh=20.0, type="BP", proto="BU"), RCOLS=c(rgb(0.2, .2, 1), rgb(.2, .2, .2))  , add=FALSE )
  {
    if(missing(FIX)) { FIX=24 }
    if(missing(dy)) { dy  = 1/18 }

    if(missing(RCOLS)) { RCOLS = c(rgb(0.2, .2, 1), rgb(.2, .2, .2),"tomato3","royalblue","forestgreen","blueviolet","tan3","lightseagreen","deeppink","cyan3","bisque3","magenta1","lightsalmon3","darkcyan") }
    

    if(missing(SCALE )) { SCALE = 0 }  ###   SCALE=0 scale by trace, !=0 scale by page

    if(missing(FILT)) { FILT = list(ON=FALSE, fl=0.05 , fh=20.0, type="BP", proto="BU") }
    
    
    h = FIX
    m1 = 0
    ry =  range(c(m1, m1+23.999/24) )
    
 adt=min(JJ$gdt)
    xa = seq(from=0, length=3600/adt, by=adt)
    mx1 = min(xa)
    mx2 = max(xa)
    
    bcol = rgb(1, .8, .8)
    gcol = rgb(.8, 1, .8)
   
 rcol = RCOLS
    
    par(mar=c(5, 4, 4, 4)+0.1,  xaxs='i', yaxs='i', lwd=0.5, bty="u")

    
    if(!add)
      {
        plot( c(0, 3600), -ry  , type='n', xpd=TRUE, axes=FALSE, xlab="Time,s", ylab="")
##  90:100
    box(col=grey(0.7) )
      }

    
    tix = rep(NA, length=h)
    
  altcol = length(RCOLS)

    cols = rep(1:altcol, length=24)


    miny = rep(NA, length(24))
    maxy = rep(NA, length(24))

    
    
    for(i in 1:24 )
      {
        
        zed = JJ$sigs[[i]]
        
        fy = zed
        
        if(FILT$ON==TRUE)
          {
            
            L = length(zed)
            ipad = ceiling(L*0.02)
            
            if(ipad>10)
              {
                
###  in this case we pad each end with the time reversed
                ##      parts of the signals
                ibeg = zed[1:ipad]
                iend = zed[(L-ipad+1):L]
                
                
                ked = c(rev(ibeg), zed, rev(iend))
           ###      print("filtering")
                if(!any(is.na(ked)))
                  {
                    fy = butfilt(ked,FILT$fl, FILT$fh , adt, FILT$type , FILT$proto )
                  }
                else
                  {
                    fy = ked
                  }

                
                
                jed =  fy[(ipad+1):(ipad+L)    ]
                
                
                fy = jed
                
                
                
              }
            
          }
        fy = fy-mean(fy, na.rm=TRUE)
        miny[i] = min(fy, na.rm=TRUE)
        maxy[i]  = max(fy, na.rm=TRUE)
        JJ$sigs[[i]] = fy
            
      }
  


bigmax = max(maxy, na.rm=TRUE)
bigmin = min(miny, na.rm=TRUE)

 if(SCALE>0) rat = abs(bigmax-bigmin)/SCALE
    
    for(i in 1:24 )
      {

        a1 = m1 + (i-1)/24
        a2 = m1 + (i)/24
        y1 = -a1
        fy = JJ$sigs[[i]]
        
        zna = JJ$zna[[i]]
        icol = RCOLS[cols[i]]
        if(length(fy)>2)
          {
            if(SCALE==0)
              {
                zee  = RPMG::RESCALE(fy,  -1,   1, miny[i], maxy[i])
              }
            else
              {
                ##    w
                zee  = RPMG::RESCALE(fy, -1 ,  1, bigmin, bigmax)       
              }
            tmean = mean(zee, na.rm=TRUE)
            
            zee = dy*(zee-tmean)
            
            zee[ zna ] = NA
            
            tix[i] = y1
            lines(c(mx1, 50) , c(y1, y1)   , col=bcol )
            
            lines(xa,y1+zee, col=icol, xpd=TRUE)
          }
        else
          {
            y1 = -a1
            lines(c(mx1, 50) , c(y1, y1)   , col=bcol )
            
          }

      }

      axis(1)
 
    days = JJ$jd
    
  ##   print(days)
    modays = getmoday(days, JJ$yr[1])

    tlocs = abs(tix[!is.na(tix)])
  
    
    labs4 = format.default(1+round(24*(tlocs  - floor(tlocs))), digits=2)
    ## axis(4, at=tix[!is.na(tix)], labels=abs(days) )
    labs2 = format.default(round(24*(tlocs  - floor(tlocs))), digits=2)

    
    axis(4, at=tix[!is.na(tix)], labels=labs4, las=1)
    axis(2, at=tix[!is.na(tix)], labels=labs2, las=1) 
    ## box()

    modays = getmoday(JJ$jd[1], JJ$yr[1])

  ###  Ltit = paste(sep="/",  JJ$yr[1], modays$mo, modays$dom)

   ###    strptime(adate, "%m/%d/%y")
idate = ISOdate(JJ$yr[1], modays$mo, modays$dom, hour = 0, min = 0, sec = 0, tz = "GMT")

  adate =   format(idate,  format = "%Y-%b-%d GMT",   tz = "GMT" )
   
    mtext(adate, 3, line=1, at=0, adj=0)

if(FILT$ON==TRUE)
  {

    fl = FILT$fl
    unitlow = "Hz"
    fh = FILT$fh
    unithi = "Hz"
    ftype= FILT$type

    if(fl<1)
      {
        fl = 1/fl
        unitlow = "s"

      }

    if(fh<1)
      {
        fh = 1/fh
        unithi = "s"

      }


    
    filttag =  paste(sep= " ", ftype,fl, unitlow, fh, unithi  )
    
    mtext(filttag, 3, line=1, at=mx2, adj=1)
    
  }


    if(SCALE>0)
      {

           mtext("Scaled by window", 1, line=3, at=0, adj=0)

      }


    
    
    invisible(list( x=xa, y=tix,  yr=JJ$yr[1],  jd=JJ$jd[1] ))
  }
