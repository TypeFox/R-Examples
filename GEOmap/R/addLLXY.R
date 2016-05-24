`addLLXY` <-
function(lats, lons, PROJ=PROJ, PMAT=NULL, col=gray(0.7), GRID=TRUE, GRIDcol=1, LABS=NULL, LABcol=1, BORDER=NULL, TICS=c(1,1), xpd=TRUE )
  {

    if(missing(col)) {col=gray(0.7)}
    if(missing(GRID)) { GRID = TRUE }
    if(missing(GRIDcol)) { GRIDcol = gray(0.7) }
    if(missing(LABS)) { LABS =NULL }
    if(missing(LABcol)) { LABcol = 1 }
    
    if(missing(BORDER)) { BORDER=NULL }
    if(missing(TICS)) { TICS=NULL }
    if(missing(xpd)) { xpd=TRUE }

    if(LABS==FALSE) { LABS = NULL }
    
  ############   Hershey Fonts ( vfont  typeface  fontindex )  degrees = de   minutes = fm   seconds = sd

   typefaces  =  c("serif","sans serif", "script", "gothic english", "serif symbol" , "sans serif symbol")

    fontindeces = c("plain", "italic", "bold", "bold italic", "cyrillic")

    typeface = typefaces[1]
    fontindex = fontindeces[1]

    vfont = c(typeface, fontindex)
    
    HL = 0
    tpoints = as.list(NA)


    if(GRID==FALSE) { GRIDcol = NULL }
    
########  latitudes

    if(!is.null(GRIDcol) )
      {
        col = GRIDcol
        for(i in lons)
          {
            xy = GLOB.XY(lats ,  rep(RPMG::fmod(i, 360), length(lats)) ,  PROJ)
            if(!is.null(PMAT))
              {
                tem = trans3d(xy$x, xy$y, rep(0, length(xy$y)) , PMAT)
                xy$x = tem$x
                xy$y = tem$y
              }

            lines(xy$x, xy$y, col=GRIDcol, lty=2, lwd=.4, xpd=xpd )
            HL = HL +1
            tpoints[[HL]] = xy
                    
          }

########  longitudes
        for(i in lats)
          {
            xy = GLOB.XY(rep(i, length(lons)), lons, PROJ  )

            if(!is.null(PMAT))
              {
                tem = trans3d(xy$x, xy$y, rep(0, length(xy$y)) , PMAT)
                xy$x = tem$x
                xy$y = tem$y
              }

            lines(xy$x, xy$y,  col=GRIDcol, lty=2, lwd=.4, xpd=xpd )
            HL = HL +1
            tpoints[[HL]] = xy

          }
      }

    if(!is.null(TICS))
      {
        addTIX(lats, lons, PMAT=PMAT, col=GRIDcol, TICS=TICS, PROJ=PROJ)
        
      }
    
###  anotations
###  longitudes

    if(!is.null(LABS))
      {
        col = LABcol
        i = lats[1]
        j = lons[-1]
        
        xy = GLOB.XY(rep(i, length(j)), j , PROJ )
        
        if(!is.null(PMAT))
          {
            tem = trans3d(xy$x, xy$y, rep(0, length(xy$y)) , PMAT)
            xy$x = tem$x
            xy$y = tem$y
          }
        
###i = lats[2]
###xy2 = GLOB.XY(rep(i, length(lons)), lons  ) 
###angs = 180*atan2(xy$y-xy2$y,xy$x-xy2$x   )/pi
        
      
        signew = sign(j)
        degs = floor(abs(j))
        mins = round(60*(abs(j)-degs))
        degs = signew*degs
        

        

        for(k in 1:length(j))
          {
            #print(paste(sep=' ', degs[k], mins[k]))

            
            if(mins[k]==0)
              {
                alab = paste(sep="", degs[k], "\\de")
                ## text(xy$x[k], xy$y[k], labels=bquote(.(degs[k])*degree)  , xpd=xpd, pos=1, col=col)
                text(xy$x[k], xy$y[k], labels=alab  , xpd=xpd,  vfont=vfont, pos=1, col=col)
                
              }
            else
              {

                alab = paste(sep="", degs[k], "\\de", mins[k], "\\fm")
               ## text(xy$x[k], xy$y[k], labels=bquote(.(degs[k])*degree ~ .(mins[k])*minute)  , xpd=xpd, pos=1, col=col)
                text(xy$x[k], xy$y[k], labels=alab,  vfont=vfont , xpd=xpd, pos=1, col=col)
                

              }
          }
        
        i = lats[-1]
        j = lons[1]
        
        xy = GLOB.XY(i,  rep(j, length(i)), PROJ )
        if(!is.null(PMAT))
          {
            tem = trans3d(xy$x, xy$y, rep(0, length(xy$y)) , PMAT)
            xy$x = tem$x
            xy$y = tem$y
          }
        
        degs = floor(i)
        mins = round(60*(i-degs))
        for(k in 1:length(i))
          {
            if(mins[k]==0)
              {
                alab = paste(sep="", degs[k], "\\de")
              ## text(xy$x[k], xy$y[k], labels=bquote(.(degs[k])*degree)  , xpd=xpd, pos=2, col=col)

                ############   Hershey Fonts (vfont)  degrees = de   minutes = fm   seconds = sd

                
                text(xy$x[k], xy$y[k], labels=alab,  vfont=vfont  , xpd=xpd, pos=2, col=col)
              }
            else
              {
                alab = paste(sep="", degs[k], "\\de", mins[k], "\\fm")
                ## text(xy$x[k], xy$y[k], labels=bquote(.(degs[k])*degree ~ .(mins[k])*minute)  , xpd=xpd, pos=2, col=col)
                 text(xy$x[k], xy$y[k], labels=alab,  vfont=vfont , xpd=xpd, pos=2, col=col)
                
                
              }
          }
        
      }

    invisible(tpoints)
    
  }

