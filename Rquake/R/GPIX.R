GPIX<-function(nh, g)
  {
    if(g$zenclick>=2)
      {
        zappa = match(g$KLICK, g$BLABS)
        col = g$colpix[which(g$pnos=="GPIX")]
        kix = RSEIS::legitpix(g$sel, g$zloc, g$zenclick)
        ypick =  kix$ypick
        ppick = kix$ppick

       ##### print(ppick)
############   proceed only if have legitimate picks
        if(length(ypick)>0)
          {

            onewpx = RSEIS::cleanWPX()
            azap = "GPIX"
            kzap = "G"
            ipick = g$sel[ypick]
            
#### print(paste(sep=" ", "DUMP YPIX", zappa, col, azap, kzap , ppick , ypick,ipick)) 
            
            for(iz in 1:length(ypick))
              {
              ####  g$NPX = g$NPX+1
              ####  Nn = names(g$WPX)
             ####  g$WPX =rbind(g$WPX, rep(NA, length(Nn)))
                
                i1 = ipick[iz]
                i2 = ypick[iz]

                ycol = g$colpix[zappa]
                if(is.na(ycol)) { ycol = rgb(0,0,1) }
                
                err = NA
                res = 0

                asec = nh$info$sec[i1] + nh$info$msec[i1]/1000 + nh$info$t1[i1] -
                  nh$info$off[i1] + ppick[iz]
                pic1 = RSEIS::recdate(nh$info$jd[i1], nh$info$hr[i1], nh$info$mi[i1],
                  asec, yr=nh$info$yr[i1])
                
                ##    
                ##   g$WPX =  pickhandler(i1=i1, ppick=ppick[iz], kzap=kzap, res=res, err=NA, ycol=ycol, onoff=2, NPX=g$NPX, WPX=g$WPX, NH=nh)
                ##  g$ADDPIX =  pickhandler(i1=i1, ppick=ppick[iz], kzap=kzap, res=res, err=NA, ycol=ycol, NPX=g$NPX, WPX=g$WPX, NH=nh)
                
                xpx =   RSEIS::setWPX(phase = kzap, col = ycol, yr =pic1$yr , jd =pic1$jd ,
                  hr = pic1$hr, mi = pic1$mi, sec = pic1$sec, dur = 0, name = nh$STNS[i1],
                  comp = nh$COMPS[i1], dispcomp =nh$COMPS[i1] , onoff = 1)

                onewpx = RSEIS::catWPX(onewpx, xpx)
                
                ## 
              }
###PLOT.ALLPX(Torigin, STNS, COMPS, WPX, PHASE=PHASE, FORCE=forcepix, cex=pcex)

            if(any(is.na(onewpx$name)))
              {
                idpx = which(is.na(onewpx$name))
                
                onewpx =  RSEIS::deleteWPX(onewpx, idpx)
              }

            g$WPX = RSEIS::catWPX(g$WPX,onewpx )

            g$NPX = length(g$WPX$sec)
             
            g$PHASE = unique( c(g$PHASE, "Y") )
          }
      }

    ## print(g$PHASE)

    
    g$zloc = list(x=NULL, y=NULL) 
    g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
  }
############################################################
