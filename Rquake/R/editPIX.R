editPIX<-function(nh, g )
  {

    if(is.null(g$WPX))
      {
        g$zloc = list(x=NULL, y=NULL)
        g$action="replot"
        invisible(list(global.vars=g))
        
      }
    else
      {
        
        twpx  = g$WPX
        if(g$zenclick>=2)
          {
            zappa = match(g$KLICK, g$BLABS)
            col = g$colpix[which(g$pnos=="GPIX")]
            kix = RSEIS::legitpix(g$sel, g$zloc, g$zenclick)
            ypick =  kix$ypick
            ppick = kix$ppick
            

            if(length(ypick)>0)
              {
                
                onewpx = RSEIS::cleanWPX()
                
                ipick = g$sel[ypick]
                

                
                for(iz in 1:length(ypick))
                  {
                    
                    i1 = ipick[iz]
                    i2 = ypick[iz]

                    ycol = g$colpix[zappa]
                    if(is.na(ycol)) { ycol = rgb(0,0,1) }
                    
                    err = NA
                    res = 0

                    asec = nh$info$sec[i1] + nh$info$msec[i1]/1000 + nh$info$t1[i1] -
                      nh$info$off[i1] + ppick[iz]


                    sta1 =  nh$STNS[i1]
                    comp1 = nh$COMPS[i1]
                    pic1 = RSEIS::recdate(nh$info$jd[i1], nh$info$hr[i1], nh$info$mi[i1],
                      asec, yr=nh$info$yr[i1])

                   #### print(data.frame(pic1))

                    ####  match the station name and find the nearest pick
                    w1 = which(twpx$name==sta1 & twpx$comp==comp1)

                    if(length(w1)<1)next
                    
                    sek = abs( RSEIS::secdifL(pic1, twpx) )

                    iw1 = which.min( sek[w1] )
                    
                    getridofit  = w1[iw1]

                    twpx = RSEIS::deleteWPX(twpx,  getridofit) 
                    ## 
                  }
              }
          }

        g$WPX = twpx

      }

    
    
    g$zloc = list(x=NULL, y=NULL)
    g$action="replot"
    invisible(list(global.vars=g))
    

  }
