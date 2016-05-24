Pick3<-function(nh, g)
  {

     kix = RSEIS::legitpix(g$sel, g$zloc, g$zenclick)
          ypick =  kix$ypick
          ppick = kix$ppick
                
          if(length(ppick)>0)
            {
            
              ipick = g$sel[ypick]

              ipick = ipick[length(ipick)]
              
             ## cat(paste(sep=" ", ypick, ipick), sep="\n")
             ## print(ipick)
              ##
              
              ma = which(!is.na(match( nh$STNS, nh$STNS[ipick])))

              
              ##########   sort so Vertical is on top and then North and East
              acomp  = nh$COMPS[ma]
              icomp = rep(0, length(acomp))
              icomp[acomp=="V"] = 1
              icomp[acomp=="N"] = 2
              icomp[acomp=="E"] = 3

              ma = ma[order(icomp)]

              
####  print(cbind(nh$STNS[ma], nh$COMPS[ma]))

              
              if(is.null(g$Pickdev))
                {
               ####   X11(width = 12, height = 7)
                  RSEIS::screens(2)
                  devl = dev.list()
                  iw =  which(g$MAINdev!=devl)
                  
                  g$Pickdev = devl[iw[1]]
                   dev.set(g$Pickdev)
                }
              else
                {
                 #### devl = dev.list()
                ####  jsc = 2-length(devl)
                ####  if(jsc>0) { X11(width = 12, height = 7); Pickdev = dev.cur() }
                  dev.set(g$Pickdev)
                }

              if(g$zenclick>2)
                {

                  pickwin = range( c(g$zloc$x[(g$zenclick-1)], g$zloc$x[(g$zenclick-2)]))
                  
                }
              else
                {
                  pickwin = g$WIN

                }
              
              PICKLAB = c("DONE", "Ppic", "Spic", "ZOOM.out","ZOOM.in", "REFRESH", "RESTORE",
                "FILT", "UNFILT", "Pinfo", "WINFO", "ROT.RT")

              PLAB=c(  "Apic",  "Pup", "Pdown", "Pnil", "AUTOP",
                "NOPIX", "EDIX", "REPIX")

              stit = nh$STNS[ma[1]]
              ##  SWP = selAPX(WPX,  nh$STNS[ma[1]], icomp=NULL )

              ##   print(data.frame(SWP))
              ##   SWP = rectifyAPX(SWP)
              ##
              ## print(SWP)

              
              newpicks = RSEIS::swig(nh, APIX=g$WPX, sel=ma, WIN=pickwin,
                STDLAB=PICKLAB ,PADDLAB=PLAB, PHASE=1   ,
                SHOWONLY = FALSE, TIT=stit)


              
              if(length(newpicks$g$WPX)>=1)
                {
                  if(!is.null(newpicks$g$WPX))
                    {
                      g$WPX = newpicks$g$WPX
                    }
                }
              ##  
              
              ##
####    print(cbind(WPX$name, WPX$comp, WPX$phase, WPX$onoff))
              g$NPX = length(g$WPX$name)
####                print(paste(sep=' ', "DONE with PICKWIN", g$NPX))
              dev.set( g$MAINdev)

              
              
            }
          
          g$zloc = list(x=NULL, y=NULL) 
          

   g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
  }
