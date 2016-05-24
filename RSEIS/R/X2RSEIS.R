X2RSEIS<-function(nh, g)
  {
#####  BUTTONDOC:X2SAC:'data extraction '
    kix = legitpix(g$sel, g$zloc, g$zenclick)
    ypick =  kix$ypick
    ppick = kix$ppick
    
    if(length(ppick)>0)
      {

        ipick = g$sel[ypick]

        ipick = ipick[length(ipick)]
        
        ## cat(paste(sep=" ", ypick, ipick), sep="\n")
        ## print(ipick)
        ##
       if(g$zenclick>2)
              {
                pickwin = range( c(g$zloc$x[(g$zenclick-1)], g$zloc$x[(g$zenclick-2)]))
                
              }
            else
              {
                pickwin = g$WIN

              }
          
        GH =  WINGH(nh,  WIN = pickwin )


        pstamp = Zdate(GH$info)
        outfile = pstamp[1]
        cat(paste("Creating RSEIS GH file:", outfile ), sep="\n")
        
        
        save(file=paste(outfile, "GH.RDATA", sep="."), GH)

        
         g$zloc = list(x=NULL, y=NULL) 
        
        g$action="donothing"
        invisible(list(RETX =NULL , global.vars=g))
      
      }
else
      {
        cat("X_R WARNING: no window or trace has been selected:", sep="\n")
        RETX=NULL
        g$zloc = list(x=NULL, y=NULL) 
        
        g$action="donothing"
        invisible(list(global.vars=g))
        
      }

  }

