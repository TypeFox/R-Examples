pADDPIX<-function(nh, g, phase)
  {
    #####  BUTTONDOC:pADDPIX:'add picks to list'
    
   zappa = match(g$KLICK, g$PADDLAB)
    azap = g$PADDLAB[zappa]
    print(paste(sep=" ", "My PICKIN", azap, zappa))

    kix = legitpix(g$sel, g$zloc, g$zenclick)
    ypick =  kix$ypick
    ppick = kix$ppick
    
###   print(paste(sep=" " , "WIN=",sloc$x))
    
###        abline(v=ppick, col=4)
    
    ipick = g$sel[ypick]


   ####################   here I take the first click - but is that right?
   ####  is this because I am forcing only one P-wave arrival?
   ####   that does not make sense.....


   
    ipick = ipick[1]

    
 ###   print(paste(sep=" ", "ZPICK=", nh$info$yr[ipick],
  ###              nh$info$jd[ipick], nh$info$hr[ipick],
 ###               nh$info$mi[ipick], "sta=",
 ###               nh$STNS[ipick], "comp=", nh$COMPS[ipick] ))

    m = match(g$STNS[ypick],g$UNIsta)
   
###  Upix[[m]]$x  = ppick
    
###   PPIX(list(x=zloc$x[zenclick-1], y=zloc$y[zenclick-1]), YN=NSEL, col=3, lab="P")
    jj = floor((g$zloc$y[g$zenclick-1])/g$du)
    
    if((g$zenclick==2))
      {
        asec = nh$info$sec[ipick]+nh$info$msec[ipick]/1000+
          nh$info$t1[ipick]-nh$info$off[ipick]+ppick[g$zenclick-1]
        err = 0.05
      }
    else
      {
        asec = nh$info$sec[ipick]+nh$info$msec[ipick]/1000+
          nh$info$t1[ipick]-nh$info$off[ipick]+ppick[g$zenclick-2]
        bsec = nh$info$sec[ipick]+nh$info$msec[ipick]/1000+
          nh$info$t1[ipick]-nh$info$off[ipick]+ppick[g$zenclick-1]
        err =  abs(bsec-asec)
      }

###########   this looks like a bug./....
 ###   print(paste(  nh$STNS[ipick],   nh$COMPS[ipick], phase))

  
    iseek = which(g$WPX$phase==phase & g$WPX$name==nh$STNS[ipick]
      &  g$WPX$comp==nh$COMPS[ipick])

###  print(paste(sep=" ", phase, nh$STNS[ipick], nh$COMPS[ipick],
###              "ISEEK",  iseek, length(iseek) ))

   onepx = cleanWPX()

   onepx$phase=phase
   
   onepx$yr=nh$info$yr[ipick]
   onepx$mo= nh$info$mo[ipick]
   onepx$dom=nh$info$dom[ipick]
   onepx$jd=nh$info$jd[ipick]
   onepx$hr= nh$info$hr[ipick]
   onepx$mi=nh$info$mi[ipick]
   onepx$col=g$specpix.col[4]
   onepx$sec=asec
   onepx$err=err
   onepx$onoff = 1
   
   
   if(length(iseek)==1)
     {
##############   replace the pick with the current pick
       wNPX = iseek
       onepx$tag = g$WPX$tag[wNPX]
       onepx$name = g$WPX$name[wNPX]
       onepx$comp = g$WPX$comp[wNPX]
       onepx$c3 = g$WPX$c3[wNPX]

       g$WPX =  replaceWPX(g$WPX, onepx, wNPX)
       
     }
   else
     {
       onepx$name=nh$STNS[ipick]
       onepx$comp=nh$COMPS[ipick]
       onepx$c3=nh$OCOMPS[ipick]
       onepx$tag=paste(sep=".",nh$STNS[ipick],  nh$OCOMPS[ipick])
###############   add a new pick to WPX list
       g$WPX = catWPX(g$WPX,onepx )
       }


   
   g$NPX = length(g$WPX$sec)
             
   g$PHASE = unique( c(g$PHASE, "Y") )
   
   
   g$NADDPIX = 3
###

    
   return(g)

  }

