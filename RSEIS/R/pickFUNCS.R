
WPIX<-function(nh, g)
  {
    #####  BUTTONDOC:WPIX:'Window picks - pairs of picks on each trace'
    #####   WPIX are introduced as pairs of pix
    if(g$zenclick>=2)
      {
        zappa = match(g$KLICK, g$BLABS)
        
        kix = legitpix(g$sel, g$zloc, g$zenclick)
        ypick =  kix$ypick
        ppick = kix$ppick
        LIX = floor(length(ypick)/2)

        
############   proceed only if have legitimate picks
        if(LIX>0)
          {
            azap = "WPIX"
            kzap = "W"
            ipick = g$sel[ypick]
            
print(paste(sep=" ", "DUMP WPIX", zappa, azap, kzap , ppick , ypick,ipick)) 
            
            for(LZ in 1:LIX)
              {
                iz = (LZ)*2-1
                
                g$NPX = g$NPX+1
                Nn = names(g$WPX)
               g$WPX =rbind(g$WPX, rep(NA, length(Nn)))
                
                i1 = ipick[iz]
                ycol = g$colpix[zappa]
                if(is.na(ycol)) { ycol = rgb(0,0,1) }
                err = NA
                res = ppick[iz+1]-ppick[iz]
                print(paste(i1, iz, ppick[iz], kzap, res, err, ycol))
                g$WPX =  pickhandler(i1=i1, ppick=ppick[iz], kzap=kzap, res=res, err=err, ycol=ycol, NPX=g$NPX, WPX=g$WPX, NH=nh)
                g$NADDPIX = g$NADDPIX+1
                
                ## 
              }
###PLOT.ALLPX(Torigin, STNS, COMPS, WPX, PHASE=PHASE, FORCE=forcepix, cex=pcex)
            
          }
      }
    g$zloc = list(x=NULL, y=NULL) 
    g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
  }

##########################
NOPIX<-function(nh, g)
  {
    #####  BUTTONDOC:NOPIX:'Turn off All Picks (set onoff to zero)' 
    g$WPX$onoff = rep(-1, length(g$WPX$onoff))
    g$zloc = list(x=NULL, y=NULL) 
    g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
  }
##########################
noPS<-function(nh, g)
  {
    #####  BUTTONDOC:NOPIX:'Turn off 1 station Picks (set onoff to zero)'
    print(g$sel)
    
    usta = unique( nh$STNS[g$sel] )
    m = which (g$WPX$name  %in% usta )
    
    
    g$WPX$onoff[m] = rep(-1, length(m) )
    g$zloc = list(x=NULL, y=NULL) 
    g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
  }
##########################


REPIX<-function(nh, g)
  {
    #####  BUTTONDOC:REPIX:'Turn Picks back on' 
    g$WPX$onoff[g$WPX$onoff==(-1)] = 0
    g$zloc = list(x=NULL, y=NULL) 
   g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
  }
##########################



##########################
FILLPIX<-function(nh, g)
  {
    #####  BUTTONDOC:FILLPIX:'Pick line spans vertical window'
    g$fillpix = !g$fillpix 
    g$zloc = list(x=NULL, y=NULL) 
   g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
  }
##########################

##########################
RIDPIX<-function(nh, g)
  {
#####  BUTTONDOC:RIDPIX:'Remove picks' 
    tol = .1
############   find picks near the clicks and remove...temporarily
    if(g$zenclick>=2)
      {
        zappa = match(g$KLICK, g$BLABS)
        col = g$colpix[which(g$pnos=="YPIX")]
        kix = legitpix(g$sel, g$zloc, g$zenclick)
        ypick =  kix$ypick
        ppick = kix$ppick

############   proceed only if have legitimate clicks
        if(length(ypick)>0)
          {
            
            ipick = g$sel[ypick]
            
####
            
            
            for(iz in 1:length(ypick))
              {
                
                i1 = ipick[iz]
                i2 = ypick[iz]
              ##  print(paste(iz, i1, i2))

                
                asec = nh$info$sec[i1]+nh$info$msec[i1]/1000+nh$info$t1[i1]-nh$info$off[i1]+ppick[iz]
                pic1 = recdate(nh$info$jd[i1], nh$info$hr[i1], nh$info$mi[i1], asec)


             ##  print(pic1)
                
                ds = abs(secdifL(g$WPX, pic1))
                

            ##    print(ds)

                if( any(ds<=tol) )
                  {
                    
                    irid = which.min(ds)

                    
                  ##   print(irid)
                    samesta = g$WPX$name[irid]==nh$STNS[i1] & g$WPX$comp[irid]==nh$COMPS[i1]
                    if(samesta) g$WPX$onoff[irid] = -1
                  }

                
                
             ##   g$WPX =  pickhandler(i1=i1, ppick=ppick[iz], kzap=kzap, res=res, err=NA, ycol=ycol, NPX=g$NPX, WPX=g$WPX, NH=nh)
               ##  g$ADDPIX =  pickhandler(i1=i1, ppick=ppick[iz], kzap=kzap, res=res, err=NA, ycol=ycol, NPX=g$NPX, WPX=g$WPX, NH=nh)
               ## g$NADDPIX = g$NADDPIX+1
                
                ## 
              }
###PLOT.ALLPX(Torigin, STNS, COMPS, WPX, PHASE=PHASE, FORCE=forcepix, cex=pcex)
            
          }
      }
    g$zloc = list(x=NULL, y=NULL) 
    g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
    

    
  }
##########################


##########################
SEEPIX<-function(nh, g)
  {
    #####  BUTTONDOC:SEEPIX:'print picks to screen'
    options(width=180)
    print(g$WPX)
    
    #### write.table(g$WPX)
    g$zloc = list(x=NULL, y=NULL) 
   g$action = "donothing"
    invisible(list(NH=nh, global.vars=g))
  }
##########################



iNEXT<-function(nh, g)
  {
    #####  BUTTONDOC:iNEXT:'Do Nothing' 
    g$zloc = list(x=NULL, y=NULL)
    g$action = "break"
    invisible(list(NH=nh, global.vars=g))
   
  }


####################################
Ppic<-function(nh, g)
  {
#####  BUTTONDOC:Ppic:'P-wave pick' 

    
    g = pADDPIX(nh, g, "P")
    nh$WPX = g$WPX
    g$zloc = list(x=NULL, y=NULL)
    g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
    
  }
#####
Spic<-function(nh, g)
  {
#####  BUTTONDOC:Spic:'S-wave pick' 

    
    g = pADDPIX(nh, g, "S")
    nh$WPX = g$WPX
    g$zloc = list(x=NULL, y=NULL)
    g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
    
  }
#####
Apic<-function(nh, g)
  {
    #####  BUTTONDOC:Apic:'Acoustic wave pick' 

    g = pADDPIX(nh, g, "A")
    nh$WPX = g$WPX
    g$zloc = list(x=NULL, y=NULL)
    g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
    
  }
####################################
POLSWITCH<-function(nh, g, dir)
  {
#####  BUTTONDOC:POLSWITCH:'Switch Polarity'
    zappa = match(g$KLICK, g$PADDLAB)
    azap = g$PADDLAB[zappa]
###     print(paste(sep=" ", "My PICK", dir, azap, zappa))

    kix = legitpix(g$sel, g$zloc, g$zenclick)
    ypick =  kix$ypick
    ppick = kix$ppick
    
    
    
    ipick = g$sel[ypick]

    ##   print(ipick)
    
    if(length(ipick)<1)
      {
        ipick = g$sel[which(nh$COMPS[g$sel]=="V" )]
      }
    ##   print(g$sel)
   ##  print(nh$STNS)
   ##  print(nh$COMPS)
    
    print(paste(sep=" ", "PICK=", nh$info$yr[ipick],
                nh$info$jd[ipick], nh$info$hr[ipick],
                nh$info$mi[ipick], "sta=", nh$STNS[ipick],
                "comp=", nh$COMPS[ipick] ))

    m = match(g$STNS[ypick],g$UNIsta)

    jj = floor((g$zloc$y[g$zenclick-1])/g$du)
    
    iseek = which(g$WPX$phase=="P" & g$WPX$name==nh$STNS[ipick] &  g$WPX$comp==nh$COMPS[ipick])
####  print(paste(sep=" ", "ISEEK",  iseek, length(iseek) ))
    
    if(length(iseek)==1)
      {
        wNPX = iseek
        
        g$WPX$pol[wNPX]=dir
        
      }
    else
      {
        print(paste(sep=" ", "NO MATCH FOUND ISEEK",  iseek, length(iseek) ))
        print("Click in a panel first, then select polarity")

      }
    return(g)

  }

#####  polarity determinations
Pup<-function(nh, g)
  {
#####  BUTTONDOC:Pup:'Up Polarity' 

    g = POLSWITCH(nh, g, "U")
    g$zloc = list(x=NULL, y=NULL)
    g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
    
  }
#####
Pnil<-function(nh, g)
  {
#####  BUTTONDOC:Pnil:'NUll Polarity'
    g = POLSWITCH(nh, g, "N")
    
    g$zloc = list(x=NULL, y=NULL)
    g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
    
  }
#####
Pdown<-function(nh, g)
  {
#####  BUTTONDOC:Pdown:'Down Polarity'
    g = POLSWITCH(nh, g, "D")
    g$zloc = list(x=NULL, y=NULL)
    g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
    
  }
####################################
