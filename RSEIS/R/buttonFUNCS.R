
ROT.RT<-function(nh, g)
  {
    #####  BUTTONDOC:ROT.RT:'Rotate seismogram to Radial Transverse' 
    ####  get azimuth from event and station information
    ###  these are stored in list elements in nh
    ####  under pickfile (event location)
    ###  the stations are stored in a list of station info
    gsta = g$sel

    asta = unique(nh$STNS[gsta])

    asta = asta[1]
    
    ev = nh$pickfile
    stn = nh$sta

   ###  print(ev)
   ###  print(stn)
    ### print(asta)
    
    msa = which(asta==stn$name)

    ma = gsta

    acomp  = nh$COMPS[ma]
    icomp = rep(0, length(acomp))
    icomp[acomp=="V"] = 1
    icomp[acomp=="N"] = 2
    icomp[acomp=="E"] = 3
    
    ma = ma[order(icomp)]
    
    atemp = cbind( nh$JSTR[[ma[1]]],  nh$JSTR[[ma[2]]],  nh$JSTR[[ma[3]]])
    
    daz = rdistaz(ev$LOC$lat, ev$LOC$lon, stn$lat[msa], stn$lon[msa] )
    baz = daz$baz
    
    rbaz = grotseis(baz, flip=FALSE)
    btemp  = atemp  %*%  rbaz
    ############   rbaz is N by 3
    #####   so btemp is N by 3 : take the 2nd and 3rd traces
    nh$JSTR[[ma[2]]] = btemp[,2]
    nh$JSTR[[ma[3]]] = btemp[,3]

##    nh$COMPS[ma[2]] = "R"
##     nh$COMPS[ma[3]] = "T"
##     nh$OCOMPS[ma[2]] = "RAD"
##     nh$OCOMPS[ma[3]] = "TRN"
     nh$KNOTES[ma[2]] = paste(asta, "RAD", sep=" ")
     nh$KNOTES[ma[3]] = paste(asta, "TRN", sep=" ")
       
    g$zloc = list(x=NULL, y=NULL) 
          

    g$action = "replace"
    invisible(list(NH=nh, global.vars=g))
  }

JustV<-function(nh, g)
  {
#####  BUTTONDOC:JustV:'Show only vertical'
    kix = legitpix(g$sel, g$zloc, g$zenclick)
    ypick =  kix$ypick
    ppick = kix$ppick

    sel = which(nh$COMPS=="V")
     isel = sel[1]
    
          Torigin = list(jd=nh$info$jd[isel], hr=nh$info$hr[isel],
            mi=nh$info$mi[isel],
            sec=(nh$info$sec[isel]+nh$info$msec[isel]/1000+nh$info$t1[isel]-nh$info$off[isel]))
       g$Torigin=Torigin
          g$sel = sel

    g$zloc = list(x=NULL, y=NULL)
    g$STNS = nh$STNS[sel]
    g$COMPS = nh$COMPS[sel]
    g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
  }


JustE<-function(nh, g)
  {
#####  BUTTONDOC:JustE:'Show only East' 

    kix = legitpix(g$sel, g$zloc, g$zenclick)
    ypick =  kix$ypick
    ppick = kix$ppick

    sel = which(nh$COMPS=="E")
     isel = sel[1]
          
          Torigin = list(jd=nh$info$jd[isel], hr=nh$info$hr[isel],
            mi=nh$info$mi[isel],
            sec=(nh$info$sec[isel]+nh$info$msec[isel]/1000+nh$info$t1[isel]-nh$info$off[isel]))
 g$Torigin=Torigin
          g$sel = sel

    g$zloc = list(x=NULL, y=NULL)
    g$STNS = nh$STNS[sel]
    g$COMPS = nh$COMPS[sel]
    g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
  }


JustN<-function(nh, g)
  {
#####  BUTTONDOC:JustN:'Show only North' 
    kix = legitpix(g$sel, g$zloc, g$zenclick)
    ypick =  kix$ypick
    ppick = kix$ppick

    sel = which(nh$COMPS=="N")
     isel = sel[1]
          
          Torigin = list(jd=nh$info$jd[isel], hr=nh$info$hr[isel],
            mi=nh$info$mi[isel],
            sec=(nh$info$sec[isel]+nh$info$msec[isel]/1000+nh$info$t1[isel]-nh$info$off[isel]))
 g$Torigin=Torigin
          g$sel = sel

    g$zloc = list(x=NULL, y=NULL)
    g$STNS = nh$STNS[sel]
    g$COMPS = nh$COMPS[sel]
    g$action = "replot"
    invisible(list(NH=nh, global.vars=g))
  }
