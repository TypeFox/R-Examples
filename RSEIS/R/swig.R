`swig` <-function(GH, sel=1:length(GH$dt), ORD=NULL, WIN=NULL, APIX=NULL, PHASE=NULL,
                  STDLAB=NULL,
                  PADDLAB=NULL, TEMPBUT=NULL, SHOWONLY=FALSE,
                  CHOP=FALSE, TIT="", pts=FALSE, forcepix=FALSE,
                  pcex=0.7, SCALE=1, ilocstyle=1,
                  velfile="", stafile="", LOC=NULL,
                  prefilt=list(fl=.2, fh=15,  type="HP", proto="BU"), filters=NULL, YAX = 1  , xtickfactor = 1  )
{


###  NEW version of PICK.GEN remove the "global variables"

###  instead pass variables through a list:  this way we
###    can have external function defined.
  

#########  velfile is a 1D velocity file
#########  stafile is a file station locations
  
  if(missing(WIN)) { WIN = NULL }
  if(missing(sel)) { sel = 1:length(GH$dt)}
  if(missing(ORD)) { ORD = NULL }


  if(missing(xtickfactor)) { xtickfactor = 1 }
    if(missing(YAX)) { YAX = 1 }
  
  
###   if(missing(APIX)) { APIX = NULL}  else { if(!exists(deparse(substitute(APIX)))) { print("WARNING: NO WPX"); APIX=NULL} }

  if(missing(APIX)) { APIX = NULL}
  
  if(missing(PHASE)) {
    PHASE = NULL
    if(!is.null(APIX))
      {
        PHASE = unique(APIX$phase)

      }

    
  }


  
  
  if(missing(SHOWONLY)) { SHOWONLY=FALSE}
  if(missing(CHOP)) { CHOP=FALSE }
  
  if(missing(STDLAB)) {
  ##  if(exists('STDLAB.DEFAULT'))
   ##     {
   ##      STDLAB =  STDLAB.DEFAULT
  ##      }
  ##    else
   ##     {
      STDLAB = c("REPLOT","DONE", "SELBUT", "PSEL","LocStyle",
        "ZOOM.out", "ZOOM.in", "LEFT", "RIGHT", "RESTORE",  "Pinfo","WINFO",
                      "Xwin",    "XTR", "SPEC", "SGRAM" ,"WLET",
        "FILT", "UNFILT", "SCALE", "Postscript", "PREV")
  }
  
  if(missing(PADDLAB)) {
   ##   if(exists('PADDLAB.DEFAULT'))
    ##    {
    ##      PADDLAB=PADDLAB.DEFAULT
     ##   }
   ##   else
    ##    {
      PADDLAB=c("YPIX", "NOPIX", "REPIX")
    

  }

  if(missing(TEMPBUT)) { TEMPBUT=NULL } 
  
  if(missing(TIT)) { TIT=NULL }
  if(missing(pts)) { pts=FALSE }
  if(missing(forcepix)) { forcepix=FALSE }
  if(missing(pcex)) { pcex = 0.7 }
  
  if(missing(velfile)) {
    if(!is.null(GH$velfile)) {velfile=GH$velfile } else { 
      velfile=NULL }
  }
  
  if(missing(stafile)) {
    if(!is.null(GH$stafile)) {stafile=GH$stafile } else {  
      stafile=NULL }
  }

  
  if(missing(LOC)) { LOC=NULL }

  if(missing(SCALE)) {  ScaleFACT = 1 } else {  ScaleFACT = SCALE }
  if(missing(prefilt)) {
    prefilt=NULL
  }

  if(missing(filters)) {
    filters = NULL
    lastfilter=NULL
  }
  else
    {
      lastfilter=1
    }

  if(missing(ilocstyle)) { ilocstyle=1 }
  
  
 

  
  if(is.logical(sel)) { sel = which(sel) } 
  
  if( is.null(sel) ) { sel = 1:length(GH$dt) }

  SEL.ORIG = sel

  SHIFT.ORIG = NULL
  ASHIFT = SHIFT.ORIG  
  
  mark = FALSE


###################   copy of GH file ..... do I need this?
  if(CHOP==TRUE)
    {
      if(!is.null(WIN))
        {
          NH = CHOP.SEISN(GH, sel , WIN=WIN)
        }
      else
        {
          
          NH = GH
        }
      
      WIN = c(0, NH$dt*length(NH$JSTR[[1]]))
      
    }
  else
    {
      
      NH = GH
    }


  if( identical(is.na(match("NOPIX",  PADDLAB)), TRUE)) { PADDLAB = c(PADDLAB, "NOPIX") }
  if( identical(is.na(match("REPIX",  PADDLAB)), TRUE)) { PADDLAB = c(PADDLAB, "REPIX") }

  if(!exists('STDLAB'))
    {
      STDLAB = c("REPLOT", "DONE",  "ZOOM.out", "RESTORE", "SavePF", 
        "PickWin", "XTR", "SPEC", "SGRAM" ,"WLET")
    }


  defaultLAB = c("REPLOT", "DONE", "PREV", "SELBUT", "PSEL", "LocStyle", "ZOOM.out", "ZOOM.in", "LEFT", "RIGHT", "RESTORE",  "Pinfo","WINFO",
                          "XTR", "SPEC", "SGRAM" ,"WLET", "FILT", "UNFILT", "SCALE", "Postscript")
   
  
  stdlab =  STDLAB

  Pickdev = NULL
  Mapdev = NULL

  tempbuttons = NULL

  BLABS = c(stdlab, PADDLAB)
  NLABS = length(BLABS)
  NOLAB = NLABS +1000



 
##  print(BLABS)
  

fixedbuttons = c("REPLOT", "DONE", "QUIT", "SELBUT" )

BLABS = unique(c(fixedbuttons, BLABS))

##  print(BLABS)
  

OTHERbuttons = c("NEXT", "PREV","HALF","S1", "S2", "MARK", "DOC", "RESTORE",
"ZOOM.out", "ZOOM.in", "LEFT",  "RIGHT", "SCALE",  "SHOWALL",
  "SHOWSEL", "saveFN", "FLIP", "TR_INFO", "Postscript","PNG",
  "AUTOP","AUTOPALL", "DETECT",  "MAP", "XTR",  "SIG",
"SPEC.old",  "ASCII",  "AMPL",   "TRNAMPL",  "SPEC", "SGRAM",  "WLET",
  "FILT","UNFILT",  "BRUNE",   "DETAIL",   "PTS", "MMARKS",   "PMOT",
  "STERNET", "GTAZI", "ENVLP", "WINFO", "Pinfo", "XCOR",
"PHLAG", "3COMP", "Predict1D",   "SavePF", "SavePIX", "LQUAKE",
 "PickWin", "Ppic", "Spic", "Apic", "Pup",   "Pdown",
     "Pnil",
"YPIX", "WPIX", "RMS","EDIX","FILLPIX",
 "NOPIX", "REPIX",
       "ADDBUTTS")


 ALLlabs = unique( c(BLABS, OTHERbuttons) )


  litecolors = c( "peachpuff2",      "darkolivegreen2", "slategray1" ,     "lightgoldenrod1",
  "darkseagreen3",   "lavenderblush2" , "slategray2" ,     "thistle1"      , 
  "cadetblue2"  ,    "lemonchiffon3"  )


  
  ##  match("", BLABS)

  RETX =  NULL


  
  somecolors = c("black", "darkmagenta", "forestgreen", "blueviolet",
    "tan3", "lightseagreen", "deeppink", "cyan3", "bisque3",
    "magenta1", "lightsalmon3", "darkcyan", "darkslateblue",
    "chocolate4", "goldenrod4", "mediumseagreen")

  APAL=c("black","darkmagenta","forestgreen","blueviolet",
    "tan3","lightseagreen","deeppink","cyan3","bisque3",
    "darkcyan","darkred","firebrick3","rosybrown4","royalblue3",
    "darkblue","red2","violetred3","springgreen3","darkorange4",
    "palevioletred3","mediumpurple3","tomato","dodgerblue1",
    "olivedrab4","yellow4","pink4")

  ##   pnos = grep("PIX", BLABS)
  pnos = c( grep("PIX", BLABS), grep("pik", BLABS))
  colabs = rep(1,length(BLABS))
  colabs[BLABS=="PickWin"] = 'red'
  colabs[BLABS=="REPLOT"] = 'red'
  
  colabs[pnos] = somecolors[seq(from=2, length=length(pnos))]
  
  colpix = somecolors[seq(from=2, length=length(pnos))]
  
  pchlabs = rep(4,length(BLABS))
  pchlabs[pnos] = seq(from=15, length=length(pnos))
  pchlabs[BLABS=="REPLOT"] = 5
  
  specpix =     c("P", "S", "A", "P1", "Y", "W")
  specpix.col = c("violetred", "darkgoldenrod", "blueviolet" , "darkmagenta", "dodgerblue1", "darkcyan")
  
  
  NSEL = length(NH$dt[sel])

  if(is.null(APIX)==TRUE)
    {
      WPX = cleanWPX()
      
      NPX = 0
    }
  else
    {
      ## print("reading in pickfile")
      ##  
      WPX = APIX
      ##  WPX = data.frame(WPX, stringsAsFactors = FALSE)
      NPX = length(WPX$sec)

     ##  print(paste(sep=' ', "read in pickfile",NPX))
      ## print(WPX)
    }
  
  RIDPIX = list()
  ADDPIX = list()
  NADDPIX = 0
  
  BRUNKOUNT = 0
  BRUNINFO = list()

  DETLKOUNT = 0
  DETLINFO = list()

  
  STNS = NH$STNS[sel]
  COMPS = NH$COMPS[sel]

###   print(STNS)
###   print(COMPS)
  
  UNIsta = unique(STNS)
  
  NUNI = length( UNIsta)

  if( identical(NH$pcol , "AUTO") |  is.null( NH$pcol )  )
    {

      
      pcols = rep(rgb(0,0,0), length(NH$dt) )
      pcols[c(grep("1", COMPS), grep("I", COMPS), grep("LD", COMPS) )] = rgb(0,.4,0)
      pcols[c(grep("4", COMPS), grep("V", COMPS), grep("Z", COMPS), grep("v", COMPS), grep("z", COMPS)   )] = rgb(0.4,0,0)
      
      pcols[c(grep("J", COMPS), grep("K", COMPS))] = rgb(0,0,0.4)
      
    }
  else
    {

      pcols = NH$pcol

      if(is.numeric(pcols))
        {
          pcols = APAL[1+((pcols-1) %% length(APAL))]
          

        }

    }

###  print(pcols)

###   want the sorting of comps to be Vertical North East always


  if(!is.null(ORD))
    {
      STNS = STNS[ORD]
      COMPS = COMPS[ORD]
    }

  ##   print("*************    check stations and comps****** ")
  ##  print(STNS)
  ##  print(COMPS)
  
  
  du = 1/NSEL

###  pix label size
  
  isel = sel[1]
  
  Torigin = list(yr=NH$info$yr[isel], jd=NH$info$jd[isel], hr=NH$info$hr[isel],
    mi=NH$info$mi[isel],
    sec=(NH$info$sec[isel]+NH$info$msec[isel]/1000+NH$info$t1[isel]-NH$info$off[isel]))

###  print(Torigin)
  
###  print(sel)
###  print(STNS)
###  print(COMPS)
### print(NH$KNOTES[sel])
###  print(NSEL)
###  print(NH$dt[sel])
###  print(pcols[sel])
  LASTwin = WIN
##########################################
###   GLOBAL VARS to be sent to subroutines and functions
  global.vars = list(sel=sel,
    SEL.ORIG=sel,
    ilocstyle = ilocstyle,
    iloccol = rgb(1,0.6, 0.6),
    ilocnum = 1,
    SHIFT.ORIG = NULL,
    ASHIFT = SHIFT.ORIG ,
    mark = FALSE,
    STDLAB =  STDLAB,
    stdlab =  STDLAB,
    PADDLAB = PADDLAB,
    SUBTIT=NA,
    TIT = TIT,
    Pickdev = NULL,
    Mapdev = NULL,
    MAINdev=NULL,
    tempbuttons = NULL,

    BLABS = BLABS ,
    NLABS = length(BLABS),
    NOLAB = NOLAB,
    
    RETX =  NULL,
    somecolors =somecolors ,
    APAL=APAL,
    pcols = pcols,
    pnos = pnos,
    colabs = colabs,
    colpix = somecolors[seq(from=2, length=length(pnos))],
    pchlabs = pchlabs,
    specpix =  specpix  ,
    specpix.col = specpix.col,
    Torigin = Torigin,
    NPX = NPX,
    NSEL = length(NH$dt[sel]),
    du=du,
    STNS=STNS,
    COMPS=COMPS,
    UNIsta = UNIsta,
    NADDPIX=NADDPIX,
    ADDPIX=ADDPIX,
    RIDPIX=RIDPIX,
    WPX=WPX, 
    PHASE=PHASE,
    forcepix=forcepix,
    fillpix=FALSE,
    srtpix=0,
    polspix=TRUE,
    pcex=pcex,
    xtickfactor = xtickfactor,
    YAX =  YAX,
 
    filters = filters,
    lastfilter = lastfilter,
    
    ScaleFACT=ScaleFACT,
    pts = pts,
    action = "continue",
    WIN =WIN,
    LASTwin = LASTwin,
    KLICK = NULL,
    thebuts = FALSE
    )



   ##################
   ################## 
  YNreplot<-function()
    {
      sel = global.vars$sel
      
      YN = PLOT.SEISN(NH, WIN=global.vars$WIN, dt=NH$dt[sel],
        sel=global.vars$sel, sfact=global.vars$ScaleFACT ,
        notes=NH$KNOTES[sel], COL=global.vars$pcols, TIT=global.vars$TIT,
        SHIFT=global.vars$ASHIFT , pts=global.vars$pts, YAX=global.vars$YAX,   xtickfactor = global.vars$xtickfactor )

      if(!is.na(global.vars$SUBTIT)) title(sub=global.vars$SUBTIT)

      YN$STNS = NH$STNS[sel]
      YN$COMPS = NH$COMPS[sel]
      YN$notes = NH$KNOTES[sel]
      
      if(global.vars$NPX>0)
        {
          
          
          swig.ALLPX(global.vars$Torigin, YN$STNS, YN$COMPS, global.vars$WPX,
                     PHASE=global.vars$PHASE,  POLS=global.vars$polspix,
                     FILL=global.vars$fillpix , FORCE=global.vars$forcepix,
                     cex=global.vars$pcex, srt=global.vars$srtpix)
##### PLOT.WPX(Torigin, STNS, COMPS, WPX, FORCE=forcepix)
        }
      invisible(YN)
    }
##################
   ################## 
  ##################   set the initial data structure (so it can be retrieved)
  
  OLDH=NH


  #####  initial filter to be applied?
  if(!is.null(prefilt))
    {
      NH = FILT.SEISN(GH, sel = 1:length(GH$JSTR), FILT = prefilt, TAPER = 0.1, POSTTAPER = 0.1)
    }


  ##############   plot the data and set the menu buttons


  
  YN = YNreplot()

  if(is.numeric(SHOWONLY)) {

    
    Sys.sleep(SHOWONLY);
    return(list(but=NULL, zloc=0, pix=0, YN=YN))

  }
  if(SHOWONLY==TRUE) {

    
    return(list(but=NULL, zloc=0, pix=0, YN=YN))

  }

  global.vars$buttoncex = 0.8
  
  buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs, cex=global.vars$buttoncex )

  global.vars$MAINdev = dev.cur()


###  Get.Screens(2)
  dev.set(global.vars$MAINdev)
 
  
  u = par("usr")
  sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
  zloc =list(x=NULL, y=NULL)
  ppick  = NA
  spick  = NA
  xpick = NA

  zenclick = length(zloc$x)


  global.vars$BLABS = BLABS
  global.vars$buttons = buttons
  global.vars$zloc = zloc
  global.vars$sloc = sloc
  global.vars$zenclick = zenclick

  


   ####  print( data.frame(global.vars$filters) )

  while(TRUE) {
      ####### start while: each mouse click is recorded and tested for what to do next
      #######
      iloc = RPMG::ilocator(global.vars$ilocnum ,COL=global.vars$iloccol ,NUM=FALSE , YN=length(global.vars$sel), style=global.vars$ilocstyle )
      Nclick = length(iloc$x)
####  cat(paste(sep=" ", zenclick, Nclick), sep="\n")
      
      if(Nclick>0)
        {
          #######  add last click to list of clicks, continue 
          zloc  = list(x=c(zloc$x,iloc$x), y=c(zloc$y, iloc$y))
          global.vars$zenclick = length(zloc$x)
          K =  RPMG::whichbutt(iloc ,buttons)
          sloc = zloc
          
          
        }
      else
        {
###  RIGHT button was clicked
          Nclick = 0
###  zenclick=zenclick+1
###   print(zenclick)
          K = 0
          global.vars$zenclick = length(zloc$x)
          if(global.vars$zenclick<1)
            {

           #######  No left mouse click was executed - stop and return to main 
 
              buttons = RPMG::rowBUTTONS(global.vars$BLABS, col=rep(grey(.8), length(global.vars$BLABS)),
                pch=rep("NULL", length(global.vars$BLABS)), cex=global.vars$buttoncex)
              title("Done, Return to Calling Program")
              
              return(list(but="None", zloc=0, pix=0))
            }
          if(global.vars$zenclick==1)
            {
              #############  replot screen
              global.vars$WIN = NULL
              
              YN = YNreplot()
              
              buttons = RPMG::rowBUTTONS(global.vars$BLABS, col=global.vars$colabs, pch=global.vars$pchlabs, cex=global.vars$buttoncex)
              u = par("usr")
              
              sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
              
              Nclick=1
              K = 0
              zloc = list(x=NULL, y=NULL)
              global.vars$zenclick = 0
              
            }
          if(global.vars$zenclick>=2)
            {
              ############  ZOOM
              global.vars$WIN  = sort(zloc$x[c(global.vars$zenclick-1, global.vars$zenclick)])
              YN = YNreplot()
              
              buttons = RPMG::rowBUTTONS(BLABS, col=global.vars$colabs, pch=global.vars$pchlabs, cex=global.vars$buttoncex)
              
              Nclick=1
              K = 0
              zloc = list(x=NULL, y=NULL)
              global.vars$zenclick = 0
            }

        }
      
###  print(paste(sep=' ',  Nclick , global.vars$zenclick) )
############   button actions

#############################################
#############################################
      if(K[Nclick] == match("SELBUT", BLABS, nomatch = NOLAB))
        {

        ##  print('SELBUT' )
          
          labs = BLABS
          
       
          tbuts = unique( c(global.vars$BLABS , ALLlabs) )

          onoff=rep(0 , length(tbuts)) 
          onoff[ match(global.vars$BLABS, tbuts  ) ] = 1
          
          dev.new()
          newlabs = SELBUT(tbuts , onoff=onoff, default=defaultLAB );
          dev.off()
          
           ####     you must have these buttons - these are forced to stay on
          newlabs = unique( c(fixedbuttons, newlabs)) 

          
          cat(as.vector(newlabs), sep=" " );  cat("\n") 
          global.vars$BLABS =newlabs
          
          pnos = c( grep("PIX", global.vars$BLABS), grep("pik", global.vars$BLABS))
          global.vars$colabs = rep(1,length(global.vars$BLABS))
          global.vars$colabs[global.vars$BLABS=="PickWin"] = 'red'
          global.vars$colabs[pnos] = somecolors[seq(from=2, length=length(pnos))]

          
          global.vars$pchlabs = rep(4,length(BLABS))
          global.vars$pchlabs[pnos] = seq(from=15, length=length(pnos))

          
          BLABS=global.vars$BLABS
          
          YN = YNreplot()
          
          buttons = RPMG::rowBUTTONS(global.vars$BLABS, col=global.vars$colabs, pch=global.vars$pchlabs, cex=global.vars$buttoncex)
          Nclick = 0
          K = 0
          zloc = list(x=NULL, y=NULL)
          next
          
        }

      if(K[Nclick] == match("REPLOT", BLABS, nomatch = NOLAB))
        {
          YN = YNreplot()
          ## print("clicked REPLOT")
          buttons = RPMG::rowBUTTONS(global.vars$BLABS, col=global.vars$colabs, pch=global.vars$pchlabs, cex=global.vars$buttoncex)
          Nclick = 0
          K = 0
          zloc = list(x=NULL, y=NULL)
        
          next;
        }
      

      
      if(K[Nclick] == match("DONE", BLABS, nomatch = NOLAB))
        {
          buttons = RPMG::rowBUTTONS(BLABS, col=rep(grey(.8), length(BLABS)), pch=rep("NULL", length(BLABS)), cex=global.vars$buttoncex)
          title("Return to Calling Program")

          if(global.vars$zenclick>1)
            {
              rd = getrdpix(zloc, global.vars$zenclick, sel, NH)
            }
          else
            {
              rd=list(PUSHED="DONE")
            }


          
          
          invisible(list(but=BLABS[K[Nclick]], zloc=zloc, pix=rd, g=global.vars ))
          break;
        }

      if(K[Nclick] == match("QUIT", BLABS, nomatch = NOLAB))
        {
          buttons = RPMG::rowBUTTONS(BLABS, col=rep(grey(.8), length(BLABS)), pch=rep("NULL", length(BLABS)), cex=global.vars$buttoncex)
          title("Return to Calling Program")


          if(global.vars$zenclick>1)
            {
              rd = getrdpix(zloc, global.vars$zenclick, sel, NH)
            }
          else
            {
              rd=list(PUSHED="QUIT")
            }
          
          
          invisible(list(but=BLABS[K[Nclick]], zloc=zloc, pix=rd, g=global.vars))
          break;
        }

       ###################   Postscript output   ###########################  
      if(K[Nclick] == match("Postscript", BLABS, nomatch = NOLAB))
        {
          print(paste(sep=' ' ,"Start swig postscript"))
          jdev = dev.cur()
          plfname = RPMG::local.file("swig","eps")
          
          ### postscript(file=plfname, horizontal=TRUE, print.it=FALSE,  onefile=FALSE)
           P = round(par('pin'))

          opar = par(no.readonly = TRUE) 

           postscript(file=plfname , width=P[1], height=P[2], paper = "special", bg=opar$bg, fg=opar$fg, horizontal=FALSE, onefile=TRUE,print.it=FALSE)

          ### par(OPAR)
          print(paste(sep=' ' ,"Doing postscript", plfname))
           YN = YNreplot()
          
          print(paste(sep=' ' ,"Done creating postscript file: ", plfname))
          dev.off()
          dev.set(jdev)
             zloc = list(x=NULL, y=NULL)
          next
        }

      

############################################################################3
############################################################################3


      
      if(K[Nclick]>0)
        {
          if(K[Nclick] != match(BLABS[K[Nclick]]  ,  fixedbuttons , nomatch = NOLAB))
            {
              ###print("You Hoo!")
              ###print(BLABS[K[Nclick]])
              if(exists(BLABS[K[Nclick]],  mode="function"))
                {
                  FUN = match.fun( BLABS[K[Nclick]] )
                  if(is.function(FUN))
                    {
                      global.vars$zloc = zloc
                      global.vars$KLICK = BLABS[K[Nclick]]
                      GL =  FUN(NH, global.vars)
                      if(!is.null(GL))
                        {
                          NGL = names(GL)
                          if(length(grep("global.vars",NGL)))  global.vars = GL$global.vars
                        }

                      
                    }
                  else
                    {
                      print(paste("No Function named", BLABS[K[Nclick]] , "Found"))
                      global.vars$action=="donothing"
                      zloc = list(x=NULL, y=NULL) 

                    }
                }
              else
                {

                  print(paste("No Function named", BLABS[K[Nclick]] , "Found"))
                  global.vars$action=="donothing"
                  zloc = list(x=NULL, y=NULL) 
                  

                }
              
              if(global.vars$action=="break")
                {
                  invisible(list(but=global.vars$BLABS[K[Nclick]], zloc=zloc,
                                 pix=global.vars$rd))
                  break
                }
              
              if(global.vars$action=="replot")
                {
                  
                  YN = YNreplot()
                  
                  buttons = RPMG::rowBUTTONS(global.vars$BLABS, col=global.vars$colabs,
                    pch=global.vars$pchlabs, cex=global.vars$buttoncex)
                  zloc = global.vars$zloc
                }
              if(global.vars$action=="donothing")
                {
                   zloc = global.vars$zloc
                }
              
              if(global.vars$action=="replace")
                {
                  if(!is.null(GL))
                    {
                       NGL = names(GL)
                       if(length(grep("NH",NGL)))     NH = GL$NH
                    }
                  YN = YNreplot()   
                  buttons = RPMG::rowBUTTONS(global.vars$BLABS, col=global.vars$colabs,
                    pch=global.vars$pchlabs, cex=global.vars$buttoncex)
                  zloc = global.vars$zloc
                }
              
              if(global.vars$action=="revert")
                {
                  if(exists("OLDH"))NH = OLDH
                 
                  YN = YNreplot()   
                  buttons = RPMG::rowBUTTONS(global.vars$BLABS, col=global.vars$colabs,
                    pch=global.vars$pchlabs, cex=global.vars$buttoncex)
                  zloc = global.vars$zloc
                }

              
              if(global.vars$action=="exit")
                {
                  print("EXIT RETX")
                  return(GL$RETX)
                 break
                }

              next

            }
        }


############################################################################3
############################################################################3
      
      

    }

### PRET = list(TPIX=TPIX, xpix=xpix,ypixA=ypixA, ypixB=ypixB,cpixa= cpixa, cpixb=cpixb, cpixc=cpixc, colpix=colpix)
###  return(PRET)
  if(global.vars$zenclick>2)
    {
      pwin = sort(c(zloc$x[global.vars$zenclick-2], zloc$x[global.vars$zenclick-1]))
    }
  else
    {
      pwin =  c( zloc$x[global.vars$zenclick-2], zloc$x[global.vars$zenclick-1])

    }
#### print(pwin)
#### print(WIN)

  
  PushI =  RPMG::whichbutt(zloc ,buttons)
####  print(zloc)
####  print(sloc)
####  print(PushI)
  
  
  PushK=NULL
  if(length(PushI)>=1)
    {
      PushK=NULL
      if(any(PushI>0)) PushK =BLABS[PushI[PushI>0]]
    }
  
  
  but=BLABS[K[Nclick]]

  WPX = global.vars$WPX



  if(any(is.na(WPX$name)))
    {
      idpx = which(is.na(WPX$name))
      
      WPX =  deleteWPX(WPX, idpx)
    }



  
  ##############  clean off undesirable picks
 ##### whirid = which( WPX$name==NA & WPX$comp==NA & WPX$phase==NA )

  
#####  WPX = deleteWPX(WPX, whirid)

#####  whirid = which( WPX$name=='99999' & WPX$comp=='9' & WPX$phase=="X" )

#####  WPX = deleteWPX(WPX, whirid)

  

 #####  print("returning from swig")
  
  RETP = list(but=but, sloc=sloc, WPX=WPX, BRUNINFO=BRUNINFO, DETLINFO=DETLINFO,  mark=mark, PUSHED=PushK, g=global.vars)
  
  

      invisible(RETP)
   
}

