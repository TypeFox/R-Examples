XSECDEMg<-function(Data, labs=NULL, pts=NULL,  nlevels=10, demo=FALSE)
  {
    if(missing(demo)) { demo = FALSE }
#######   this program illustrates how to use the interactive program RPMG
    
#####  to make the rainbow color palette named like other palettes in R
    
  ###  this needs to be defined
    #####  rainbow.colors<-function(n){ return( rainbow(n) ) }

    if(missing(labs)) {
      labs = c("DONE", "REFRESH","CONT","PTS", 
     "XSEC","PS" )
    }


    if(is.null(labs))
      {
      labs = c("DONE", "REFRESH","CONT","PTS", 
     "XSEC","PS" )
      }
    if(missing(nlevels)) { nlevels = 10  }

    nx = dim(Data)[1]
    ny = dim(Data)[2]
    
    jx = attr(Data, 'dx')*seq(from=0, to=nx-1)
    jy = attr(Data, 'dy')*seq(from=0, to=ny-1)

    
    iseclab  = 0
    secmat = NULL
    ncol = 100
    GENsec = vector(mode="list")
    mycol.colors<-function(n)
      {
        n2 = floor(n/2)
        return(c(topo.colors(n2), terrain.colors(n2) ))   
      }

    
 TPALS =    c("rainbow", "topo.colors", "terrain.colors", "heat.colors", "tomo.colors", "mycol.colors")
    APALS = c("rainbow", "topo",         "terrain",        "heat",       "tomo",        "mycol" )
    labs = c(labs,APALS)

    
  ##  TPALS = c("rainbow", "topo", "terrain", "tomo",  "mycol")
    colabs = rep(1, length=length(labs))
    pchlabs = rep(0,length(labs))
 ##    FUN = match.fun(TPALS[1])
##     pal = FUN(ncol)
    pal = RPMG::Gcols(plow=0, phi=0,  N=100, pal=TPALS[1])
    
    image(jx, jy, Data, col=pal , asp=1)
    buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
    NLABS = length(labs)
    NOLAB = NLABS +1000  ## some large number
    
    if(demo==TRUE)  return(NULL)
    iloc = locator(1, type='p')
    zloc = iloc

    Nclick = length(iloc$x)
    if(is.null(zloc$x)) { return(NULL) }
    K =  RPMG::whichbutt(zloc , buttons)
    sloc = zloc
    CONT.FLAG = FALSE
    PTS.FLAG = FALSE
    
    XSEC.FLAG = FALSE
    PS.FLAG =  FALSE

    NSEC = 0
    
    
    while(TRUE)
      {
        ############   button actions

        ###########   quit and break loop
        if(K[Nclick] == match("DONE", labs, nomatch = NOLAB))
          {
            break;
          }

        ###########   refresh the screen
        if(K[Nclick] == match("REFRESH", labs, nomatch = NOLAB))
          {
            zloc = list(x=NULL, y=NULL)
          }

        
        if( length(which(K[Nclick] == match(APALS, labs, nomatch = NOLAB)))>0 )
          {
            J = match(labs[K[Nclick]] ,  APALS   )
            
            ##  FUN = match.fun(TPALS[J])
             ##  pal = FUN(NCOL)
            
            pal = RPMG::Gcols(plow=0, phi=0,  N=100, pal=TPALS[J])
          
        
            zloc = list(x=NULL, y=NULL)
          }

        ################  add contours
        if(K[Nclick] == match("CONT", labs, nomatch = NOLAB))
          {
            CONT.FLAG = !CONT.FLAG
            zloc = list(x=NULL, y=NULL)
          }
        if(K[Nclick] == match("PTS", labs, nomatch = NOLAB))
          {
            PTS.FLAG = !PTS.FLAG
            zloc = list(x=NULL, y=NULL)
          }

        ################  make postscript file
        if(K[Nclick] == match("PS", labs, nomatch = NOLAB))
          {
            PS.FLAG = !PS.FLAG
            zloc = list(x=NULL, y=NULL)
          }

################  cut cross sections
        if(K[Nclick] == match("XSEC", labs, nomatch = NOLAB))
          {
            n = length(zloc$x)
            if(n<3) {
              print("Need 2 locations on the plot to take the cross section:")
              print("Please try again")
              zloc = list(x=NULL, y=NULL)
              next
            }
            x1 = zloc$x[n-2]
            y1 = zloc$y[n-2]
            x2 = zloc$x[n-1]
            y2 = zloc$y[n-1]

            XSEC.FLAG = TRUE
            print(c(x1,y1, x2,y2))
            iseclab = iseclab + 1
            LAB = LETTERS[iseclab]
            secmat = rbind(secmat, c(x1, y1, x2,  y2))
          GENsec[[iseclab]] =  GETXprofile(jx, jy, Data, myloc=list(x=c(x1, x2), y=c(y1, y2)), LAB=LAB, PLOT=TRUE)
            
            zloc = list(x=NULL, y=NULL)
          }

        if(K[Nclick] > 0)
          {
            if(PS.FLAG) {
              P = round(par('pin'), digits=2); 
              postscript(file="RPMGdemo.eps"  , width=P[1], height=P[2],
                         paper = "special", horizontal=FALSE, onefile=TRUE,print.it=FALSE)
            }
            image(jx, jy, Data, col=pal , asp=1)
            if(CONT.FLAG) contour(x=jx, y=jy, Data, nlevels=nlevels ,  add=TRUE)
            if(PTS.FLAG) points(pts$x, pts$y, pch=25, bg='gold', fg='blue' )
            
            if(XSEC.FLAG) 
              {
                segments(secmat[,1],secmat[,2],secmat[,3],secmat[,4])
                seclabs = LETTERS[1:iseclab]
                text(secmat[,1],secmat[,2], labels= seclabs, pos=3)
                text(secmat[,3],secmat[,4], labels= paste(sep="",seclabs, "'") , pos=3)


              }
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
            if(PS.FLAG) {
              dev.off();
              cat("the postscript file is in: RPMGdemo.eps", sep="\n")
              PS.FLAG =  FALSE
            }


          }
        else
          {
###  in case the plot was resized with asp=1, need to replot the buttons
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
          }


        iloc = locator(1,type='p')
##### print(iloc)
        zloc  = list(x=c(zloc$x,iloc$x), y=c(zloc$y, iloc$y))
        Nclick = length(iloc$x)
        if(is.null(zloc$x)) { return(sloc) }
        K =  RPMG::whichbutt(iloc , buttons)
##### print(K)   
      }

    invisible(GENsec)
    
  }

