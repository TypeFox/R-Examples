`SNET.drive` <-
function(intempmat, pmolabs=c("Vertical", "North", "East"), STAMP="")
  {

    if(missing(STAMP)) { STAMP = " " }
    if(missing(pmolabs)) {pmolabs=c("Vertical", "North", "East")  }
    
     TPALS = c( "tomo.colors", "JBLACK", "rainbow", "topo.colors", "terrain.colors", "JGRAY")
     APALS = c( "tomo", "JBLACK", "rainbow", "topo", "terrain", "JGRAY")

    #########  these are the default colors...
   ### TPALS = c(  "rainbow", "topo.colors", "terrain.colors")
   ### APALS = c(  "rainbow", "topo", "terrain")

    
    ADDBUTS = c("More" )
  

    
    labs = c("DONE",  "LINES", "Postscript", APALS, ADDBUTS )
    NLABS = length(labs)
    NOLAB = NLABS +1000
###  FUN = match.fun(TPALS[1])
     pal = RPMG::Gcols(plow=0, phi=0,  N=100, pal=TPALS[1])
    scale.def = 0
   
    colabs = c(rep(1,2) , rep(2, length(APALS) ), rep(4,length(ADDBUTS) ))
    pchlabs = c(rep(1,2) , rep(2, length(APALS) ), rep(4,length(ADDBUTS) ))
 
    gridon = FALSE
    ADDLINES = FALSE
  
    NSEL = 1

  ###  X11()
###  

    temp = list(x=intempmat[,1], y=intempmat[,2], z=intempmat[,3])

    
    sx = partmotnet(temp, STAMP=STAMP, LINES=ADDLINES, COL=pal)

    buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)





        zloc = list(x=NULL, y=NULL)
    sloc = zloc

    
    while(TRUE)
      {


iloc = RPMG::ilocator(1, COL=rgb(1,0.8, 0.8), NUM=FALSE , YN=1, style=-1)
 Nclick = length(iloc$x)
           
          if(Nclick>0)
            {
              zloc  = list(x=c(zloc$x,iloc$x), y=c(zloc$y, iloc$y))
              zenclick = length(zloc$x)
              K =  RPMG::whichbutt(iloc ,buttons)
              sloc = zloc
            }
          else
            {
              Nclick = 0
              K = 0
              buttons = RPMG::rowBUTTONS(labs, col=rep(grey(.8), length(labs)), pch=rep("NULL", length(labs)))
              title("Return to Calling Program")
              break;
            }
     
       


        
        if(K[Nclick] == match("DONE", labs, nomatch = NOLAB))
          {
           
            buttons = RPMG::rowBUTTONS(labs, col=rep(grey(.8), length(labs)), pch=rep("NULL", length(labs)))
            title("Return to Calling Program")
        
            
            break;
          }
        if(zenclick == 1 &  K[Nclick]==0 )
          {
            ###  replot
            ###print(K[Nclick])
            
           sx = partmotnet(temp, STAMP=STAMP, LINES=ADDLINES, COL=pal)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
           zloc = list(x=NULL, y=NULL) 
           
          }

      ####################   POSTSCRIPT  ##################
      
        if(K[Nclick] == match("Postscript", labs, nomatch = NOLAB))
        {

          print("Start postscript plot.ts")
          plfname = RPMG::local.file("pmotnet","eps")
          jdev = dev.cur()
          RPMG::jpostscript("pmot")
          sx = partmotnet(temp, STAMP=STAMP, LINES=ADDLINES, COL=pal)
           print("Done creating postscript")
          dev.off()
          dev.set(jdev)
          zloc = list(x=NULL, y=NULL) 
        }
        
        if(K[Nclick]==match("LINES", labs, nomatch = NOLAB))
          {
            ADDLINES=!ADDLINES
            sx = partmotnet(temp, STAMP=STAMP, LINES=ADDLINES, COL=pal)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
            zloc = list(x=NULL, y=NULL) 
          }


              #######  source("drivers.R"); save.image()

        
        if(K[Nclick]==match("More", labs, nomatch = NOLAB))
          {
           
            sx = partmotnet(temp, STAMP=STAMP, LINES=ADDLINES, COL=pal)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
          }
        if( length(which(K[Nclick] == match(APALS, labs, nomatch = NOLAB)))>0 )
          {
            J = match(labs[K[Nclick]] ,  APALS   )
            
            ##FUN = match.fun(TPALS[J])
            ##  pal = FUN(NCOL)
            pal = RPMG::Gcols(plow=0, phi=0,  N=100, pal=TPALS[J])
            
            sx = partmotnet(temp, STAMP=STAMP, LINES=ADDLINES, COL=pal)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)

              zloc = list(x=NULL, y=NULL) 
          }

       
         
         


      }

    print("DONE with partmotnet")
    
    
  }

