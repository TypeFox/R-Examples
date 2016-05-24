`SPECT.drive` <-
function(Xamp, DT=0.008, NEW=TRUE, STAMP=NULL, freqlim=c(0, 20, 0, 20), winparams=c(4096,256, 204 ) )
  {

    if(missing(DT)) { DT=1 }
    if(missing(NEW)) { NEW=TRUE }
    if(missing(STAMP)) { STAMP=NULL }

    
    TPALS = c("rainbow", "topo.colors", "terrain.colors", "heat.colors", "tomo.colors")
    APALS = c("rainbow", "topo",           "terrain",      "heat",       "tomo")
    ADDBUTS = c("RAW", "LOG", "SQRT", "Grid", "INFO", "DOT",  "WLEN", "PARMS" , "SAVE" )
    
    NCOL = 100
    
    labs = c("DONE", "Postscript", APALS, ADDBUTS )
    NLABS = length(labs)
    NOLAB = NLABS +1000
  ###  FUN = match.fun(TPALS[1])
     
   ###  pal = FUN(NCOL)

    pal = RPMG::Gcols(plow=5, phi=0,  N=100, pal=TPALS[1])
    scale.def = 0
   
    colabs = c(rep(1,2) , rep(2, length(APALS) ), rep(4,length(ADDBUTS) ))
    pchlabs = c(rep(1,2) , rep(2, length(APALS) ), rep(4,length(ADDBUTS) ))
 
    gridon = FALSE
 
    NSEL = 1

    if(NEW==TRUE)
      {
       ### get(getOption("device"))(width=15, height=10)
        dev.new(width=15, height=10)
       ### X11(width=15, height=10)
      }
###  
###  fh is one half  the nyquist  frequency
    if(missing(freqlim))
      {
        fl=0
        fh=0.25*(1/DT)
        flshow = 0
###  fh is 60% the fh  frequency
        fhshow = round(0.6*fh)
      }
    else
      {
        fl=freqlim[1]
        fh=freqlim[2]
        flshow = freqlim[3]
        fhshow = freqlim[4]
      }

   ###  multi = 1

    ###  here need to choose a length for the roving window.
    ###  I used to use 2 second windows by default, but here I think we need to
    ###  be more creative - how about 2% of the total length of the record?
##  tsecs = DT*(length(Xamp)*.02)
    ### No, now lets do 256 per window, approximately


    ###  need 3 parameters:  FFT length
#################   time window length and time Overlap


    if(missing(winparams))
      {
        
        TWOSEC = 256
        NS = floor(TWOSEC)
        NOV = floor(TWOSEC-.2*TWOSEC)

        
        Nfft=4096

      }
    else
      {
        TWOSEC = winparams[2]
        NS = floor(TWOSEC)
        NOV = floor(winparams[3])
        Nfft=floor(winparams[1])
      }

    
    DEV = evolfft(Xamp, DT , Nfft=Nfft, Ns=NS , Nov=NOV,  fl=fl, fh=fh  )

    if(is.null(DEV)) {
       emsg = c("ERROR: illegal call to evolfft.", "Exit Program SPECT.drive")
        cat(emsg, sep="\n")
        return(NULL)
    }
       
    
    PE = plotevol(DEV, log=scale.def, fl=flshow, fh=fhshow, col=pal, ygrid=gridon, STAMP=STAMP)
    
    buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)

    pstyle=3
   
     zloc = list(x=NULL, y=NULL)
    sloc = zloc

  ####   print(paste("Button",K, Nclick,K[Nclick] ))





    
    while(TRUE)
      {


    iloc = RPMG::ilocator(1, COL=rgb(1,0.8, 0.8), NUM=FALSE , YN=1, style=0)
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

        if(K[Nclick] == match("SAVE", labs, nomatch = NOLAB))
          {
            buttons = RPMG::rowBUTTONS(labs, col=rep(grey(.8), length(labs)), pch=rep("NULL", length(labs)))
            title("Return to Calling Program")
            invisible(DEV)
            break;
          }



    
      
      ####################   POSTSCRIPT  ##################
        if(K[Nclick] == match("Postscript", labs, nomatch = NOLAB))
        {

          print(paste(sep=' ' ,"Starting postscript file plotevol"))
           jdev = dev.cur()
           plfname = RPMG::local.file("sgram","eps")
           P = round(par('pin'))

           postscript(file=plfname , width=P[1], height=P[2], paper = "special", horizontal=FALSE, onefile=TRUE,print.it=FALSE)

          ## postscript(file=plfname, horizontal=TRUE, print.it=FALSE,  onefile=FALSE)
           PE = plotevol(DEV, log=scale.def, fl= flshow, fh=fhshow, col=pal, ygrid=gridon, STAMP=STAMP)
           cat(paste(sep=' ' ,"Done creating postscript file: ", plfname), sep="\n")

          dev.off()
          dev.set(jdev)
          zloc = list(x=NULL, y=NULL) 
        }

       
        if( length(which(K[Nclick] == match(APALS, labs, nomatch = NOLAB)))>0 )
          {
            J = match(labs[K[Nclick]] ,  APALS   )
            
            ##FUN = match.fun(TPALS[J])
            ##  pal = FUN(NCOL)
            pal = RPMG::Gcols(plow=5, phi=0,  N=100, pal=TPALS[J])
            
              PE = plotevol(DEV, log=scale.def, fl= flshow, fh=fhshow, col=pal, ygrid=gridon, STAMP=STAMP)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
            zloc = list(x=NULL, y=NULL) 
          }


###     "RAW", "LOG", "SQRT"

        if(K[Nclick]==match("RAW", labs, nomatch = NOLAB))
          {
            scale.def = 0
         
              PE = plotevol(DEV, log=scale.def, fl= flshow, fh=fhshow, col=pal, ygrid=gridon, STAMP=STAMP)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
             zloc = list(x=NULL, y=NULL) 
        }
        if(K[Nclick]==match("LOG", labs, nomatch = NOLAB))
          {
            scale.def = 1
           PE =  plotevol(DEV, log=scale.def, fl= flshow, fh=fhshow, col=pal, ygrid=gridon, STAMP=STAMP)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
             zloc = list(x=NULL, y=NULL) 
        }
        if(K[Nclick]==match("SQRT", labs, nomatch = NOLAB))
          {
            scale.def = 2
            PE = plotevol(DEV, log=scale.def, fl=flshow, fh=fhshow, col=pal, ygrid=gridon, STAMP=STAMP)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
             zloc = list(x=NULL, y=NULL) 
        }

         if(K[Nclick]==match("Grid", labs, nomatch = NOLAB))
          {
            
            gridon = !gridon
            
            
           PE =  plotevol(DEV, log=scale.def, fl= flshow, fh=fhshow, col=pal, ygrid=gridon, STAMP=STAMP)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
             zloc = list(x=NULL, y=NULL) 
        }


         if(K[Nclick]==match("WLEN", labs, nomatch = NOLAB))
           {

             A = readline(prompt ="Type in the window multiplier:")
             multi=as.numeric(A)
            
             
             NS = multi*TWOSEC
             NOV = multi*(TWOSEC-10)
             
             
             DEV = evolfft(Xamp,DT , Nfft=Nfft, Ns=NS , Nov=NOV,  fl=fl, fh=fh  )
  
            PE = plotevol(DEV, log=scale.def, fl=flshow, fh=fhshow, col=pal, ygrid=gridon, STAMP=STAMP)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
              zloc = list(x=NULL, y=NULL) 
        }

         if(K[Nclick]==match("PARMS", labs, nomatch = NOLAB))
           {


             prom = paste(sep=" ","Current 7 Parameters: Nfft Ns Nov  fl fh flshow fhshow\n=",Nfft, NS, NOV , fl, fh, flshow, fhshow, "\nINPUT> " )  
             A = readline(prompt =prom)

             
             
             RL= as.numeric(unlist(strsplit(split=" ", A)))
             Nfft = RL[1]
             NS = RL[2]
        
             NOV = RL[3]
             fl  = RL[4]
             fh = RL[5]
             flshow = RL[6]
             fhshow = RL[7]

             
             DEV = evolfft(Xamp,DT , Nfft=Nfft, Ns=NS , Nov=NOV,  fl=fl, fh=fh  )
             if(!is.null(DEV))
               {
                 PE =  plotevol(DEV, log=scale.def, fl=flshow, fh=fhshow, col=pal, ygrid=gridon, STAMP=STAMP)
               }
             buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
              zloc = list(x=NULL, y=NULL) 
        }

        if(K[Nclick]==match("INFO", labs, nomatch = NOLAB))
          {

            ###   print(zloc)
            
            whyat = min(PE$y)+(diff(range(PE$y)))*(zloc$y[1:(Nclick-1)]-min(PE$why))/(diff(range(PE$why)))
            exat  = min(PE$x)+(diff(range(PE$x)))*(zloc$x[1:(Nclick-1)]-min(PE$x))/(diff(range(PE$x)))
            rDUMPLOC(list(x=exat, y=whyat))
             zloc = list(x=NULL, y=NULL) 
            
           ###  print(exat)
           ###  print(whyat)            
          }

        if(K[Nclick]==match("DOT", labs, nomatch = NOLAB))
          {
            if(pstyle!=(-1))
              {
                pstyle = -1
              }
            else
              {
                pstyle=3
              }

             zloc = list(x=NULL, y=NULL)
            
          }


        
      
      
        ###  print(paste(sep=' ', "scale.def=", scale.def))
     

      }

    cat("DONE with Sgram", sep="\n")
    
    
  }

