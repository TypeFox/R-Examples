`wlet.drive` <-
function(Xamp, DT=0.008, STAMP=NULL)
  {

      
 
    if(missing(DT)) { DT=1 }
     if(missing(STAMP)) { STAMP=NULL }


    
   TPALS = c("rainbow", "topo.colors", "terrain.colors", "heat.colors", "tomo.colors")
    APALS = c("rainbow", "topo", "terrain", "heat", "tomo")
    ADDBUTS = c("RAW", "LOG", "SQRT", "INFO","BOXZ", "CLIMB", "ZBOUNDS" )
   NCOL = 100
    labs = c("DONE", "Postscript", "ContPS",  APALS,  ADDBUTS)
    NLABS = length(labs)
    NOLAB = NLABS +1000

    pfreqs = c(0.5, 1, 2,3,4,5, 10, 14)
   ###    FUN = match.fun(TPALS[1])

   ###    pal = FUN(NCOL)

    pal = RPMG::Gcols(plow=5, phi=0,  N=100, pal=TPALS[1])
    
    scale.def = 3
    colabs = c(rep(1,2) , rep(2, length(APALS) ), rep(4,length(ADDBUTS) ))
    pchlabs = c(rep(1,2) , rep(2, length(APALS) ), rep(4,length(ADDBUTS) ))
 
    NSEL = 1

    ### get(getOption("device"))(width=15, height=10)
       dev.new(width=15, height=10)
       
   ###  X11(width=15, height=10)
###  
    WOUT =  wlet.do(Xamp, DT, noctave=7, zscale=scale.def,  col=pal, STAMP=STAMP)
    print("Finished calculation for  wlet.do")

    faha = WOUT$baha
    PE = WOUT$PE
    pwlet2freqs(faha$noctave, faha$nvoice, DT, flip=faha$flip, pfreqs , plot=TRUE)

    
    buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
    
   
    PEAXBoxes = NULL
 ####    iloc = RPMG::ilocator(1, COL=rgb(1,0.8, 0.8), NUM=FALSE , YN=NSEL, style=1)
    
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
           #### PE =plotwlet(faha, Xamp, DT, zscale=scale.def,  col=pal, STAMP=STAMP)
           #### pwlet2freqs(faha$noctave, faha$nvoice, DT, flip=faha$flip,  pfreqs , plot=TRUE)

            buttons = RPMG::rowBUTTONS(labs, col=rep(grey(.8), length(labs)), pch=rep("NULL", length(labs)))
            title("Return to Calling Program")
            break;
          }
      
        if(K[Nclick] == match("ZBOUNDS", labs, nomatch = NOLAB))
          {
            print("Click pairs of points in window, end with middle mouse")

            zloc = plocator(COL=rgb(1,0.8, 0.8), NUM=FALSE , YN=NSEL, style=-1)
            print(zloc)
            numc = ncol(faha$img)
            numf = nrow(faha$img)
            
            szi = seq(from=1, to=(length(zloc$x)-1), by=2)
            
            MYx = rep(NA, length(szi))
            MYy = rep(NA, length(szi))
            j = 1

            
            zi = length(szi)
              
                
                rect(zloc$x[zi], zloc$y[zi],zloc$x[zi+1],zloc$y[zi+1], border=2)
                
                why1 = round(numc*(zloc$y[zi]-min(PE$why))/(diff(range(PE$why))))
                ex1  = round(numf*(zloc$x[zi]-min(PE$x))/(diff(range(PE$x))))
                
                why2 = round(numc*(zloc$y[zi+1]-min(PE$why))/(diff(range(PE$why))))
                ex2  = round(numf*(zloc$x[zi+1]-min(PE$x))/(diff(range(PE$x))))
                
                
                MAT = faha$img[ex1:ex2, why1:why2  ]
                w = which.max(MAT)
            wmin = which.min(MAT)
                c1 = col(MAT)[w]+why1
                r1 = row(MAT)[w]+ex1

                PE$x[r1]
                PE$y[c1]
                points(PE$x[r1], PE$why[c1])
                abline(v=PE$x[r1])

                
                MYx[j] = PE$x[r1]
                MYy[j] = PE$why[c1]
                j = j+1
                
              
            rDUMPLOC(list(x=MYx, y=MYy))
            
            print(paste(sep=" ", "MAX=", MAT[w], "MIN=", MAT[wmin]))
            
            dbounds = range(MAT)
            print(dbounds)

            
         ###   if(FALSE) {   ### this code is removed
         ###   RIM = range(faha$img, na.rm=TRUE)

            
        ###    print("Type in the upper and lower Z Bounds for the plot")
         ###   
         ###   JEAD = readline(prompt="ZBounds: ")
            
         ###   dbounds = as.numeric(unlist(strsplit(JEAD, split=" ")))
         ###   if(length(dbounds)<2)dbounds=NULL
        ###  }

            zloc = list(x=NULL, y=NULL) 
            
            
          }
      
          

      ####################   POSTSCRIPT  ##################
      
        if(K[Nclick] == match("Postscript", labs, nomatch = NOLAB))
          {
            print(paste(sep=' ' ,"Start postscript file plotwlet"))
            jdev = dev.cur()
            plfname = RPMG::local.file("wlet","eps")


            P = round(par('pin'))

            postscript(file=plfname , width=P[1], height=P[2], paper = "special", horizontal=FALSE, onefile=TRUE,print.it=FALSE)


            
           ##  postscript(file=plfname, horizontal=TRUE, print.it=FALSE,  onefile=FALSE)

            
            PE =plotwlet(faha, Xamp, DT, zscale=scale.def, col=pal,  ygrid=FALSE, STAMP=STAMP)
            pwlet2freqs(faha$noctave, faha$nvoice, DT, flip=faha$flip, pfreqs  , plot=TRUE)
            
            print(paste(sep=' ' ,"Done creating postscript file:", plfname))

            
            dev.off()
            dev.set(jdev)
            zloc = list(x=NULL, y=NULL) 
            
          }

        if(K[Nclick] == match("ContPS", labs, nomatch = NOLAB))
          {
            print(paste(sep=' ' ,"Start postscript file plotwlet"))
            jdev = dev.cur()
            plfname = RPMG::local.file("wletc","eps")

            
            dbounds = range(faha$img )

            P = round(par('pin'))

            postscript(file=plfname , width=P[1], height=P[2], paper = "special",
                       horizontal=FALSE, onefile=TRUE,print.it=FALSE)


            
           ##  postscript(file=plfname, horizontal=TRUE, print.it=FALSE,  onefile=FALSE)

            if(is.null(PEAXBoxes) )
              {
                PEAX=NULL
              }
            else
              {
                PEAX=PEAXBoxes

                 print("bef cont: PEAX:")
                print(PEAX)
              }
            PE =contwlet(faha, Xamp, DT,   clev=0.75,  NLEV=20, zscale=scale.def, zbound=dbounds,
              col=col, ygrid=FALSE, PEAX=PEAX)
            pwlet2freqs(faha$noctave, faha$nvoice, DT, flip=faha$flip, pfreqs  , plot=TRUE)
            
            print(paste(sep=' ' ,"Done creating postscript file:", plfname))

            
            dev.off()
            dev.set(jdev)

             zloc = list(x=NULL, y=NULL) 
          }

        
       
        if( length(which(K[Nclick] == match(APALS, labs, nomatch = NOLAB)))>0 )
          {
            J = match(labs[K[Nclick]] ,  APALS   )
            
            ##  FUN = match.fun(TPALS[J])
             ##  pal = FUN(NCOL)
            
            pal = RPMG::Gcols(plow=5, phi=0,  N=100, pal=TPALS[J])
 
           PE = plotwlet(faha, Xamp, DT, zscale=scale.def, col=pal,  ygrid=FALSE, STAMP=STAMP )
            pwlet2freqs(faha$noctave, faha$nvoice, DT, flip=faha$flip, pfreqs , plot=TRUE)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
             zloc = list(x=NULL, y=NULL) 
          }

        if(K[Nclick]==match("RAW", labs, nomatch = NOLAB))
          {
            scale.def = 1
            PE =plotwlet(faha, Xamp, DT, zscale=scale.def,  col=pal, STAMP=STAMP)
            pwlet2freqs(faha$noctave, faha$nvoice, DT, flip=faha$flip,  pfreqs , plot=TRUE)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
             zloc = list(x=NULL, y=NULL) 
          }
        if(K[Nclick]==match("LOG", labs, nomatch = NOLAB))
          {
            scale.def = 2
            PE =plotwlet(faha, Xamp, DT, zscale=scale.def,  col=pal, STAMP=STAMP)
            pwlet2freqs(faha$noctave, faha$nvoice, DT, flip=faha$flip,  pfreqs , plot=TRUE)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
             zloc = list(x=NULL, y=NULL) 
          }
        if(K[Nclick]==match("SQRT", labs, nomatch = NOLAB))
          {
            scale.def = 3
            PE =plotwlet(faha, Xamp, DT, zscale=scale.def,  col=pal, STAMP=STAMP)
            pwlet2freqs(faha$noctave, faha$nvoice, DT, flip=faha$flip,  pfreqs , plot=TRUE)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
             zloc = list(x=NULL, y=NULL) 
          }
        
       if(K[Nclick]==match("CLIMB", labs, nomatch = NOLAB))
          {
            PE = plotwlet(faha, Xamp, DT, zscale=scale.def,  col=pal, STAMP=STAMP)
            pwlet2freqs(faha$noctave, faha$nvoice, DT, flip=faha$flip,  pfreqs , plot=TRUE)
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)

            clHOW = crc(faha$img, nbclimb=1000)
            cfHOW = cfamily(clHOW)
            
            SH  = cfHOW$ordered
            
            SH[SH==0] = NA
            
            image(SH, add=TRUE, col=tomo.colors(100))
            
            zloc = list(x=NULL, y=NULL) 
            
          }
        
        if(K[Nclick]==match("INFO", labs, nomatch = NOLAB))
          {

            points(zloc$x[1:(zenclick-1)], zloc$y[1:(zenclick-1)], pch =3)
            
            whyat = min(PE$y)+(diff(range(PE$y)))*(zloc$y[1:(zenclick-1)]-min(PE$why))/(diff(range(PE$why)))
            exat  = min(PE$x)+(diff(range(PE$x)))*(zloc$x[1:(zenclick-1)]-min(PE$x))/(diff(range(PE$x)))

            
            CLICKEDpts = list(x=exat, y=whyat)
            rDUMPLOC(CLICKEDpts)
             zloc = list(x=NULL, y=NULL) 
          }

        if(K[Nclick]==match("BOXZ", labs, nomatch = NOLAB))
          {  
            zloc = plocator(COL=rgb(1,0.8, 0.8), NUM=FALSE , YN=NSEL, style=-1)
            print(zloc)
            numc = ncol(faha$img)
            numf = nrow(faha$img)

            
            szi = seq(from=1, to=(length(zloc$x)-1), by=2)
            MYx = rep(NA, length(szi))
            MYy = rep(NA, length(szi))
            j = 1
            for(zi in szi)
              {
                
                rect(zloc$x[zi], zloc$y[zi],zloc$x[zi+1],zloc$y[zi+1], border=2)

              
                
                why1 = round(numc*(zloc$y[zi]-min(PE$why))/(diff(range(PE$why))))
                ex1  = round(numf*(zloc$x[zi]-min(PE$x))/(diff(range(PE$x))))

                why2 = round(numc*(zloc$y[zi+1]-min(PE$why))/(diff(range(PE$why))))
                ex2  = round(numf*(zloc$x[zi+1]-min(PE$x))/(diff(range(PE$x))))

                
                MAT = faha$img[ex1:ex2, why1:why2  ]
                w = which.max(MAT)
                c1 = col(MAT)[w]+why1
                r1 = row(MAT)[w]+ex1

                PE$x[r1]
                PE$y[c1]
                points(PE$x[r1], PE$why[c1])
                abline(v=PE$x[r1])

                
                MYx[j] = PE$x[r1]
                MYy[j] = PE$why[c1]
                j = j+1
                
              }

            PEAXBoxes = list(x=MYx, y=MYy)
            
            rDUMPLOC(PEAXBoxes)
             zloc = list(x=NULL, y=NULL)
            
          }


         
        
        
      }
    
    print("DONE with Wlet")
    
    
  }

