viewCHAC <-
function( DBnov , gstas, gcomps,sched, stas, buts='GPIX', preFILT=list() , replot=TRUE , kind=2, Iendian=1, BIGLONG=FALSE)
  {

    
   if(missing(buts))
     {
       buts = c("fspread", "ReSet", "SGRAM")

     }
   if(missing(preFILT))
     {
       preFILT = list(ON=FALSE, fl=1/2 , fh=8, type="BP", proto="BU")


     }

   RSEIS::screens(2)
   dev.set(2)

   IDB = RSEIS::infoDB(DBnov)
   
   BUTLAB = c("REPLOT", "DONE",  "PREV", "HALF", "PSEL", 
     "ZOOM.out", "ZOOM.in", "RESTORE", "Pinfo", 
     "FILT", "UNFILT", "GPIX", "SavePF", "SaveCSV",  "editPIX" , "CONTPF", "PickWin",
     "Postscript")

    NSCHED = length(sched)
   AT2shift = sched[2] - sched[1] 
    i = 1 
    while( i>0 & i<= NSCHED  ) 
      {

        
        at1 = sched[i]
        at2 = at1+AT2shift
        ##
        cat(paste(i, at1) , sep="\n" )
        GH = RSEIS::Mine.seis(at1, at2, DBnov , gstas, gcomps, kind=kind, Iendian=Iendian, BIGLONG=BIGLONG)
        GH$sta = stas
        GH$IDB = IDB
     ##   GH$pickfile = PF
     ##   GH$vel = vel
        
        if(preFILT$ON==TRUE) GH = RSEIS::FILT.SEISN(GH,  FILT=preFILT, TAPER=0.1, POSTTAPER=0.1)

        hord = which(GH$COMPS=="V")

        ##############################    SWIG
        gret = RSEIS::swig(GH, sel=hord, STDLAB=BUTLAB, PADDLAB=buts)
        ##############################    SWIG

        
        if(gret$but=="ReSet")
          {

            
            gstas = gret$g$RESET$gs
            gcomps = gret$g$RESET$gc
            next
          }
        
        if(gret$but=="PREV")
          {
            i = i - 1
            next
          }
        if(gret$but=="HALF")
          {
            sched = sched+AT2shift/(2)
            next
          }

         if(gret$but == "QUIT")
           {
             break
           }

####  save any picks to a file on disc for later use
     ####   print(data.frame(gret$g$WPX))


    if( identical(legitWPX(gret$g$WPX),0)  ) {

    ####   print("No left over picks")
      
    }
    else
      {
        if(length(gret$g$WPX$sec)>=2)
          {
          ###  print(paste("Still have unsaved picks",length(gret$g$WPX$sec)) )
            twpx = gret$g$WPX
          ###   print(data.frame(twpx))
            
          ###  whirid =   gret$g$WPX$name=='99999' &  gret$g$WPX$comp=='9' &  gret$g$WPX$phase=="X"
           ##########  99999  is what is used in PickWin indicating no pick
        ###    nona = is.na(twpx$tag) |  whirid
            
         ###   if(length(nona)>0)
          ###    {
          ###      twpx = RSEIS::deleteWPX(twpx, nona)
          ###    }
###  print(data.frame(twpx))
           ##    print(twpx)
             ###print("going to Qrangedatetime")
            A1T = Qrangedatetime(twpx)
            s1 = RSEIS::secdifL(A1T$min,  twpx)


           ### print(data.frame(twpx))


             ### print("going to INITpickfile")
            PF =  INITpickfile(stas=stas, src=NULL, WPX=twpx)

            ##
            if(replot==TRUE)
              {
                ords1 = order(s1)

                osta =  twpx$name[ords1]

                ohoh = list(dist=s1[ords1] , name=osta)
                
                jord =  RSEIS::seisorder(GH, ohoh, VNE="V")
                

                gzoom1 = RSEIS::secdifL(  GH$info, A1T$min  )
                gzoom2 = RSEIS::secdifL(  GH$info, A1T$max  )
                
                WIN=c(min(gzoom1-10), max(gzoom2)+30)

               
                GH$pickfile = PF
                
                hret = RSEIS::swig(GH, sel=jord, APIX=twpx, WIN=WIN, STDLAB=BUTLAB, PADDLAB=buts)
                
              }
           PFoutput(PF, stas = stas, sol = NULL, format = 2)
          }
          }
        
###  Control-c  will break the loop,
#### but Q-button is
###  a gentler way
       
        i = i + 1 
      }
  }
