viewseis <-
function( DBnov , gstas, gcomps,sched, stas, buts='GPIX', replot=TRUE , kind=0, Iendian=1, BIGLONG=FALSE)
  {
   if(missing(buts))
     {
       buts = c("fspread", "ReSet")

     }


   IDB = RSEIS::infoDB(DBnov)
   
   BUTLAB = c("REPLOT", "DONE",  "PREV", "HALF", "PSEL", 
     "ZOOM.out", "ZOOM.in", "RESTORE", "Pinfo", 
     "FILT", "UNFILT", "GPIX", "SavePF", "RQ" , "CONTPF", "PickWin",
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
        
        ## FH = RSEIS::FILT.SEISN(GH,  FILT=preFILT, TAPER=0.1, POSTTAPER=0.1)

        hord = which(GH$COMPS=="V")
        
        gret = RSEIS::swig(GH, sel=hord, STDLAB=BUTLAB, PADDLAB=buts)
        
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
        
        if(length(gret$g$WPX$sec)>=2)
          {
      
            twpx = gret$g$WPX
            whirid =   gret$g$WPX$name=='99999' &  gret$g$WPX$comp=='9' &  gret$g$WPX$phase=="X"

            
            ##########  99999  is what is used in PickWin indicating no pick

           
            nona = is.na(twpx$tag) |  whirid
            
            if(length(nona)>0)
              {
                twpx = RSEIS::deleteWPX(twpx, nona)
              }

            A1T = Qrangedatetime(twpx)
            s1 = RSEIS::secdifL(A1T$min,  twpx)

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

###            dout = c(A1T$min$yr, A1T$min$jd, A1T$min$hr, A1T$min$mi,  A1T$min$sec)
            
###            fout1 = filedatetime(dout, 0)

###            foutcvs = paste(fout1,"wpx", sep="." )
###            fout2 = paste(fout1,"RDATA", sep="." )

###            save(file=fout2, twpx)
 ###            write.csv(twpx, file =foutcvs, row.names = FALSE)
###  read.csv("foo.csv")

           PFoutput(PF, stas = stas, sol = NULL, format = 2)

          }
        
###  Control-c  will break the loop, but this is
###  a gentler way
       
        i = i + 1 
      }
  }
