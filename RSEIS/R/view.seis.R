view.seis<-function(aday, ihour, inkhour, SAVEFILE, days, DB, usta, acomp,  STDLAB =c("QUIT",  "NEXT", "PREV", "HALF") , TZ=NULL  )
  {
    ####     aday   = index of which day to use in vector days
    ####     ihour     = hour to start
    ####     inkhour   = increment in hours for viewing panel
    ####     SAVEFILE  = file to save window picks in
    ####     days   =list of years and days to select from days=list()
    ########## for example: days = list(jd=c(365,366, 1,2,3,4), yr=c(2008, 2008, rep(2009, times=4)))

    ####     DB     = data base list of file names and start-times and durations
    ####     usta   = stations to select
    ####     acomp  =  components to select
    ####     STDLAB  labels for prescribed buttons

   if(missing(STDLAB ))
     {
       STDLAB = c("QUIT",  "NEXT", "PREV", "HALF",  "WPIX", "zoom out",
         "refresh", "restore",  "SPEC", "SGRAM" ,"WLET", "FILT", "UNFILT",
         "Pinfo", "WINFO","Postscript")
     }

   if(missing(TZ)) { TZ=NULL }

   if(is.null(attr(DB, "origyr")))  attr(DB, "origyr") = min(DB$yr, na.rm=TRUE)

   if(!is.list(days) )
     {


       
       print("NO days list")
       print("using:")
       
       
       udays = unique(paste(DB$yr, DB$jd))
       sdays =  as.numeric( unlist( strsplit(udays, split=" ") ) )


       ye = sdays[seq(from=1, to=length(sdays), by=2)]
       d  = sdays[seq(from=2, to=length(sdays), by=2)]
       o = order(ye+d/366)

       days = list(yr = ye[o], jd=d[o])

       print("using:")
       print(days)

     }
   else{
     if(length(days$jd)<1 | length(days$yr)<1 )
       {

         udays = unique(paste(DB$yr, DB$jd))
         sdays =  as.numeric( unlist( strsplit(udays, split=" ") ) )
         
         
         ye = sdays[seq(from=1, to=length(sdays), by=2)]
         d  = sdays[seq(from=2, to=length(sdays), by=2)]
         o = order(ye+d/366)
         
         days = list(yr = ye[o], jd=d[o])
         
         
       }


     

   }


   
   
   eday = EPOCHday(days$yr[aday], jd=days$jd[aday], origyr =  attr(DB, "origyr") )
           
    while( ihour >= 0 & ihour < 24)
      {
##########  set the time window slot:  start time (at1) and  end time (at2)

        
        at1 = eday$jday+(ihour)/24
       ####  print(at1)
        cat(paste(days$yr[aday], days$jd[aday],ihour,  eday$jday,  at1), sep="\n")
        
        
        at2 =  at1+inkhour/24

##################   grab data from DB list, glue together station/component pairs
        GH = Mine.seis(at1, at2, DB, usta, acomp)
        GH$TZ=TZ   

        if(length(GH$STNS)<1)
          {
            ihour = ihour+inkhour
            next
          }

        

########   show wiggles
        KOUT = swig(GH, STDLAB=STDLAB)


        if(length(KOUT$WPX$yr)>1) { save.wpix(KOUT, fn= SAVEFILE) }
        
#####  depending on what came out of swig do different things
######   remember to save the picks in a file for later use

        if(KOUT$but == "QUIT") {
          break  }
        if(KOUT$but == "NEXT") {
          ihour = ihour+inkhour;
          next }
        if(KOUT$but == "PREV") {
          ihour = ihour-inkhour;
          next   }
        if(KOUT$but == "HALF") {
          ihour = ihour+inkhour/2;
          next   }

### if  these buttons are not engaged, then increment and move to time slot
        ihour = ihour+inkhour
      }

  }
