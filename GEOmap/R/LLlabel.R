LLlabel<-function(DD, dir=1, ksec=-1 )
  {
    if(missing(dir))
      { dir = 1 }
    if(missing(ksec))
      { ksec=-1 }

    if(is.character(dir))
      {
        if(any(grepl("ew", dir, ignore.case = TRUE)))
          {
            idir = 1

          }
        if(any(grepl("lon", dir, ignore.case = TRUE)))
          {
            idir = 1

          }
        if(any(grepl("ns", dir, ignore.case = TRUE)))
          {
            idir = 2

          }
        if(any(grepl("lat", dir, ignore.case = TRUE)))
          {
            idir = 2

          }

        dir = idir

      }


    #####    make sure all the lons are rectified

    DD = RPMG::fmod(DD,  360)
    DD[DD>180] =  DD[DD>180]-360
    N = length(DD)

    retstrings = rep(" ", length=N)

    for(i in 1:N)
      {

        sn = sign(DD[i])

        
        dms1 = dms(abs(DD[i]) )

        
       
        
        if(ksec>0)
          {
            string = paste(sep="", format(dms1$d), "\\de", format(dms1$m), "\\fm", format(dms1$s, digits=ksec), "\\sd")
            
          }
        else
          {

            asec = round(dms1$s)

             if(asec==60) { dms1$m=dms1$m+1 ;    dms1$s =0; asec= 0 }

            
            amin = round(dms1$m)
            
            if(amin==60) {  dms1$d=dms1$d+1 ; dms1$m=0  }

           
            string = paste(sep="", format(dms1$d), "\\de", format(dms1$m), "\\fm", format(asec), "\\sd")
            
            if(asec==0)
              { 
                string = paste(sep="", format(dms1$d), "\\de", format(dms1$m), "\\fm")
              }
            if(asec==0 & amin==0)
              { 
                string = paste(sep="", format(dms1$d), "\\de")
              }

          }



        
        if(dir==1)
          {
            if(sn<0) { adir="W" }
            else
              { if(sn==0)adir=" " else adir="E" }

                
            string = paste(sep="",string, adir)

            retstrings[i] = string


          }

       if(dir==2)
          {
            if(sn<0) { adir="S" }
            else
              { if(sn==0)adir=" " else adir="N" }

                
            string = paste(sep="",string, adir)

            retstrings[i] = string


          }



        
        
      }
    

    invisible(retstrings)

  }
