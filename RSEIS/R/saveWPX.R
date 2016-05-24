saveWPX<-function(twpx, destdir="." )
  {
    ##########  given a pick data frame, save
    ##########   to disk, based on minimum time
    ######  write the file in the destdir directory (folder)
    
    RDATES = rangedatetime(twpx)

    ########### check  if the twpx dataframe is empty:
    if(RDATES$max$yr == 0 & RDATES$min$yr == 0)
      {
        return()
      }
  
    fout1 = filedatetime(RDATES$min, 0)

    fout2 = paste(fout1,"RDATA", sep="." )

    fout3 = paste(destdir, fout2, sep="/")

    save(file=fout2, twpx)

    invisible(fout2)
  }
