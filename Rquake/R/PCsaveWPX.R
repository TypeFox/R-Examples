PCsaveWPX<-function(twpx, destdir="." )
  {
##########  given a pick data frame, save
##########   to disk, based on minimum time
######  write the file in the destdir directory (folder)
   
    
    if( identical(legitWPX(twpx),0)  )
      {

        cat("No legitimate picks", sep="\n")
      }
    else
      {

        RDATES = Qrangedatetime(twpx)

########### check  if the twpx dataframe is empty:
        if(RDATES$max$yr == 0 & RDATES$min$yr == 0)
          {
            return()
          }
        
        fout1 = PCfiledatetime(RDATES$min, 0)

        fout2 = paste(fout1,"RDATA", sep="." )

        if( !identical(destdir, ".") )
          {
            fout3 = paste(destdir, fout2, sep="/")
          }
        else
          {

            fout3 = fout2
          }

        save(file=fout3, twpx)

        invisible(fout2)
      }


    invisible(NULL)
  }
