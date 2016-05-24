`xport.file.header` <-
  function( cDate=Sys.time(), mDate=cDate, sasVer="7.00", osType="Unknown" )
  {
    .C("fill_file_header",
       cDate = xport.dateFMT(cDate),           # Creation date
       mDate = xport.dateFMT(mDate),           # Modification date
       sasVer = toupper(as.character(sasVer)), # SAS version number
       osType = as.character(osType), # Operating System (can include lowercase)
       PACKAGE="SASxport"
       )

    .Call("getRawBuffer", PACKAGE="SASxport")

  }

