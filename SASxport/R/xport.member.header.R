`xport.member.header` <- 
function(dfName, cDate=Sys.time(), mDate=cDate, sasVer="7.00", osType="Unknown",
         dfLabel="", dfType="" )
  {
    .C("fill_member_header",
       dfName = toupper(as.character(dfName)), # Name of data set
       sasVer = toupper(as.character(sasVer)), # SAS version number
       osType = as.character(osType),          # Operating System (can include lowercase)
       cDate  = xport.dateFMT(cDate),          # Creation date
       mDate  = xport.dateFMT(mDate),          # modification date
       dfLabel= as.character(dfLabel),         # Data set label
       dfType = as.character(dfType),          # Data set type
       PACKAGE="SASxport"       
       )

    .Call("getRawBuffer", PACKAGE="SASxport")
  }

