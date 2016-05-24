pedcheck <-
function(df){
# pedcheck()  -  check Id, SId, DId in dataframe df are valid for dmm
#             -  check all SId's are male, all DID's female
  errcount=0
# Id
  d <- diff(as.numeric(df$Id)) 
  if(as.numeric(df$Id[1]) != 1) {
    cat("Id's must start at 1:\n")
    errcount <- errcount + 1
  }
  if(d[1] != 1) {
    cat("Id's must be an arithmetic sequence:\n")
    errcount <- errcount + 1
  }
  if(length(unique(d)) != 1) {
    cat("Id's must be unique:\n")
    errcount <- errcount + 1
  }
# SId
  if(any(is.na(match(df$SId[!is.na(df$SId)],df$Id)))){
    cat("SId's must occur as an Id in the dataframe:\n")
    errcount <- errcount + 1
  }
# DId
  if(any(is.na(match(df$DId[!is.na(df$DId)],df$Id)))){
    cat("DId's must occur as an Id in the dataframe:\n")
    errcount <- errcount + 1
  }
# Sex of SID's
  sids <- unique(df$SId)[!is.na(unique(df$SId))]
  sidsex <- df[sids,"Sex"]
  if(any(sidsex != sidsex[1])){
    cat("All SId's must be male:\n")
    errcount <- errcount + 1
  }
# Sex of DID's
  dids <- unique(df$DId)[!is.na(unique(df$DId))]
  didsex <- df[dids,"Sex"]
  if(any(didsex != didsex[1])){
    cat("All DId's must be female:\n")
    errcount <- errcount + 1
  }
  return(errcount)
}
