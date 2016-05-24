printList <-
function(toPrint = letters[1:3]
                      , finalSepWord = "and"
                      , midSep = ","
                      ){
  # set final separator
  finalSep <- ifelse(length(toPrint)>2,
                     paste(midSep, " ",finalSepWord," ",sep=""),
                     paste(" " ,finalSepWord," ",sep=""))
  if(length(toPrint) == 1){
    out <- toPrint
  } else{
    out <- paste(
      paste(toPrint[1:(length(toPrint)-1)],collapse = paste(midSep," ",sep = "")),
      toPrint[length(toPrint)],sep = finalSep
    )
  }
  return(out)
}
