"ujfnev" <-
function(fnev="hu") {

  #############################################################
  # 
  # TITLE:  ujfnev()
  # AUTHOR: SANDOR KABOS
  # DATE:   20 AUG, 2003
  # CALLS:  
  # NEEDS:  
  # NOTES:  GENERATE A UNIQUE NUMBER FOR AN OUTPUT FILENAME
  #			  ENSURING THAT IT DOES NOT ALREADY EXIST
  #			CURRENTLY, THESE NUMBERS ARE CAPPED AT 999
  #			THIS NUMBER IS APPENDED TO THE FILE BASENAME
  #			  PROVIDED BY fnev
  #
  #############################################################
  
  if(fnev=="") {
    return("")
  }
  for (i in 1:999) {
    FNEV <- paste(sep="", fnev, i, ".txt")
    if(!file.exists(FNEV)) {
      return(FNEV)
    }
  }
}

