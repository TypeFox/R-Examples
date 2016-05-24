`focleg` <-
function(i)
  {

 #####    FKINDS = c("strikeslip",
  #####    "rev-obl strk-slp",
  #####    "oblique reverse",
  #####    "reverse",
  #####    "norm-oblq strkslp",
  #####    "oblq norm",
  #####    "normal")

    FKINDS =  c("STRIKESLIP" ,       "REV-OBL STRK-SLP" , "OBLIQUE REVERSE",
 "REVERSE" ,          "NORM-OBLQ STRKSLP", "OBLQ NORM",
 "NORMAL")

    return(FKINDS[i])
  }

