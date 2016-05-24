# -----------------------------------------------------
# AB : 09/2/2006
# FUNCTION:
# Read the polygon file when it is in format 1
# ARGUMENTS:
# filename:  pathname of the file which contains 
# the pogygons coordinates
# delim: character separator on the file
# RETURN:
# a list of class 'listpoly': each component is an object 'poly'
# which contains one polygon coordinates.
# -----------------------------------------------------
readpoly1 <- function(filename, delim=" ") {
# -----------------------------------------------------
  # Read the whole file
  lu = readLines(filename)

  npoly = as.integer(lu[1])
  
  retour = list()
  
  k=1
  for (i in 1:npoly) {
    k=k+1
    # Read the identificator and the x-coordinates
    ligne = strsplit(lu[k],delim)[[1]]
    ligne = ligne[ligne !=""]
    id = ligne[1]
    xlu =  as.double(ligne[2:length(ligne)])
    k=k+1
    # Read the identificator and the y-coordinates
    ligne = strsplit(lu[k],delim)[[1]]
    ligne = ligne[ligne !=""]
    if (ligne[1] != id) {
      stop(paste ("On the lines ", k-1, " and", k,
                  "the identifiers are respect. ", id,
                  " and", ligne[1], "\n"))
    }
    ylu =  as.double(ligne[2:length(ligne)])
    coo = matrix(c(xlu, ylu), ncol=2)
    dimnames(coo) = list(NULL, c("xcoord", "ycoord"))
    retour[[id]] = coo
  } # end npoly
# Add the class and the calling command:
  class(retour) <- "listpoly"
  attr(retour, "call") <- match.call()
  return(retour)
}
