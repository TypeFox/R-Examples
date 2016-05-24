# -----------------------------------------------------
# AB : feb 2006
# FUNCTION:
# Read the polygon file when it is in format 2
# ARGUMENTS:
# filename:  pathname of the file which contains 
# the pogygons coordinates
# delim: character separator on the file
# RETURN:
# a list of class 'listpoly': each component is an object 'poly'
# which contains one polygon coordinates.
# -----------------------------------------------------
readpoly2 <- function(filename, delim="\t") {
# -----------------------------------------------------
  #  Read the whole file
  lu = readLines(filename)
  
  npoly = as.integer(lu[1])

  retour = list()
  
  k=1
  for (i in 1:npoly) {
    k=k+1
    #  Read the identificator, number of vertices, area
    id = strsplit(lu[k],delim)[[1]][1]
    id = id[id !=""]
    k=k+1
    xlu =  as.double(strsplit(lu[k],delim)[[1]])
    k=k+1
    ylu =  as.double(strsplit(lu[k],delim)[[1]])
    coo = matrix(c(xlu, ylu), ncol=2)
    dimnames(coo) = list(NULL, c("xcoord", "ycoord"))
    retour[[id]] = coo
  } # end npoly
# Add the class and the calling command:
  class(retour) <- "listpoly"
  attr(retour, "call") <- match.call()
  return(retour)
}
