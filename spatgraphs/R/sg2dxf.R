#' sg to dxf format
#'
#' @param g sg object
#' @param x pattern object used for computing g
#' @param file filename for output
#'
#' @export

sg2dxf <- function(g, x, file) {
  if (!inherits(g, "sg"))
    stop("'g' should be of class sg")
  if (substr(file, nchar(file)-3, nchar(file))!=".dxf")
    file<-paste(file, ".dxf", sep="")

  x <- sg_parse_coordinates(x)

  text<-paste("SECTION\n2\nENTITIES\n")

  d <- ncol(x)
  if(d > 3) stop("Support only for 2D or 3D patterns at them moment.")
  ## The main part of the file
  #The points
  for(i in 1:g$N) {
    text<-paste(text,"0\nPOINT\n8\n0\n10\n",x[i,1],"\n20\n",x[i,2],"\n")
    if(d==3)text<-paste(text,"\n30\n",x[i,3],"\n")
  }
  #The lines
  for(i in 1:x$N)
  {
    for(j in x$edges[[i]])
    {
      text<-paste(text,"0\nLINE\n8\n1\n10\n",x[i,1],"\n20\n",x[i,2],"\n")
      if(d==3)text<-paste(text,"\n30\n",x[i,3],"\n")
      text<-paste(text,"11\n",x[j,1],"\n21\n",x[j,2],"\n")
      if(d==3)text<-paste(text,"\n31\n",x[j,3],"\n")

    }
  }

  text<-paste(text, "0\nENDSEC\n  0\nEOF\n")

  ## write the file
  cat(text, file=file)
}
