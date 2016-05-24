writePRJ<-function(spobj,filename) {
  if(!is.na(proj4string(spobj))) {
    if(toupper(substr(filename,nchar(filename)-3,nchar(filename)))!=".PRJ" & filename!="") filename<-paste(filename,".prj",sep="")
    cat(showWKT(proj4string(spobj)),file=filename)
  } else stop("The object has no CRS specified !")
}