assign("export.point",
function(point.obj, filename) {

  if (!inherits(point.obj,"point")) stop('point.obj must be of class, "point".\n')

  cat(paste(paste(names(point.obj),collapse=","),'\n',collapse=""), file=filename)

  line <- NULL
  for (i in 1:(length(point.obj$x))) {
    for (j in 1:length(names(point.obj))) {
      line <- paste(line,(round(point.obj[[j]],digits=7))[i], collapse=",",sep=",")
    }
    cat(paste(line,'\n',collapse="",sep=""),file=filename,append=TRUE)
  }

})
