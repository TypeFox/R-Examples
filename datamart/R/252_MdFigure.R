#' MdFigure
#' 
#' Internal Class representing figures in MdReport
#' 
#' @examples
#' getSlots("MdReport")
#'
#' @seealso \code{\link{mdfigure}}
#'
#' @name MdFigure-class
#' @rdname MdFigure-class
setClass(
    Class="MdFigure", 
    representation=representation(
      name="character",
      resource="character",
      xdata="Xdata"
    ),
    contains="Target"
)

#' Constructor for MdFigure
#'
#' Internal function to create an MdFigure object.
#'
#' @param name        name of the figure, default ''
#' @param xdata       call for creating the figure
#' @param resource    name of the resource
#' @param clss        class name, default 'MdFigure'
#'
#' @return MdFigure
#' @rdname MdFigure-class
mdfigure <- function(name, xdata, resource=name, clss="MdFigure") {
  new(clss, name=name, xdata=xdata, resource=resource)
}

#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,MdFigure,DirectoryLocation-method
setMethod(
  f="put",
  signature=c(target="MdFigure", where="DirectoryLocation"),
  definition=function(target, where, ...) {
    #' set up device
    fname <- target@name
    if(fname=="") fname <- basename(tempfile(pattern="Rplot", tmpdir="", fileext=""))
    fname <- paste(fname, "png", sep=".")
    fname <- file.path(as.character(where), fname)
    png(filename=fname, width=640, height=480, units="px")
    query(target@xdata, target@resource)
    dev.off()
    return(fname)
  }
)

#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,MdFigure,MemoryLocation-method
setMethod(
  f="put",
  signature=c(target="MdFigure", where="MemoryLocation"),
  definition=function(target, where) {
    ifile <- put(target, new("DirectoryLocation", path=tempdir()))
    ext <- "png"
    tfile <- tempfile()
    on.exit(unlink(tfile))
    base64::encode(file.path(tempdir(), basename(ifile)), tfile)
    sprintf("data:image/%s;base64,\n%s", ext, paste(readLines(tfile), collapse = "\n"))
  }
)

