##' @include methods.R
NULL

##' dialog for file and directory selection
##'
##' @param text initial text
##' @param type type of browser: to open a file, to save a file or to select a directory
##' @param initial.filename Suggested file name
##' @param initial.dir initial directory. If a filename is given, and is not an absolute name, this will be prepended. If filename given initial directory will be taken from that.
##' @param filter A filter specifiation. This can be a named character vector of file extensions or something toolkit specific.
##' Here are some examples:
##' \itemize{
##' \item{character}{c("csv"="csv","txt"="txt")}
##' \item{RGtk2}{ Something like
##'   \preformatted{
##' list("All files" = list(patterns = c("*")),
##'      "R files" = list(patterns = c("*.R", "*.Rdata")),
##'      "text files" = list(mime.types = c("text/plain"))
##'     )
##'   }
##' }
##' \item{tcltk}{}
##' \item{Qt}{}
##' }
##' @param multi Logical. Allow multiple files to be selected?
##' @inheritParams gwidget
##' @return returns filename(s) or \code{character(0)} if no selection.
##' @export
gfile <- function(
                  text = "",
                  type = c("open", "save", "selectdir"),
                  initial.filename = NULL,
                  initial.dir=getwd(),  # initial directory 
                  filter = list(),
                  multi=FALSE, ...,
                  toolkit=guiToolkit()){

  deprecated_args <- list(initialfilename=c("Renamed to initial.filename"),
                          handler=c("No handler argument, use return value instead","(handler(gfile), not gfile(handler)"),
                          action=c("Removed as argument with handler")
                          )
  check_deprecated(deprecated_args, ...)

  type <- match.arg(type)
  
  if(!is.null(initial.filename)) {
    if(!file.exists(initial.filename))
      initial.filename <- paste(initial.dir, initial.filename, sep=.Platform$file.sep)
    else
      intial.dir <- dirname(initial.filename)
  }
  
  val <- .gfile (toolkit,
                 text=text, type=type, initial.filename=initial.filename,
                 initial.dir = initial.dir,
                 filter=filter, multi=multi, ...
                 )
  Encoding(val) <- "UTF-8"
  as.character(val)                     # coerce to character
}


##' generic for toolkit dispatch
##'
##' @export
##' @rdname gfile
.gfile <-  function(toolkit,
                    text = "", type = c("open", "save", "selectdir"),
                    initial.filename = NULL,
                    initial.dir = getwd(),
                    filter = list(), multi=FALSE,
                    ... )
           UseMethod( '.gfile' )



##' constructor for file/directory selection widget
##'
##' Basically an entry box instance with a button to initiate \code{gfile}.
##' @param quote quote output
##' @return Returns an object of class \code{gFilebrowse}. This should
##' inherit the methods  of \code{gedit} instances.
##' @export
##' @rdname gfile
gfilebrowse <- function (
                         text = "Select a file...",
                         type = c("open", "save", "selectdir"),
                         initial.filename = NULL,
                         initial.dir = getwd(),
                         filter = list(),
                         quote = TRUE,
                         handler=NULL, action=NULL,
                         container = NULL, ..., toolkit = guiToolkit()) {


  
  obj <- .gfilebrowse (toolkit,
                       text=text, type=type,
                       initial.filename=initial.filename,
                       initial.dir=initial.dir,
                       filter=filter,
                       quote=quote,
                       handler=handler, action=action, container=container, ...)
  check_return_class(obj, "GFileBrowse")


  obj
}

##' generic for toolkit dispatch
##'
##' @rdname gfile
##' @export
.gfilebrowse <- function(toolkit,
                         text = "Select a file...",
                         type = c("open", "save", "selectdir"),
                         initial.filename = NULL,
                         initial.dir = getwd(),
                         filter = list(),
                         quote=TRUE,
                         handler=NULL, action=NULL,
                         container = NULL, ...)
           UseMethod( '.gfilebrowse' )
