##' @include guiComponents.R

##' dialog for file and directory selection
##'
##' @exports
##' @param text initial text
##' @param type type of browser: to open a file, to save a file or to select a directory
##' @param initialfilename Suggested file name
##' @param filter A filter specifiation. This is toolkit specific
##' @param multi Logical. Allow multiple files to be selected?
##' @param handler called when a file or files have been selected
##' @param action passed to handler
##' @param ... ignored
##' @param toolkit toolkit
##' @return returns filename(s) or \code{NA} if no selection. 
##' 
gfile <- function(
                  text = "", type = c("open", "save", "selectdir"),
                  initialfilename = NULL,
                  filter = list("All files" = list(patterns = c("*")), "R files" = list(patterns = c("*.R",          "*.Rdata")), "text files" = list(mime.types = c("text/plain"))          ),
                  multi=FALSE,
                  handler = NULL, action = NULL, ... ,
                  toolkit=guiToolkit()){
  widget =  .gfile (toolkit,
    text=text, type=type, initialfilename=initialfilename,
    filter=filter, multi=multi, handler=handler, action=action ,...
    )
}


##' generic for toolkit dispatch
##'
##' @export
##' @rdname gfile
setGeneric( '.gfile' ,
           function(toolkit,
                    text = "", type = c("open", "save", "selectdir"),
                    initialfilename = NULL,
                    filter = list("All files" = list(patterns = c("*")), "R files" = list(patterns = c("*.R",          "*.Rdata")), "text files" = list(mime.types = c("text/plain"))
                      ),
                    handler = NULL, action = NULL, ... )
           standardGeneric( '.gfile' ))


##' class for a widget to select a file
setClass("gFilebrowse",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )


##' constructor for file/directory selection widget
##'
##' Basically a \code{gedit} instance with a button to initiate \code{gfile}.
##' @export
##' @param text Instructional text
##' @param type type of dialog (see \code{\link{gfile}})
##' @param quote Do we quote value 
##' @param container parent container
##' @param ... passed to \code{add} method of parent
##' @param toolkit toolkit
##' @return Returns an object of class \code{gFilebrowse}. This should
##' inherit the methods  of \code{gedit} instances.
gfilebrowse <- function (
                         text = "Select a file...", type = "open", quote = TRUE, 
                         container = NULL, ..., toolkit = guiToolkit()) {
  widget <- .gfilebrowse (toolkit,
                          text=text, type=type, quote=quote, container=container, ...)
  obj <- new('gFilebrowse',widget=widget,toolkit=toolkit) 
  return(obj)
}

##' generic for toolkit dispatch
##' @alias gfilebrowse
setGeneric(".gfilebrowse",
           function(toolkit,
                    text = "Select a file...", type = "open", quote = TRUE, 
                    container = NULL, ...)
           standardGeneric( '.gfilebrowse' ))
