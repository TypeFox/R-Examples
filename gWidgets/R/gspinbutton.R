##' @include guiComponent.R
##'
##'
##'
##' RangeSelector Class
setClass("guiComponentRangeSelector",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

## svalue, svalue<- [, [<-

##' spinbutton class
setClass("gSpinbutton",
         contains="guiComponentRangeSelector",
         prototype=prototype(new("guiComponentRangeSelector"))
         )

##' Spinbutton constructor
##'
##' @export
gspinbutton =function(
  from = 0, to = 10, by = 1,
  length.out = NULL, along.with=NULL,
  value = from, digits = 0,
  handler = NULL, action = NULL, container = NULL, ... ,
  toolkit=guiToolkit()){

  ## mostly from seq.default
  if (!missing(along.with)) {
    length.out <- length(along.with)
  } else if(!missing(length.out)) {
    len <- length(length.out)
    if (!len) 
      stop("argument 'length.out' must be of length 1")
    if (len > 1L) {
      warning("first element used of 'length.out' argument")
      length.out <- length.out[1L]
    }
    length.out <- ceiling(length.out)
  }
  if(!is.null(length.out)) {
    ## set up by to be for length,out
    by <- (to - from)/(length.out[1] - 1)
  }


  widget <- .gspinbutton (toolkit,
                          from=from, to=to, by=by, value=value, digits=digits,
                          handler=handler, action=action, container=container ,...
                          )
  obj <- new( 'gSpinbutton',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gspinbutton
setGeneric( '.gspinbutton' ,
           function(toolkit,
                    from = 0, to = 10, by = 1, value = from, digits = 0,
                    handler = NULL, action = NULL, container = NULL, ... ) standardGeneric( '.gspinbutton' ))
