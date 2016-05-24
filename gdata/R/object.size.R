###------------------------------------------------------------------------
### What: Print object size in human readable format - code
###------------------------------------------------------------------------
object.size <- function(...)
{
  structure(sapply(list(...),
                   utils::object.size),
                   class=c("object_sizes", "numeric"))
}

print.object_sizes <- function(x,
                               quote=FALSE,
                               humanReadable=getOption("humanReadable"),
                               standard="IEC",
                               units,
                               digits=1,
                               width=NULL,
                               sep=" ",
                               justify = c("right", "left"),
                               ...)
{
    print(format(x,
                 humanReadable=humanReadable,
                 standard=standard,
                 units=units,
                 digits=digits,
                 width=width,
                 sep=sep,
                 justify=justify),
          quote=quote,
          ...)


    invisible(x)
}

format.object_sizes <- function(x,
                                humanReadable=getOption("humanReadable"),
                                standard="IEC",
                                units,
                                digits=1,
                                width=NULL,
                                sep=" ",
                                justify = c("right", "left"),
                                ...)
{
    if( !missing(units) )
       {
           if (units=="bytes")
               paste(x, "bytes")
           else
               humanReadable(x,
                             standard=standard,
                             units=units,
                             digits=digits,
                             width=width,
                             sep=sep,
                             justify=justify
                             )
       }
    else if( is.null(humanReadable) || humanReadable==FALSE )
        paste(x, "bytes")
    else
        humanReadable(x,
                      standard=standard,
                      units=units,
                      digits=digits,
                      width=width,
                      sep=sep,
                      justify=justify)

}



is.object_sizes <- function(x) inherits(x, what="object_sizes")

as.object_sizes <- function(x)
{
  if(!is.numeric(x) || any(x<0)) stop("'x' must be a positive numeric vector")

  class(x) <- c("object_sizes", "numeric")
  x
}

c.object_sizes <- function(..., recursive=FALSE)
{
  x <- NextMethod()
  if(is.numeric(x)) class(x) <- c("object_sizes", "numeric")
  x
}

###------------------------------------------------------------------------
### object.size.R ends here
