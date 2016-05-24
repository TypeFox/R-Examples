# $Id: locus.R 1337 2008-04-30 00:54:56Z warnes $

getlocus  <- function(x,...)
{
  if(is.locus(x))
    return(x)
  else if(!is.null(x["locus"]))
        return(x["locus"])
  else if(!is.null(attr(x,"locus")))
       return(attr(x,"locus"))
  else
    NULL
}

getmarker <- getgene <- getlocus

locus  <- function(name, chromosome, arm=c("p","q","long","short",NA),
                   index.start=NULL, index.end=NULL)
  {
    
    object  <-  list()

    if(!missing(name))
      object$name  <- name
    
    if(!missing(chromosome))
      object$chromosome <- chromosome
    
    if(!missing(arm))
      {
        arm  <- match.arg( arm )
        object$arm  <- switch( arm, p="p", q="q", long="p", short="q")
      }
    if(!missing(index.start))
      object$index.start  <- index.start
    if(!missing(index.end))
      object$index.end  <- index.end
    
    class(object)  <- "locus"
    return(object)
  }


gene  <-  function(name, chromosome, arm=c("p","q","long","short"),
                   index.start, index.end=NULL)
{
  object  <- locus(name, chromosome, arm, index.start, index.end)
  class(object)  <- c("gene","locus")
  object
}


marker <- function(name, type,
                   locus.name, bp.start, bp.end=NULL, relative.to=NULL,
                   ...
                   )
{
  if(is.locus(locus.name))
      object <- locus.name
  else
    object  <-  locus(locus.name, ...)

  if(!missing(name))
    object$marker.name  <- name

  if(!missing(type))
    object$type  <- type

  if(!missing(bp.start))
    object$bp.start  <- bp.start

  if(!missing(bp.end))
    object$bp.end  <- bp.end

  if(!missing(relative.to))
    object$relative.to  <- relative.to
  
  class(object)  <- c("marker","locus")
  object
}

is.locus  <- function(x)
    inherits(x, "locus")

is.gene  <- function(x)
    inherits(x, "gene")

is.marker  <- function(x)
    inherits(x, "marker")



as.character.locus  <- function(x,...)
  {
    loc <- paste( x$chromosome, x$arm, x$index.start, sep="" )
    if( !is.null(x$index.end ) && x$index.start != x$index.end )
      loc  <- paste(loc, "-", x$index.end, sep="")
    loc
  }

as.character.gene  <- function(x,...)
  as.character.locus(x,...)

as.character.marker  <- function(x,...)
  {
    loc  <- as.character.locus(x)
    loc  <- paste(loc, ":", x$bp.start, sep="")
    if(!is.null(x$bp.end)) loc  <-  paste(loc, "-", x$bp.end, sep="")
    loc
  }

print.locus  <-  function(x,...)
  {
    cat("Locus: ", x$name, " (", as.character.locus(x), ")\n", sep="" )
  }

print.gene  <-  function(x,...)
  {
    cat("Gene: ", x$name, " (", as.character.locus(x), ")\n", sep="" )
  }

print.marker  <- function(x,...)
  {
    cat("Marker: ", paste(x$name,":",x$marker.name,sep=""),
        " (", as.character.marker(x), ")\tType: ",x$type,"\n", sep="" )
  }


"locus<-" <- function(x,value)
  {
    attr(x,"locus") <- value
    x
  }


"marker<-" <- "gene<-" <-  get("locus<-")
