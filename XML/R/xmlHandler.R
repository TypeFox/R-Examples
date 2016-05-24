xmlHandler <- 
function() {
  data <- list()
  startElement <- function(name, atts,...) {
    if(is.null(atts))
      atts <- list()
    data[[name]] <<- atts
  }
  text <- function(x,...) {
    cat("MyText:",x,"\n")   
  }
  comment <- function(x,...) {
    cat("comment", x,"\n")
  }
  externalEntity <- function(ctxt, baseURI, sysId, publicId,...) {
    cat("externalEntity", ctxt, baseURI, sysId, publicId,"\n")
  }
  entityDeclaration <- function(name, baseURI, sysId, publicId,notation,...) {
    cat("externalEntity", name, baseURI, sysId, publicId, notation,"\n")
  }

  foo <- function(x,attrs,...) { cat("In foo\n")}
  return(list(startElement=startElement, getData=function() {data},
               comment=comment, externalEntity=externalEntity,
                entityDeclaration=entityDeclaration,
                text=text, foo=foo))
}
