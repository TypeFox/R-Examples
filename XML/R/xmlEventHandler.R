xmlEventHandler <- 
function() {
  con <- xmlOutputDOM()

  startElement <- function(name, atts,...) {
    con$addTag(name, attrs=atts, close=FALSE)
  }
  endElement <- function(name) {
    con$closeTag(name)
  }
  text <- function(x,...) {
    con$addNode(xmlTextNode(x))  
  }
  comment <- function(x,...) {
    xmlCommentNode(x)
  }
  externalEntity <- function(ctxt, baseURI, sysId, publicId,...) {
    cat("externalEntity", ctxt, baseURI, sysId, publicId,"\n")
  }
  entityDeclaration <- function(name, baseURI, sysId, publicId, notation, ...) {
    # just writing to screen at this point.
    cat("externalEntity", name, baseURI, sysId, publicId, notation,"\n")
  }

 processingInstruction <- function(sys, value) {
   con$addNode(xmlPINode(sys, value))
 }

 list(startElement=startElement, endElement=endElement,
      processingInstruction=processingInstruction,
      text=text,
      comment=comment,
      externalEntity=externalEntity,
      entityDeclaration=entityDeclaration,
      dom=function(){con})
}

