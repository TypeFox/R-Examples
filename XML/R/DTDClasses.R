#
# Some methods for the DTD classes, similar in spirit
# to those in XMLClasses
#
#    print()
#
#
#
# XMLSystemEntity
# XMLEntity
# XMLElementDef
# XMLSequenceContent
# XMLOrContent
# XMLElementContent
# XMLAttributeDef
#


print.XMLElementDef <-
function(x, ...)
{
 cat("<!ELEMENT", x$name," ")
 print(x$contents)
 cat(">\n")
 if(length(x$attributes)) {

 cat("<!ATTLIST ", x$name,"\n")
  for(i in x$attributes) {
    cat("\t")
    print(i)
    cat("\n")
  }
  cat(">\n")
 }
}


print.XMLElementContent <-
function(x, ...)
{
 if(names(x$type)[1] == "PCData") {
   cat(" ( #PCDATA ) ")
   return()
 }
 cat("(")
 cat(x$elements)
 cat(")",switch(names(x$ocur)[1],Once="", "One or More"="+","Zero or One"="?","Mult"="*")) 
}


print.XMLOrContent <-
function(x, ...)
{
 n <- length(x$elements)
 cat("( ")
 for(i in 1:n) {
   print(x$elements[[i]])
   if(i < n)
    cat(" | ")
 }
 cat(" )")
}

print.XMLSequenceContent <-
function(x, ...)
{
 cat("( ")
 n <- length(x$elements)
 for(i in 1:n) {
    print(x$elements[[i]])
    if(i < n)
        cat(", ")
 }
 cat(" )")
}


print.XMLAttributeDef <-
function(x, ...)
{
 if(names(x$defaultType)[1] != "Implied")
   dflt <- paste("\"", x$defaultValue,"\"",collapse="",sep="")
 else
  dflt <- ""

 cat(x$name, xmlAttributeType(x), xmlAttributeType(x, TRUE), dflt)
}

xmlAttributeType <-
function(def, defaultType = FALSE)
{

 if(defaultType == FALSE & names(def$type)[1] == "Enumeration") {
   return( paste("(",paste(def$defaultValue,collapse=" | "),")", sep=" ", collapse="") )
 }

 switch(ifelse(defaultType, names(def$defaultType)[1], names(def$type)[1]),
         "Fixed" = "#FIXED",
         "CDATA" = "CDATA",
         "Implied" = "#IMPLIED",
         "Required" = "#REQUIRED",
         "Id" = "#ID",
         "IDRef" = "#IDREF",
         "IDRefs" = "#IDREFS",
         "Entity" = "#ENTITY",
         "Entities" = "ENTITIES",
         "NMToken" = "#NMTOKEN",
         "NMTokens" = "#NMTOKENS",
         "Enumeration" = "",
         "Notation" = "",
         "<BROKEN>"
       )
}


print.XMLEntity <-
function(x, ...)
{
 cat("<!ENTITY %", x$name,paste("\"", x$value,"\"",sep="",collapse=""), ">\n")
}




xmlAttrs.XMLElementDef <-
function(node, ...)
{
 node$attributes
}


if(useS4) {
  setGeneric("xmlAttrs", function(node, ...) standardGeneric("xmlAttrs"))
  setMethod("xmlAttrs", "XMLElementDef", xmlAttrs.XMLElementDef)
}
