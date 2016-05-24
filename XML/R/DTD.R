dtdIsAttribute <-
function(name, element, dtd)
{
 if(!inherits(element,"XMLElementDef")) {
   element <- dtdElement(as.character(element), dtd)
 }

# return(!is.na(amatch(name, names(element$attributes))))
 return(!is.na(match(name, names(element$attributes))))
}

dtdValidElement <-
#
# checks whether an XML element named `name'
# can be inserted into an element named `within'
# as defined in the specific DTD, optionally
# specifying the position the `name' element would
# be added.
#
# Ideally, this would be used when writing to an XML stream
# (doesn't exist in R or S, yes).
# The stream would monitor the currently open tags
# (as a stack) and would be able to test whether a new 
# insertion was valid.

function(name, within, dtd, pos=NULL)
{

 el <- dtdElement(within, dtd)
 if(is.null(el))
     stop(paste("No such element \"",within,"\" in DTD",sep="", collapse=""))

 return(dtdElementValidEntry(el, name,pos=pos))
}

dtdElementValidEntry <-
function(element, name, pos=NULL)
{
 UseMethod("dtdElementValidEntry", element) # , name, pos)
}

dtdElementValidEntry.XMLElementDef <-
function(element, name, pos=NULL)
{
 return(dtdElementValidEntry(element$contents,name,pos=pos))
}

dtdElementValidEntry.XMLOrContent <-
function(element, name, pos=NULL)
{
 for(i in element$elements) {
   if(dtdElementValidEntry(i, name, pos=pos))
     return(TRUE)
 }

 return(FALSE)
}

dtdElementValidEntry.XMLElementContent <-
function(element, name, pos=NULL)
{
 # if there are no sub-element types, then can't be here.
 # Might check this is a PCDATA by looking at the type.
 if(is.null(element$elements)) {
  return(FALSE)
 }

 return( any(element$elements == name) )
}

dtdElementValidEntry.character <-
function(element, name, pos=NULL)
{
 return(element == name)
}

dtdElementValidEntry.XMLSequenceContent <-
function(element, name, pos=NULL)
{
 if(!is.null(pos)) {
   tmp <- element$elements[[as.integer(pos)]]
   if(!is.null(tmp))
      return(dtdElementValidEntry(tmp))
   else
     return(FALSE)
 }

 for(i in element$elements) {
   if(dtdElementValidEntry(i, name)) {
     return(TRUE)
   }
 }

 return(FALSE)
}

xmlContainsEntity <-
#
# Determine if a particular entity is defined
# within the DTD.
#
function(name, dtd)
{
 return(!is.na(match(name,dtd$entities)))
}

xmlContainsElement <-
#
# Determine if a particular entity is defined
# within the DTD.
#
function(name, dtd)
{
 return(!is.na(match(name,dtd$element)))
}


dtdEntity <-
#
# Retrieves the specified entity from the DTD definition.
# Uses the `dtd$entitities' list.
#
function(name, dtd)
{
 dtd$entities[[name]]
}

dtdElement <-
#
# Retrieves the specified element from the DTD definition.
# Uses the `dtd$elements' list.
function(name, dtd)
{
 dtd$elements[[name]]
}
