xmlSerializeHook =
function(x)
{
   if(inherits(x, c("XMLInternalDocument", "XMLInternalElementNode")))
     c(as(x, "character"), class(x)[1])
   else if(is(x, "XMLNodeSet"))
     c(sapply(x, as, "character"), class(x)[1])   
   else
     NULL
}

xmlDeserializeHook =
function(x)
{
  if(length(x) == 2) {
    if(x[2] == "XMLInternalElementNode")
      xmlRoot(xmlParse(I(x[1])))
    else if(x[2] == "XMLNodeSet") {
        #XXX we should put these into the same document, but it is hard to make this sensible.
      structure(lapply(x, function(x) xmlRoot(xmlParse(I(x)))), "XMLNodeSet")
    } else if(x[2] == "XMLInternalDocument")
       xmlParse(I(x[1]))
    else
       stop("Not sure how to handle ", x[2])
  } else
    xmlParse(I(x))    
}

