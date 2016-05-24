"[<-.XMLNode" <-
function(x, i, value)
{
  x$children[i] <- value
  if(!is.character(i)) {
       names(x$children)[i] =
               if(inherits(value, "XMLNode"))
                 xmlName(value)
               else
                 sapply(value, xmlName)
  }
 x
}


"[[<-.XMLNode" <-
function(x, i, value)
{
  x$children[[i]] <- value
  if(!is.character(i)) {
       names(x$children)[i] =
               if(inherits(value, "XMLNode"))
                 xmlName(value)
               else
                 sapply(value, xmlName)
  }  
 x
}


append.xmlNode <-
function(to, ...)
{
 UseMethod("append")
}

append.XMLNode <-
function(to, ...)
{
 args <- list(...)
 if(!inherits(args[[1]], "XMLNode") && is.list(args[[1]]))
   args <- args[[1]]
    
 idx <- seq(length(to$children) + 1, length=length(args))

 args = addNames(args)
 
 if(is.null(to$children))
   to$children <- args
 else  {
   to$children[idx] <- args
   names(to$children)[idx] <- names(args)
 }

 to
}

if(FALSE) {
xmlAddChild <-
function(node, child) {
  node$children <- append(node$children, list(child))
  names(node$children) <- sapply(node$children,xmlName)
  node
}
}
