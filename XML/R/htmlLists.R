setGeneric("readHTMLList",
          function(doc, 
                    trim = TRUE, elFun = xmlValue,
                      which = integer(), ...)
             standardGeneric("readHTMLList"))


setMethod("readHTMLList",
           "character",
          function(doc, 
                    trim = TRUE, elFun = xmlValue,
                     which = integer(), encoding = character(), ...) {
             readHTMLList(htmlParse(doc, encoding = encoding), trim, elFun, which, ...)
          })


setMethod("readHTMLList",
           "HTMLInternalDocument",
          function(doc, 
                    trim = TRUE, elFun = xmlValue,
                     which = integer(), ...) {
            lists = getNodeSet(doc, "//ol | //ul | //dl")
            if(length(which))
               lists = lists[which]
            ans = lapply(lists, readHTMLList, trim = trim, elFun = elFun)
            if(length(which) == 1)
              ans[[1]]
            else
              ans
          })

setMethod("readHTMLList",
           "XMLInternalNode",
          function(doc, 
                    trim = TRUE, elFun = xmlValue,
                     which = integer(), ...) {

            if(xmlName(doc) == "dl")
                return(readHTMLDefinitionList(doc, trim, elFun))

            
            ans = unname(sapply(xmlChildren(doc)[!xmlSApply(doc, is, "XMLInternalTextNode")], elFun))

            if(trim) 
              ans = unname(sapply(ans, function(x) if(is.character(x)) trim(x) else x))

            ans
          })

readHTMLDefinitionList =
function(node, trim = TRUE, elFun = xmlValue)
{
  kids = xmlChildren(node)
  structure(sapply(kids[names(node) == "dd"], elFun),
            names = sapply(kids[names(node) == "dt"], elFun))
}
