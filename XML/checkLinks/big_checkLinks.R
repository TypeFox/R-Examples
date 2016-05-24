library(XML)

checkHTMLDocLinks =
function(fileName)
{  
   doc = htmlParse(fileName, error = function(...){})
   href = unlist(getNodeSet(doc, "//a/@href"))
   ng = grep("^(#|/|http|ftp|mailto)", href, invert = TRUE, value = TRUE)

   ex = file.exists(ng)
   ans = list(href = character(), object = character())
   if(!all(ex)) {
     ans$href = ng[!ex]
   }

   
   refs = unlist(getNodeSet(doc, "//object/@data|//embed/@src"))
   ex = file.exists(refs)
   if(!all(ex))
     ans$object = refs[!ex]

   internal = gsub("^#", "", grep("^#", href, value = TRUE))
   anchors = unlist(getNodeSet(doc, "//*/@id"))

   ans$internal = setdiff(internal, anchors)

   src = unlist(getNodeSet(doc, "//*/@src"))

   lsrc = grep("^(http|ftp)", src, invert = TRUE, value = TRUE)
   ex = file.exists(lsrc)
   ans$src = lsrc[!ex]

   structure(ans, class = "MissingLinks")
}
