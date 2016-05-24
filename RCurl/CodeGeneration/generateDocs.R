curl = htmlParse("http://curl.haxx.se/libcurl/c/curl_easy_setopt.html")
els = getNodeSet(curl, "//p[@class='level0' and ./a[@name]]")
topics = sapply(els, function(x) xmlGetAttr(x[[1]], "name"))

i = !sapply(topics, is.null)
els = els[i]
names(els) = topics[i]

##########################


ps = getNodeSet(curl, "//table//h1/following-sibling::p")

ps = getNodeSet(curl, "//table//h1/following-sibling::p[./a][1]")
first = getNodeSet(curl, "//table//h1/following-sibling::p[@class='level0' and ./a[position() = 1]]")


p.level = sapply(ps, xmlGetAttr, "class")

idx = cumsum(p.level == "level0")


#############################################################
#
# This works, somewhat inlegenantly.
# We find all the nodes that identify a parameter. Then we walk
# the sibling list until we find the next one or a section/grouping
# header.

first = getNodeSet(curl, "//table//h1/following-sibling::p[@class='level0' and ./span[starts-with(string(.), 'CURLOPT_')]]")
names(first) = tolower(gsub("^CURLOPT_", "", sapply(first, function(x) xmlValue(x[["span"]]))))

isSectionHeading =
function(node)
{
 (names(node) == 'a' && xmlValue(node) == "" ) || (xmlName(node) == "h2" && xmlGetAttr(node, "class") == "nroffsh")
}

f =
function(node)
{
  sib = getSibling(node)
  txt = character()
  while(!is.null(sib) && xmlGetAttr(sib, "class", "") != "level0" && !isSectionHeading(sib)) {
     txt = c(txt, xmlValue(sib))
     sib = getSibling(sib)
  }
  paste(txt, collapse = "\n")
}

docs = sapply(first, f)

docs = gsub("char *",  "string", docs, fixed = TRUE)

docText = gsub("CURLOPT_([A-Z_]+)", "\\\\code{\\L\\1}", docs, perl = TRUE)

x = sprintf("\\item{%s}{%s}", names(docText), docText)

cat("\\itemize{", x, "}", "\n", sep = "\n")



       

