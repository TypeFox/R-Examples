checkInternalLinks =
function(doc)
{  
   to = getNodeSet(doc, "//a/@href")
   to = grep("^#", to, value = TRUE)

   ids = getNodeSet(doc, "//a/@name")
   to[ is.na(to %in% ids) ]
}

checkExternalLinks =
function(doc)
{
   to = getNodeSet(doc, "//a/@href")
   to = grep("^#", to, value = TRUE, invert = TRUE)

}

checkImages =
function(doc)
{
  hf = xpathSApply(doc, "//img", xmlGetAttr, "src")
  hf[!file.exists(hf)]
}
