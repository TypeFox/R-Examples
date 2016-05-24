checkLinks =
function(doc, local = TRUE)
{
  doc = htmlParse(doc)
  links = sapply(getNodeSet(doc, "//a[@href]"), xmlGetAttr, "href")
  if(local)
    links = grep("^(http|mailto)", links, invert = TRUE, value = TRUE)

  cur = getwd()
  on.exit(setwd(cur))
  setwd( dirname(docName(doc)))
  links[!file.exists(links)]
}