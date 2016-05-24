library(XML)

opts = xpathSApply(doc, "//p[@class='level0']/a[@name]/following-sibling::span[@class='nroffip']/text()", xmlValue)
opts = tolower(gsub("_", ".", gsub("CURLOPT_", "", opts)))

library(RCurl)
actual = listCurlOptions()
