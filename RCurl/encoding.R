library(XML)
xmlParse("<doc/>", asText = TRUE)

myParse =
  function(x)
     parse(text = x)

foo =
function(x) {
 browser()
 iconv(x, "", "UTF8")
}

foo =
function(x) {
 tt = readLines(file(x, encoding = "UTF8"))
 browser()
 tt
}


library(RCurl)
url <- "http://search.twitter.com/search.json?q=%E5%95%A4%E9%85%92&result_type=recent&rpp=1&page=1"
doc <- getURLContent(URLencode(url), .encoding = "UTF-8")
Encoding(doc)

library(RJSONIO)
j <- fromJSON(doc)
tweet <- j$results[[1]]$text
