library(RCurl)
sURL <- "http://www.nseindia.com/marketinfo/ipo/ipo_pastissues.jsp?year=ALL#"
tt = getURLContent(sURL, binary = TRUE, useragent = "Mozilla/5.0 (Macintosh; U; Intel Mac OS X 10.6;en-US; rv:1.9.2.12) Gecko/20101026 Firefox/3.6.12")
txt = paste(rawToChar(tt, TRUE), collapse = "")
doc = htmlParse(txt, asText = TRUE)

length(getNodeSet(doc, "//a"))



