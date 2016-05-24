library(XML)
url<-paste("http://allegro.pl/listing.php/search?category=15821&sg=0&p=",1:5,"&string=facebook",sep="")
o = readHTMLTable(url[1], stringsAsFactors = FALSE, which = 5)
o = readHTMLTable(url[1], which = 5, encoding = "UTF-8")
