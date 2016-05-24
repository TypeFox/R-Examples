#Example of use of text in UTF-8 encoding
library(RGIFT)


Sys.setlocale(locale="es_ES.iso-8859-1")
txt<-"\xbfEs Madrid la capital de Espa\xf1a?"

#Get Local enconding:
#Option 1: Set locale
#locale<-"es_ES.iso-8859-15"
#Option 2: Get locale from OS
locale<-Sys.getlocale("LC_CTYPE")

encod<-strsplit(locale, "\\.|@")[[1]][2]


Sys.setlocale(locale="es_ES.UTF-8")

txt2<-iconv(txt, encod, "UTF8")

sink("utf8.txt")
GIFTTF(txt2, TRUE)
sink()

Sys.setlocale(locale=locale)

