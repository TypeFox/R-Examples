### Susumu Tanimura <aruminat@gmail.com>
### Time-stamp: <2011-02-17 12:26:43 umusus>
###

sjis2utf8 <- function(x, CP932=TRUE){
  if(!is.character(x)) x <- as.character(x)
  from <- ifelse(CP932,"CP932","SHIFT_JIS")
  return(iconv(x,from=from,to="UTF-8"))
}

eucjp2utf8 <- function(x){
  if(!is.character(x)) x <- as.character(x)
  return(iconv(x,from="EUC-JP",to="UTF-8"))
}

jis2utf8 <- function(x){
  if(!is.character(x)) x <- as.character(x)
  return(iconv(x,from="ISO-2022-JP",to="UTF-8"))
}


