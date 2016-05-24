zen2han <- function(s){
    ## if(localeToCharset()[1] == "CP932" && substitute(s) != "zenkaku"){
    ##     s <- iconv(unlist(s), from = "CP932", to = "UTF-8")
    ## }
    if(Encoding(s) != "UTF-8")  s <- iconv(s, from = "", to = "UTF-8")
    s <- paste(s, sep='', collapse='')
    y <- sapply(unlist(strsplit(s, split = "")), function(x){
        i <- utf8ToInt(x)
        if(i >= 65281 && i <= 65374){
            return(intToUtf8(i - 65248))
        }else{
            return(x)
        }
    })
    return(paste(y, collapse = ""))
}
## Old function
## zen2han <-
## function(x){
##   paste2 <- function(x,...){paste(x,...,sep='',collapse='')}
##   zen <- paste(unlist(zenkaku),collapse='')
##   han <- paste2(paste2(0:9),paste2(letters),paste2(LETTERS))
##   return(chartr(zen,han,x))
## }

