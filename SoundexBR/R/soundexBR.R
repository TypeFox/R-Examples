#' @encoding UTF-8
#' @title Phonetic-Coding For Portuguese
#' @description This package provides an algorithm for decoding names into phonetic codes, as they are pronounced in Portuguese. The goal is for homophones to be encoded to the same representation so that they can be matched despite minor differences in spelling. The algorithm mainly encodes consonants; a vowel will not be encoded unless it is the first letter. The soundex code resultant consists of a four digits long string composed by one letter followed by three numerical digits: the letter is the first letter of the name, and the digits encode the remaining consonants. Similar consonants are assigned to the same value, for instance, \samp{D}, \samp{T} are each encoded as the number \samp{3}. Therefore, it can be used for indexing names by their sound proximity compared to alphabetical order.
#' @name SoundexBR
#' @docType package
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#' @import parallel
NULL


#' @useDynLib SoundexBR
#' @export
soundexBR <- function(term, BR=TRUE, useBytes = FALSE) {
  x <- as.character(term)
  x <- accent(x)
  if(BR){
    x<- .ptsound(x)
  }else{
    x
  }
  stopifnot(is.logical(useBytes))
  if (!useBytes){ x <- char2int(x)
    r <- .Call("R_soundexBR", x)
    if (!useBytes) int2char(r) else r
  }
}
NULL



#' @encoding UTF-8
#' @title  Interger to Character
#'
#' @description Declares integer inputs as UTF-8 output.
#'
#' @param x A vector consisting of whole numbers.
#' @seealso \link{char2int}, \link[base]{iconv}.
#'
#' @keywords Encoding
#' @examples
#' int2char(c(84,104,105,115,32,115,104,111,117,108,
#' 100, 32, 98, 101, 32, 119, 104, 111, 108, 101,
#' 32, 110, 117, 109, 98, 101, 114, 115))
#'
#' @export
int2char <- function(x) {
  enc2native(sapply(x, intToUtf8)) # declares its output as UTF-8
}
NULL


#' @encoding UTF-8
#' @title  Character to Integer
#'
#' @description Declares characters to integer.
#'
#' @param x A vector consisting of characters.
#' @seealso \link{int2char}, \link[base]{iconv}.
#'
#' @keywords Encoding
#' @examples
#' 	char2int("This should be whole numbers")
#'
#' @export
char2int <- function(x){
  lapply(enc2utf8(x),utf8ToInt) # declares its output as
}
NULL




#' @encoding UTF-8
#' @title Adjusts for common mispelled words
#'
#' @description This is a helper function for soundexBR and not intended to standalone.
#'
#' @param x a vector of words.
#'
#' @export
.ptsound <- function(x){
N<-nchar(x);
x<-ifelse(substr(x,1,2)=="WA",sub("W","V",x),x);
x<-ifelse(substr(x,1,2)=="Wa",sub("W","V",x),x);
x<-ifelse(substr(x,1,2)=="wa",sub("W","V",x),x);
x<-ifelse(substr(x,1,1)=="H",substr(x,2,N),x);
x<-ifelse(substr(x,1,1)=="h",substr(x,2,N),x);
x<-ifelse(substr(x,1,2)=="KA", sub("K","C",x),x);
x<-ifelse(substr(x,1,2)=="Ka", sub("K","C",x),x);
x<-ifelse(substr(x,1,2)=="ka", sub("K","C",x),x);
x<-ifelse(substr(x,1,2)=="KO", sub("K","C",x),x);
x<-ifelse(substr(x,1,2)=="Ko", sub("K","C",x),x);
x<-ifelse(substr(x,1,2)=="ko", sub("K","C",x),x);
x<-ifelse(substr(x,1,2)=="KU",sub("K","C",x),x);
x<-ifelse(substr(x,1,2)=="Ku",sub("K","C",x),x);
x<-ifelse(substr(x,1,2)=="ku",sub("K","C",x),x);
x<-ifelse(substr(x,1,1)=="Y",sub("Y","I",x),x);
x<-ifelse(substr(x,1,1)=="y",sub("Y","I",x),x);
x<-ifelse(substr(x,1,2)=="CE",sub("C","S",x),x);
x<-ifelse(substr(x,1,2)=="Ce",sub("C","S",x),x);
x<-ifelse(substr(x,1,2)=="ce",sub("C","S",x),x);
x<-ifelse(substr(x,1,2)=="CI",sub("C","S",x),x);
x<-ifelse(substr(x,1,2)=="Ci",sub("C","S",x),x);
x<-ifelse(substr(x,1,2)=="ci",sub("C","S",x),x);
x<-ifelse(substr(x,1,2)=="GE",sub("G","J",x),x);
x<-ifelse(substr(x,1,2)=="Ge",sub("G","J",x),x);
x<-ifelse(substr(x,1,2)=="ge",sub("G","J",x),x);
x<-ifelse(substr(x,1,2)=="GI",sub("G","J",x),x);
x<-ifelse(substr(x,1,2)=="Gi",sub("G","J",x),x);
x<-ifelse(substr(x,1,2)=="gi",sub("G","J",x),x);
}
NULL