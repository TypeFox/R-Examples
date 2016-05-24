



#'  Locate position of first occurence of a pattern
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{CHAR.INDEX} Function. 
#'
#' @usage computeChar_index (x,pattern = NULL, split = 0)
#' @param x atomic character or character vector.
#' @param pattern atomic character to look for.
#' @param split atomic numeric. Number of parts pattern to divide to.  
#'
#' @return Numeric. Position of the first occurence of the \code{pattern}.
#' @details \code{computeChar_index(x="Hello user", pattern="user")} Returns a number indicating the character position of the first occurrence of the first letter of 'user' in 'Hello user'. The optional third argument, \code{split}, is a number, which must be a positive integer used to divide 'user' into separate strings. The \code{split} value must be able to divide the \code{pattern} string without remainder.
#'
#' @author Bastian Wiessner
#' @seealso \code{\link{str_locate_all}}
#' @importFrom stringr str_length str_locate_all str_sub str_locate
#' @keywords internal
#' @examples 
#'
#' xpssCompute(x="there is no letter in word", fun="computeChar_index", pattern="letter")
#'
#' xpssCompute(x="there is no letter in word", fun="computeChar_index", pattern="string")
#' 
#' xpssCompute(x=c("fruits are sweet", "fruits are sour","fruits are salty"), 
#' fun="computeChar_index", pattern="sweet")
#' @keywords internal 
#' @export 

computeChar_index <- function(x,pattern = NULL, split = 0){
  options(warn=-1)
  if(!(is.numeric(split))){
    stop("argument split has to be numeric")
  }
  if(is.null(x) || is.null(pattern)){
    stop("argument x and pattern are empty.")
  }
  if(!(is.character(x)) || !(is.character(pattern))){
    stop("argument x and pattern has to be character")
  }
  out <- x
  if(split > 0) {
    if((str_length(pattern) %% split) == 0) {
      split_length <- str_length(pattern) / split
      start_string <- 1
      end_string <- split_length
      splitpattern <- vector()
      for(i in 1:split_length) {
        splitpattern[[i]] <- str_sub(pattern, start = start_string, end = start_string+split-1)
        start_string <- start_string+split
      }
      for(i in 1:length(x)){
        for(j in 1:length(splitpattern)){
          if(sum(as.numeric(grepl(pattern=splitpattern[j],x=x[[i]])))>=1){
              if(min(str_locate(x[[i]],splitpattern[j]),na.rm=T)<out[[i]]){
                out[[i]] <- min(str_locate(x[[i]],splitpattern[j]))  
              }
              
          }  
        }
      }
      out <- as.numeric(out)
      if(is.na(unique(out[[1]]))){
        pos <- which(is.na(out))
        out[pos] <- 0
      }
    } else {
      stop("Rest has to be 0")
    }
  } else {
    out <- unlist(str_locate_all(string=x,pattern=pattern))
    if(length(out) == 0){
      out <- NA
    } else {
      out <- min(out)
    }
  }
  options(warn=0)
  return(out)
}



#' Expand strings on the left
#'
#'Helper Function for xpssCompute. R Implementation of the SPSS \code{CHAR.LPAD} Function. 
#'
#' @usage computeChar_lpad (x, length = NULL, fill = NULL)
#' @param x atomic character or character vector.
#' @param length atomic numeric. Number of characters x is to be filled on the left. 
#' @param fill atomic character. String which should be append on the left side.
#'
#' @return String, left-padded by length \code{lengthpa}.   
#' 
#' @details The value of length represents the number of characters and must be a positive integer. If the optional argument \code{fill} is not specified, \code{x} is padded with blank spaces.
#' @author Bastian Wiessner
#' @seealso \code{\link{paste0}}
#' @importFrom stringr str_length 
#' @keywords internal
#' @examples 
#'
#' xpssCompute(x="My Friend",fun="computeChar_lpad",length=15)
#' 
#' xpssCompute(x="My Friend",fun="computeChar_lpad",length=15, fill="Hello ")
#'
#' @export 

computeChar_lpad <- function(x, length = NULL, fill = NULL){
  if(is.null(length)){
    stop("argument length is empty.")
  }
  if(!(is.numeric(length))){
    stop("argument length has to be numeric.")
  }
  out <- x
  if(is.null(fill)){
    
    replicates <- length-str_length(x)
    for(i in 1:replicates){
      fill <- paste0(" ",fill)
    }
    out <- paste0(fill,out)
  }else{
    for(i in 1:length(x)){
      for(j in 1:trunc(length/str_length(fill))){
        if(str_length(out[i])<length)
        out[i] <- paste0(fill,out[i])  
      }
      if(str_length(out[i])>length){
        out[i] <- substr(out[i],start=0,stop=length)
      }
    }
  
  }

  return(out)
}



#' Length of a string in characters
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{CHAR.LENGTH} Function. 
#'
#' @usage computeChar_length (x)
#' @param x atomic character or character vector.
#'
#' @return Numeric. Returns the length of \code{x} in characters, with any trailing blanks removed.
#' @author Bastian Wiessner
#' @importFrom stringr str_length str_trim
#' @seealso \code{\link{str_length}}
#' @keywords internal
#' @examples 
#'
#'
#' xpssCompute(x="            please remove trailing blanks",fun="computeChar_length")
#'
#' xpssCompute(x="please remove trailing blanks",fun="computeChar_length")
#' 
#' xpssCompute(x=c("please","remove","trailing","blanks"),fun="computeChar_length")
#' @keywords internal
#' @export 


computeChar_length <- function(x){
  x <- str_trim(x)
  out <- str_length(x)
  return(out)
}



#' Byte per character or sign
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{CHAR.MBLEN} Function. 
#'
#' @usage computeChar_mblen (x,pos = NULL)
#' @param x character vector, or a vector to be coerced to a character vector. Giving a factor is an error.
#' @param pos position of character or sign the number of bytes return to. 
#'
#'
#' @return Numeric. Returns the number of bytes in the character at position \code{pos} in \code{x}.
#' @details Important for Asian languages, where a character can fill more than one byte.
#' @author Bastian Wiessner
#' @seealso \code{\link{nchar}}
#' @keywords internal
#' @examples 
#'  
#'#Returns 1 cause the letter "R" fills one byte.
#' xpssCompute(x="R is great!", fun="computeChar_mblen", pos=2)
#' 
#' @keywords internal
#' @export 


computeChar_mblen <- function(x,pos = NULL){
  if(is.null(pos)){
    stop("argument pos is empty.")
  }
  if(!(is.numeric(pos))){
    stop("argument pos has to be numeric")
  }
  out <- substr(x,start=pos,stop=pos)
  out <- nchar(x=out,type="bytes")  
  return(out)
}



#' Locate position of last occurence of a pattern
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{CHAR.RINDEX} Function. 
#'
#' @usage computeChar_rindex (x,pattern = NULL, split = 0)
#' @param x atomic character or character vector.
#' @param pattern atomic character to look for.
#' @param split numeric. Number of parts pattern to divide to.  
#'
#' @return Numeric. Position of the last occurence of the \code{pattern}
#' @details \code{computeChar_index(x="Hello user", pattern="user")} Returns a number indicating the character position of the last occurrence of the first letter of 'user' in 'Hello user'. The optional third argument, \code{split}, is a number, which must be a positive integer used to divide 'user' into separate strings. The \code{split} value must be able to divide the \code{pattern} string without remainder.
#' @author Bastian Wiessner
#' @seealso \code{\link{str_locate_all}}
#' @importFrom stringr str_length str_sub str_locate
#' @keywords internal
#' @examples 
#'
#' xpssCompute(x="At the end i'm looking for a good end", fun="computeChar_rindex", pattern="end")
#'
#' @export 


computeChar_rindex <- function(x,pattern = NULL, split = 0){
  options(warn=-1)
  if(!(is.numeric(split))){
    stop("argument split has to be numeric")
  }
  if(is.null(x) || is.null(pattern)){
    stop("argument x and pattern are empty.")
  }
  if(!(is.character(x)) || !(is.character(pattern))){
    stop("argument x and pattern has to be character")
  }
  out <- vector(mode="numeric",length=length(x))
  if(split > 0) {
    if((str_length(pattern) %% split) == 0) {
      split_length <- str_length(pattern) / split
      start_string <- 1
      end_string <- split_length
      splitpattern <- vector()
      for(i in 1:split_length) {
        splitpattern[[i]] <- str_sub(pattern, start = start_string, end = start_string+split-1)
        start_string <- start_string+split
      }
      for(i in 1:length(x)){
        for(j in 1:length(splitpattern)){
          if(sum(as.numeric(grepl(pattern=splitpattern[j],x=x[[i]])))>=1){
            if(max(str_locate_all(x[[i]],splitpattern[j])[[1]],na.rm=T)>out[[i]]){
              out[[i]] <- max(str_locate_all(x[[i]],splitpattern[j])[[1]],na.rm=T)
            }
          }  
        }
      }
      out <- as.numeric(out)
      if(is.na(unique(out[[1]]))){
        pos <- which(is.na(out))
        out[pos] <- 0
      }
    } else {
      stop("Rest has to be 0")
    }
  }else {
    result <-  str_locate_all(string=x,pattern=pattern)    
    for(j in 1:length(out)){
      if(length(result[[j]][nrow(result[[j]])]) == 0){
        out[[j]] <- NA
      } else {
        out[[j]] <-  result[[j]][nrow(result[[j]])]
      }
    }    
    out <- unlist(out)
  }
  return(out)
}


#' Expand strings on the right
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{CHAR.RPAD} Function. 
#'
#' @usage computeChar_rpad (x, length = NULL, fill = NULL)
#' @param x atomic character or character vector.
#' @param length numeric. Number of characters x is to be filled on the right. 
#' @param fill optional. String which x should be filled with.
#'
#' @return String, right-padded by length \code{length}.   
#' 
#' @details The value of length represents the number of characters and must be a positive integer. If the optional argument \code{fill} is not specified, \code{x} is padded with blank spaces.
#' @author Bastian Wiessner
#' @seealso \code{\link{paste0}}
#' @importFrom stringr str_length 
#' @keywords internal
#' @examples 
#'
#' xpssCompute(x="My Friend", fun="computeChar_rpad", length=15)
#' xpssCompute(x="My Friend", fun="computeChar_rpad", length=15, fill = "Hello")
#' 
#'
#' @export 


computeChar_rpad <- function(x, length = NULL, fill = NULL){
  if(is.null(length)){
    stop("argument length is empty.")
  }
  if(!(is.numeric(length))){
    stop("argument length has to be numeric.")
  }
  out <- x
  if(is.null(fill)){
    
    replicates <- length-str_length(x)
    for(i in 1:replicates){
      fill <- paste0(" ",fill)
    }
    out <- paste0(out,fill)
  }else{
    for(i in 1:length(x)){
      for(j in 1:trunc(length/str_length(fill))){
        if(str_length(out[i])<length)
          out[i] <- paste0(out[i],fill)  
      }
      if(str_length(out[i])>length){
        out[i] <- substr(out[i],start=0,stop=length)
      }
    }
    
  }
  return(out)
}


#'  creates a substring
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{CHAR.SUBSTR} Function. 
#'
#' @usage computeChar_substr (x,pos = NULL, length= NULL)
#' @param x input character vector.
#' @param pos atomic numeric. Indicates the start of the substring.
#' @param length atomic numeric. Specifies the length of the substring.
#'
#' @return String. Returns an shortened string.
#' @author Bastian Wiessner
#' @seealso \code{\link{substr}}
#' @keywords internal
#' @examples 
#'
#' xpssCompute(x=fromXPSS, variables="V1", fun="computeChar_substr", pos = 2, length=3)
#'
#' @export 

computeChar_substr <- function(x,pos = NULL, length= NULL){
  
  if(is.null(pos)){
    stop("pos cant be null")
  }
  out <- substr(x,start=pos,stop=length(x))
  if(!(is.null(length))){
    out <- substr(out,start=0,stop=length)
  }
  return(out)
}





#' computeConcat (x, sep = "")
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{CONCAT} Function. 
#'
#' @usage computeConcat (x, sep = "")
#' @param x input data.
#' @param sep atomic character element which splits the input data. Default is "".
#'
#'
#' @return String. Returns a string that is the concatenation of all its arguments.
#' @author Bastian Wiessner
#' @seealso \code{\link{paste0}}
#' @keywords internal
#' @examples 
#' 
#' data(fromXPSS)
#' xpssCompute(x = fromXPSS, variables = c("V1","V2"), fun="computeConcat")
#'
#' @export 

computeConcat <- function(x,sep=""){
  x <- as.data.frame(x)
  out <- character()
  for(i in 1:nrow(x)){ 
    temp <- x[[1]][i]
    temp <- paste(temp,sep=sep)
    for(j in 2:ncol(x)){
      if(ncol(x)==j){
        temp <- paste0(temp,x[[j]][i])  
      } else{
        temp <- paste0(temp,x[[j]][i])  
      }
    }
    out[i] <- temp
  }
  out <- unlist(out)
  return(out)
}



#' Number of bytes in a string
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{LENGTH} Function. 
#'
#' @usage computeLength (x)
#' @param x character or character vector.
#' @return Numeric. Returns the length of \code{x} in bytes, including all trailing blanks. 
#'
#' @author Bastian Wiessner
#' @seealso \code{\link{nchar}}
#' @keywords internal
#' @examples 
#'
#' xpssCompute(x="trailing blanks matter",fun="computeLength")
#' xpssCompute(x="trailing blanks matter                      ",fun="computeLength")
#' xpssCompute(x="                      trailing blanks matter",fun="computeLength")
#' 
#' @export 


computeLength <- function(x){
  out <- nchar(x,type="bytes")
  return(out)
}



#' Convert upper-case letters to lower-case 
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{LOWER} Function. 
#'
#' @usage computeLower (x)
#' @param x a character vector, or an object that can be coerced to character by as.character.
#'
#' @return String. Returns unput with uppercase letters changed to lowercase. The argument can be a string variable or a value. 
#' @author Bastian Wiessner
#' @seealso \code{\link{tolower}}
#' @keywords internal
#' @examples 
#'
#' xpssCompute(x="MAKE ME small PLEASE", fun="computeLower")
#'
#' @export 


computeLower <- function(x){
  out <- tolower(x)
  return(out)
}



#' Trims string on left side
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{LTRIM} Function. 
#'
#' @usage computeLtrim (x,trim = NULL)
#' @param x input character vector
#' @param trim single quoted character or an expression that yields a single character. 
#'
#'
#' @return String. Returns the input string removed by \code{trim} on the left. If trim is not specified, leading blanks are removed. 
#' @author Bastian Wiessner
#' @importFrom stringr str_locate str_trim
#' @seealso \code{\link{str_trim}}
#' @keywords internal
#' @examples 
#'
#' xpssCompute(x="           inputstring", fun="computeLtrim")
#'
#' @export 

computeLtrim <- function(x,trim = NULL){
  if(!(is.null(trim))){
    if(str_length(trim)>1){
      stop("trim argument is limited to one value")
    }  
  }
  if(is.null(trim)) {
    out <- str_trim(x,side="left")  
  } else {
    out <- character()
    x <- str_split(x,pattern = " ")
    for(i in 1:length(x[[1]])){
      if(length(str_locate_all(x[[1]][[i]],trim)[[1]])>0){
        if(min(str_locate_all(x[[1]][[i]],trim)[[1]][,1]) == 1){
          temp <- substring(x[[1]][[i]],2,nchar(x[[1]][[i]]))     
          out <- paste(out,temp)
        } else {
          out <- paste(out,x[[1]][[i]])
        }
      } else {
        out <- paste(out,x[[1]][[i]])
      }
    }
  }
  return(out)
}




#' Return the input data without removing trailing blanks
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{NTRM} Function. 
#'
#' @usage computeNtrim (x)
#' @param x input character vector
#' @return String. Returns the values of the input data
#' @author Bastian Wiessner
#' @keywords internal
#' @examples 
#'
#' xpssCompute(x=fromXPSS,variables="V1", fun="computeNtrim")
#'
#' @export 

computeNtrim <- function(x){
  return(x)
}




#' Replace matched pattern in a string
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{REPLACE} Function. 
#'
#' @usage computeReplace (x,pattern=NULL,match=NULL, count = NULL)
#' @param x input vector 
#' @param pattern pattern to look for
#' @param match string tp replace with
#' @param count number of occurencees of pattern to replace
#' @return String. Returns a character vector
#' @author Bastian Wiessner
#' @seealso \code{\link{str_replace}} \code{\link{str_replace_all}} \code{\link{str_sub<-}}
#' @importFrom stringr str_replace_all str_locate_all str_sub<-
#' @keywords internal
#' @examples 
#'
#' xpssCompute(x="Makfts sftnsft", fun="computeReplace",pattern="ft",match="e")
#'
#' @export 

computeReplace <- function(x,pattern=NULL,match=NULL, count = NULL){
  if(is.null(match) || is.null(pattern))
  {
    stop("pattern and match cant be null")
  }
  if(is.null(count)){
    out <-  str_replace_all(string=x,pattern=pattern,replacement=match)      
  }else {
    position <- str_locate_all(string=x,pattern)
    for(i in 1:count){
      str_sub(string=x,start=position[[1]][,1][[i]],end=position[[1]][,2][[i]]) <- match
    }
    out <- x
  }
  return(out)
}



#' Trims string on right side
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{RTRIM} Function. 
#'
#' @usage computeRtrim (x,trim = NULL)
#' @param x input character vector
#' @param trim single quoted character or an expression that yields a single character. 
#'
#'
#' @return String. Returns the input string removed by \code{trim} on the right. If trim is not specified, trailing blanks are removed. 
#' @author Bastian Wiessner
#' @seealso \code{\link{str_trim}}
#' @importFrom stringr str_trim str_locate
#' @keywords internal
#' @examples 
#'
#' xpssCompute(x="inputstring          ", fun="computeLtrim")
#'
#' @export 

computeRtrim <- function(x,trim = NULL){
  if(!(is.null(trim))){
    if(str_length(trim)>1){
      stop("trim argument is limited to one value")
    }  
  }
  if(trim == " ") {
    out <- str_trim(x,side="right")  
  } else {
    out <- character()
    x <- str_split(x,pattern = " ")
    for(i in 1:length(x[[1]])){
      if(length(str_locate_all(x[[1]][[i]],trim)[[1]])>0){
        if(max(str_locate_all(x[[1]][[i]],trim)[[1]][,1]) == str_length(x[[1]][[i]])){
          temp <- substring(x[[1]][[i]],1,nchar(x[[1]][[i]])-1)     
          out <- paste(out,temp)
        } else {
          out <- paste(out,x[[1]][[i]])
        }
      } else {
        out <- paste(out,x[[1]][[i]])
      }
    }
  }
  return(out)
}




#' Truncating strings
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{STRUNC} Function. 
#'
#' @usage computeStrunc (x, length = NULL)
#' @param x input character vector.
#' @param length length x is truncated to
#'
#'
#' @return String. Returns \code{x} truncated to \code{length} (in bytes) and trimmed of remaining blanks. 
#' @author Bastian Wiessner
#' @importFrom stringr str_sub
#' @seealso \code{\link{str_sub}}
#' @keywords internal
#' @examples 
#'
#' xpssCompute(x="Hello    ", fun="computeStrunc", length=7)
#'
#' @export 


computeStrunc <- function(x, length = NULL){
  if(is.null(length)){
    stop("length cant be null")
  } else{
    if(!(is.numeric(length))){
      stop("length has to be numeric")
    }
  }  
  out <- str_sub(x,start=1,end=length)
  out <- str_trim(string = out,side = "right")
  return(out)
}



#' Convert lower-case letters to upper-case 
#'
#'  Helper Function for xpssCompute. R Implementation of the SPSS \code{UPCASE} Function. 
#'
#' @usage computeUpcase (x)
#' @param x a character vector, or an object that can be coerced to character by as.character.
#'
#'
#'
#' @return String. Returns unput with lowercase letters changed to uppercase. The argument can be a string variable or a value. 
#' @author Bastian Wiessner
#' @seealso \code{\link{toupper}}
#' @keywords internal
#' @examples 
#'
#' xpssCompute(x="make it big", fun="computeUpcase")
#' 
#' @export 


computeUpcase <- function(x){
  out <- toupper(x)
  return(out)
}


