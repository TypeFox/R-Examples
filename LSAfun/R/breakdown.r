#### Transform words for usability

#' @export
breakdown <- function(x){
  
  x <- tolower(x)
  
  ## Umlaute
  
  x <- gsub(x=x,pattern="\xe4",replacement="ae")
  x <- gsub(x=x,pattern="\xf6",replacement="oe")
  x <- gsub(x=x,pattern="\xfc",replacement="ue")  
  
  ## Accents
  
  x <- gsub(x=x,pattern="\xe0",replacement="a")
  x <- gsub(x=x,pattern="\xe1",replacement="a")
  x <- gsub(x=x,pattern="\xe2",replacement="a")
  
  x <- gsub(x=x,pattern="\xe8",replacement="e")
  x <- gsub(x=x,pattern="\xe9",replacement="e")
  x <- gsub(x=x,pattern="\xea",replacement="e")
  
  x <- gsub(x=x,pattern="\xec",replacement="i")
  x <- gsub(x=x,pattern="\xed",replacement="i")
  x <- gsub(x=x,pattern="\xee",replacement="i")
  
  x <- gsub(x=x,pattern="\xf2",replacement="o")
  x <- gsub(x=x,pattern="\xf3",replacement="o")
  x <- gsub(x=x,pattern="\xf4",replacement="o")
  
  x <- gsub(x=x,pattern="\xf9",replacement="u")
  x <- gsub(x=x,pattern="\xfa",replacement="u")
  x <- gsub(x=x,pattern="\xfb",replacement="u")
 
  
  x <- gsub(x=x,pattern="\xdf",replacement="ss")
  
  ## Convert to ASCII
  
  x <- iconv(x,to="ASCII//TRANSLIT")
  
  ## Punctation, Numbers and Blank lines
  
  x <- gsub(x=x,pattern="[[:punct:]]", replacement=" ")
  x <- gsub(x=x,pattern="[[:digit:]]", replacement=" ")
  x <- gsub(x=x,pattern="\n", replacement=" ")
  x <- gsub(x=x,pattern="\"", replacement=" ")

  
  return(x)  
}

