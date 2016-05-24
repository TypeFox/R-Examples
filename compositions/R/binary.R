# A function to format a number in binary digits
# x the number
# mb the number of binary digits
binary <- function(x,mb=max(maxBit(x,g)),g=2) {
  if( is.character(x) ) x<-unbinary(x)
  if( g==2 )
    do.call(paste,c(sep="",lapply(mb:0,function(i) ifelse(bit(x,i,g=2),"1","0"))))
  else{
    toDigit <- function(x) c(0:9,LETTERS)[x+1]
    do.call(paste,c(sep="",lapply(mb:0,function(i) toDigit(bit(x,i,g=g)))))
  }
}
# Converts a binary character string to a number
unbinary <- function(x,g=2) {
  if( is.numeric(x) )
    return(x)
  nc =nchar(x)
  D = max(max(nchar(x)),2)
  asDigit <- (if(g==2) function(x) as.logical((match(x,c("0","1","F","T"))-1)%%2) else function(x) as.integer((match(x,as.character(c(0:9,LETTERS,10:19,letters)[c(1:g,37:(36+g))]))-1)%%g))
  c(sapply(1:D,function(i) ifelse(i<=nc,asDigit(substring(x,i,i))*g^(nc-i),rep(0,length(x))))%*%rep(1,D)) 
}

# a function to extract a bit from a binary number
# either given as number or as character string
# x the number or string (may be vectors)
# b the bit to be extracted (may be a vector)
bit <- function(x,b,g=2) UseMethod("bit")                       
bit.numeric   <- function(x,b=0:maxBit(x,g),g=2)  {
  erg <- sapply(b,function(b) (x%/% (g^b) %% g ))
  structure((if(g==2) as.logical else as.integer)(erg),dim=dim(erg))
}
bit.character <- function(x,b=0:maxBit(x,g),g=2)  {
  nc = nchar(x)
  asDigit <- (if(g==2) function(x) as.logical((match(x,c("0","1","F","T"))-1)%%2) else function(x) as.integer((match(x,as.character(c(0:9,LETTERS,10:19,letters)[c(1:g,37:(36+g))]))-1)%%g))
  erg <- sapply(b,function(b) ifelse(b<nc,substring(x,nc-b,nc-b),"0"))
  structure(asDigit(erg),dim=dim(erg))
}
"bit<-" <- function(x,b,g=2,value) UseMethod("bit<-",x)                       
"bit<-.character" <- function(x,b=0:maxBit(x,g),g=2,value)  {
  if( length(x) > 1){
    for(i in 1:length(x))
      bit(x[i],b=b,g=g)<-value
    return(paste(x,sep="",collapse=""))
  }
  erg <- bit(x)
  erg[b+1]<-value
  erg<-ifelse(is.finite(erg),erg,0)
  toDigit <- function(x) c(0:9,LETTERS)[x+1]
  paste(sapply(erg,toDigit),collapse="",sep="")
}
"bit<-.numeric" <- function(x,b=0:maxBit(x,g),g=2,value)  {
  if( length(x) > 1){
    for(i in 1:length(x))
      bit(x[i],b=b,g=g)<-value
    return(x)
  }
  erg <- bit(x)
  erg[b+1]<-value
  erg<-ifelse(is.finite(erg),erg,0)
  c( erg %*% 2^(0:(length(erg)-1)))
}

# maxBit The maximum bit used in a binary number  
maxBit <- function(x,g=2) UseMethod("maxBit")
maxBit.numeric <- function(x,g=2) ceiling(log(x+1,g))-1
maxBit.character <- function(x,g=2) max(nchar(x))-1
bitCount <- function(x,mb=max(maxBit(x,g)),g=2) {
  c(bit(x,0:mb) %*% rep(1,mb+1))
}

#%AND% <- function(x,y) UseMethod("%AND%",x)
#%OR% <- function(x,y) UseMethod("%OR%",x)
## A bitwise AND function
#"%AND%.default" <- function(x,y) {
#  mb <- max(max(maxBit(x),maxBit(y)),2)
#  if( is.numeric(x) && is.numeric(y) ) {
#    # Use C-Code
#    sapply(0:mb, function(i) (bit(x,i) & bit(y,i))*2^i )
#  } else {
#    do.call(paste,c(sep="",lapply(mb:0,function(i) ifelse((bit(x,i) & bit(y,i)),"1","0"))))
#  }
#}


# A bitwise OR function
#"%OR%.default" <- function(x,y) {
#  mb <- max(max(maxBit(x),maxBit(y)),2)
#  if( is.numeric(x) && is.numeric(y) ) {
#    # Use C-Code
#    sapply(0:mb, function(i) (bit(x,i) | bit(y,i))*2^i )
#  } else {
#    do.call(paste,c(sep="",lapply(mb:0,function(i) ifelse((bit(x,i) | bit(y,i)),"1","0"))))
#  }
#}

# A bitwise NOT function
# x the number
# mb the number of bits to be processes
#bitwiseNOT <- function(x,mb=max(maxBit(x))) {
#  if( is.numeric(x) ) {
#                                        # Use C-Code
#    sapply(0:mb, function(i) (!bit(x,i))*2^i )
#  } else {
#    do.call(paste,c(sep="",lapply(mb:0,function(i) ifelse(!bit(x,i),"1","0"))))
#  }
#  
#}
                        
# A bitwise OR of many element
gsi.orSum <- function(...,g=2) {
  l <- list(...)
  if( length(l) > 1 )
    return(Recall(sapply(l,gsi.orSum)))
  x <- l[[1]]
  mb <- max(maxBit(x,g))
  do.call(paste,c(sep="",lapply(mb:0,function(i) ifelse(any(bit(x,i)),"1","0"))))
}
                       
# A bitwise AND of many element
#gsi.andSum <- function(...) {
#  l <- list(...)
#  if( length(l) > 1 )
#    return(andSum(sapply(l,andSum)))
#  x <- l[[1]]
#  mb <- max(maxBit(x))
#  do.call(paste,c(sep="",lapply(mb:0,function(i) ifelse(all(bit(x,i)),"1","0"))))
#}


# gives a list of bit setted in a binary number
# x the number (my be a vector)
# mb the maximum Bit to be used
whichBits <- function(x,mb=max(maxBit(x,g)),g=2,values=c(TRUE)) {
  if( length(x) > 1 )
    return(lapply(x,Recall,mb=mb))
  #unlist(lapply(0:mb,function(i) if(bit(x,i) %in% value) i else c()))
  (0:mb)[bit(x,0:mb) %in% values]
}

binary2logical <- function(x,mb=max(maxBit(x,g)),g=2,values=c(TRUE)) {
  if( length(x) > 1 )
    return(lapply(x,Recall,mb=mb))
  #unlist(lapply(0:mb,function(i) if(bit(x,i) %in% value) i else c()))
  bit(x,0:mb) %in% values
}

#takeIf <- function(c,x,y) {
#y <- ifelse(c,y,y)
#y[c]<-x
#y
#}
