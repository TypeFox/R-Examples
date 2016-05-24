classIndex <- function(x){
  if(is.null(x)) return(1)
  if(is(x, 'logical'))return(2)
  if(is(x, 'integer'))return(3)
  if(is(x, 'numeric'))return(4)
  if(is(x, 'complex'))return(5)
  if(is(x, 'raw'))return(6)
  if(is(x, 'character'))return(7)
  8
}

index2class <- function(i, otherCharacter=TRUE){
##
## 1.  1 <= i <= 8
##
  if(i<1)stop('i = ', i, '; must be positive')
  if(i>8)stop('i = ', i, '; must be at most 8')
##
## 2.  integer?  
##
  if((i%%1)!=0)stop('i = ', i, '; must be an integer')
##
## 3.  other? 
##
  if(i==8){
    if(otherCharacter)return('character')
    return('other')
  }
## 
## 4.  1 <= i <= 7
##
  c('NULL', 'logical', 'integer', 'numeric', 
    'complex', 'raw', 'character')[i]
}

