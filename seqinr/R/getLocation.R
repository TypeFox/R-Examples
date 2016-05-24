#
# To get the location of subsequences from an ACNUC server
#

getLocation <- function(object, ...) UseMethod("getLocation")

getLocation.list <- function(object, ...)
  lapply(seq_len(length(object)), function(i) getLocation(object[[i]], ...))

getLocation.default <- function(object, ...)
  stop(paste("no getLocation method for objects of class:", class(object)))

getLocation.SeqAcnucWeb <- function(object, ..., socket = autosocket()){
  getLocationSocket <- function( socket, name){

   writeLines(paste("isenum&name=",name,sep=""),socket,sep="\n")
         res = readLines( socket , n=1 )
         number = parser.socket(res)[1] 
   
         writeLines(paste("readsub&num=",number,sep=""),socket,sep="\n")
         res2 = readLines( socket , n=1 ) 
         rr = parser.socket(res2)
  
   # Test si subsequence           
     
   l=list() 
         if(as.numeric(rr[5]) != 0){
     warning("It's a parent sequence\n")
     return( NA )
    }
         else {
    i=1
    writeLines(paste("readext&num=",rr[6],sep=""),socket,sep="\n")    
    res3 = readLines( socket , n=1 )
    r = parser.socket(res3)
    l[[i]] = as.numeric(c(r[3],r[4]))
    n=r[5] 
  }
        while(as.numeric(n) != 0){
    i=i+1
    writeLines(paste("readext&num=",n,sep=""),socket,sep="\n")   
    res4 = readLines( socket , n=1 )
    rrr = parser.socket(res4)
    l[[i]] = as.numeric(c(rrr[3],rrr[4]))
    n=rrr[5]
      }
  return(l)
} 

  unlist(getLocationSocket(socket, name = object))
}

getLocation.qaw <- function(object, ...) getLocation(object$req, ...)

getLocation.logical <- function (object, ...)
  object # so that NA is returned for virtual lists
