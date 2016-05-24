# msplit <- function( vec, m ) {
#   V <- length(vec)
#   SV <- ceiling(V/m)
#
#   ranges <- list()
#   cur <- 1
#   for( i in 1:m ) {
#     if( cur<=V ) {
#       curEnd <- cur+SV-1
#       if( curEnd>V )
#         curEnd <- V
#       ranges[[i]] <- vec[ cur:curEnd ]
#       cur <- cur + SV
#     }
#   }
#   return( ranges )
# }

msplit <- function(vec, m){
  N <- length(vec)
  G <- ceiling(N/m)
  numSmall <- G*m - N
  numRegular <- m - numSmall

  curIndex <- 1
  ranges <- list()
  for(i in 1:m){
    curSize <- G-1
    if(numRegular > 0){
      numRegular <- numRegular - 1
      curSize <- G
    }

    curEnd <- curIndex + curSize - 1
    if(curEnd > N){
      warning("Bug in msplit(), curEnd > N; fixed, but still a bug.")
      curEnd <- N
    }

    ranges[[i]] <- vec[curIndex:curEnd]
    curIndex <- curEnd + 1
  }

  return(ranges)
}

## DEBUGGING ONLY
#for(i in 1:11)
#  print(msplit(1:(100-i), 10))
#print(msplit(1:45,12))
#print(msplit(13:45,12))
#print(msplit(1:760,8))

