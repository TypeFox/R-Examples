# Internal function:
# Genearte all possible subsequences of the input vector. Require all items in the vector to be unique
#
# param vin input vector that contains unique items
# param ignoreLenOneSubseq Boolean. If true, all length one subequences are not counted as common subsequences
# return a list of vectors, each a unique subsequence of the input vector

GeneratePossibleSubsequences <- function(vin, ignoreLenOneSubseq=FALSE){
  #test
  #vin=A
  #i=2
  #j=2

  lout<-list()

  #i is the length of the subsequence extracted
  #j is the starting position of the subsequence in vin

  if (ignoreLenOneSubseq) {
    for(i in 2:length(vin)){
      for(j in 1:(length(vin)-(i-1))){
        lout[[length(lout)+1]] <- vin[j:(j+i-1)]
      }
    }
  } else {
    for(i in 1:length(vin)){
      for(j in 1:(length(vin)-(i-1))){
        lout[[length(lout)+1]] <- vin[j:(j+i-1)]
      }
    }
  }
  lout
}
