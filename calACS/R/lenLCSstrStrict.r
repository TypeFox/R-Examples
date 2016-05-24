# Internal function:
# Calculate the length of the longest common subsequence between 2 input strings
# param strA First string
# param strB Second string
# param sep Delimiter separating each items in a sequence
# param dropFirstItem Boolean. If true, the first item in each sequence is excluded from counting all subsequences
# return a vector containing the length of each common subsequence


lenLCSstrStrict <- function(strA, strB, sep="-", dropFirstItem=FALSE){
  #     #test
  #     strA <- 'aa-b-c'
  #     strB <- 'aa-b-c-d'
  #     sep='-'
  #     dropFirstItem=FALSE
  #     i=4
  #     i=5
  #     i=6

  lenACS <- NULL #count the empty set as one common subsequence

  A <- unlist(strsplit(strA, sep))
  B <- unlist(strsplit(strB, sep))

  if (dropFirstItem){
    A <- A[-1]
    B <- B[-1]
  }

  listSubseq <- GeneratePossibleSubsequences(A)

  for(i in 1:length(listSubseq)){
    if(is.subvector(subvec = listSubseq[[i]],vec = A) && is.subvector(subvec = listSubseq[[i]], vec = B)){
      lenACS[length(lenACS)+1] <- length(listSubseq[[i]])
    }
  }
#
#   if (is.na(lenACS) || is.null(lenACS)){
#     0
#   }
#   else{
#     lenACS[length(lenACS)+1] <- 0
#     max(lenACS)
#   }

  if(is.null(lenACS)){
    lenACS <- 0
  }
  max(lenACS)
}

#lenACS <- c(1,2)
#lenACS <- NULL
#max(lenACS, 0)

