# Internal function:
# Calculate the total number of all common subsequences between 2 strings
# Subsequences cannot be interrupted by any item, i.e. q-w is not a common subsequence of q-w-e-r and q-e-w-r due to the interrupting 'e' in the second sequence
#
# param strA First string
# param strB Second string
# param sep Delimiter separating each items in a sequence
# param dropFirstItem Boolean. If true, the first item in each sequence is excluded from counting all subsequences
# param ignoreLenOneSubseq Boolean. If true, all length one subequences are not counted as common subsequences
# param ignoreLenZeroSubseq Boolean. If true, the length zero subsequence (empty set) is not counted as a common subsequence
# return The total number of all common subsequences as an integer

calACSstrStrict <- function(strA, strB, sep="-", dropFirstItem=FALSE, ignoreLenOneSubseq=FALSE, ignoreLenZeroSubseq=FALSE){
#   #test
#   strA <- 'aa-b-c'
#   strB <- 'aa-b-c-d'
#   sep='-'
#   dropFirstItem=FALSE
  if (ignoreLenZeroSubseq) {
    countACS <- 0 #count the empty set as one common subsequence
  } else {
    countACS <- 1
  }


  A <- unlist(strsplit(strA, sep))
  B <- unlist(strsplit(strB, sep))

  if (dropFirstItem){
    A <- A[-1]
    B <- B[-1]
  }

  listSubseq <- GeneratePossibleSubsequences(A, ignoreLenOneSubseq)

  for(i in 1:length(listSubseq)){
    if(is.subvector(subvec = listSubseq[[i]],vec = A) && is.subvector(subvec = listSubseq[[i]], vec = B)){
      countACS=countACS+1
    }
  }
  countACS
}


