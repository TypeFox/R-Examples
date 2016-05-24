# Internal function:
# Calculate the total number of all common subsequences between 2 strings.
# Subsequences can be interrupted by items, i.e. q-w is a common subsequence of q-w-e-r and q-e-w-r
#
# param strA First string
# param strB Second string
# param sep Delimiter separating each items in a sequence
# param dropFirstItem Boolean. If true, the first item in each sequence is excluded from counting all subsequences
# return The total number of all common subsequences as an integer


calACSstrLoose <- function(strA, strB, sep="-", dropFirstItem=FALSE){
#   #test
#   strA <- 'aa-b-c'
#   strB <- 'aa-b-c-d'
#   sep='-'
#   dropFirstItem=FALSE

  A <- unlist(strsplit(strA, sep))
  B <- unlist(strsplit(strB, sep))

  if (dropFirstItem){
    A <- A[-1]
    B <- B[-1]
  }

    N <- array(data=1,dim=c(length(A)+1,length(B)+1))
    for (i in 2:(length(A)+1)){
      for (j in 2:(length(B)+1)){
        ifelse( identical(A[i],B[j]), N[i,j] <- N[i-1,j-1]*2, N[i,j] <- N[i-1,j]+N[i,j-1]-N[i-1,j-1])
      }
    }
    N[i,j]
}


