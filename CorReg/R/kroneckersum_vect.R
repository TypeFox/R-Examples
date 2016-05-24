# ' Kronecker sum of two vectors
# ' @param A vector of length n
# ' @param B vector of length m
# ' 
kroneckersum_vect<-function(A,B){
  return(kronecker(A,rep(1,times=length(B)))+kronecker(rep(1,times=length(A)),B))
}