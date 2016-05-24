# Arguments:
# x - a matrix of ranks
# ranksums - column sums of the current block of "B" (complete)
#  rows that are used to direct the imputation
# Arow - the row of x to be completed
# maxcon - whether to complete the row maximally consistent (TRUE) or
#  maximally inconsistent (FALSE) with the ranksums
# Return value:
# If there are no tied ranksums, the matrix x with imputed values for row Arow
# If there are tied ranksums, a list of matrices with one value filled in each
# row Arow

fillArow<-function(x,ranksums=NA,Arow,maxcon=TRUE) {
 # all possible ranks
 allranks<-1:dim(x)[2]
 # ranks missing in this row
 misranks<-allranks[!allranks%in%x[Arow,]]
 # how many NAs are in this row
 rowholes<-is.na(x[Arow,])
 cranksums<-ranksums[rowholes]
 if(length(cranksums) > length(unique(cranksums))) {
  # find a tie
  scranksums<-sort(cranksums)
  cindex<-1
  while(scranksums[cindex]!=scranksums[cindex+1] &&
   cindex < length(scranksums)) cindex<-cindex+1
  # logical where TRUE is the conjunction of the tied ranksum(s) and NAs
  thismiss<-(ranksums == scranksums[cindex]) & rowholes
  # there is at least one tie, calculate the number of permutations needed
  # and create a list with that length
  nties<-sum(cranksums == cranksums[cindex])
  nperms<-gamma(nties+1)
  newx<-vector("list",nperms)
  fillorders<-permute(1:nties)
  thismisranks<-misranks[nties]
  for(perm in 1:nperms) {
   newx[[perm]]<-x
   # fill in the ranks for this set of tied ranksums
   newx[[perm]][Arow,thismiss]<-thismisranks[fillorders[perm,]]
  }
  x<-newx
 }
 else {
  # there are no ties, just fill and return
  fillorder<-order(ranksums[rowholes])
  if(!maxcon) fillorder<-rev(fillorder)
  x[Arow,rowholes]<-misranks[fillorder]
 }
 return(x)
}


# Arguments:
# x - a matrix of ranks
# maxcon - whether to fill ranks maximally (TRUE) or minimally
# consistent with the filled rows
# Return value:
# If all of the rows listed in Arows have already been filled,
# x is returned. This prevents the function from filling rows
# in matrices that are not in the current block of rows.
# x with one row filled if there is only one row with the
# smallest number of missing values.
# Otherwise, a list of copies of x, the length of which is
# equal to the number of permutations of row indices.
# each copy will have the row filled indicated by the first
# index in that permutation.
# If there are also ties in the ranksums for missing column
# values, each matrix will be permuted again by fillArow

fillArows<-function(x,maxcon=TRUE) {
 LWargs<-getLWargs(x)
 if(!is.null(LWargs)) {
  if(is.matrix(LWargs$Arows)) {
   # create the same number of matrices as rows in Arows
   # and send them off to fillArow to fill the first row
   # of each. Forget about the other rows, they'll happen
   newx<-vector("list",LWargs$nArows)
   for(newmat in 1:LWargs$nArows)
    newx[[newmat]]<-fillArow(x,LWargs$ranksums,LWargs$Arows[newmat,1],maxcon)
   x<-newx
  }
  else x<-fillArow(x,LWargs$ranksums,LWargs$Arows,maxcon)
 }
 return(x)
}
