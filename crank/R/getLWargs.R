# Argument:
# matrix of ranks, possibly with missing values
# Return value:
# ranksums - if x contains no missing values, NULL
# OR if x does contain missing values, a list containing:
# the column sums of all complete rows ("B" rows)
# Arows - the index of the row with the smallest number of missing values OR
# a matrix of the permutations of the row indices of the rows with the
# smallest number of missing values
# nArows - the number of permutations of the Arows

getLWargs<-function(x) {
 dimx<-dim(x)
 # get the default "no Brows" ranksums here
 ranksums<-rep(1,dimx[2])
 matholes<-is.na(x)
 if(sum(matholes)) {
  rowholes<-rowSums(matholes)
  # get the Brow indices
  Brows<-which(rowholes == 0)
  if(length(Brows) > 1) ranksums<-colSums(x[Brows,])
  if(length(Brows) == 1) ranksums<-x[Brows,]
  # row(s) with smallest number of NAs
  minholes<-min(rowholes[rowholes > 0])
  # indices of above
  Arows<-which(rowholes == minholes)
  # if it's not a single number
  if(length(Arows) > 1) {
   nArows<-gamma(length(Arows)+1)
   # it's a matrix of permutations
   Arows<-permute(Arows)
  }
  else nArows<-1
 }
 # if no NAs, there's nothing to do
 else return(NULL)
 return(list(ranksums=ranksums,Arows=Arows,nArows=nArows,Brows=Brows))
}
