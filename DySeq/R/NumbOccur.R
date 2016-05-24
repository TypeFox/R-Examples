#'NumbOccur
#'
#'Returns the Number of occurences for sequences with different lenghts.
#'Outside the context of sequence analysis this means that for each frequency of one specifics integer will be computed. 
#'
#'@param x Dataframe or matix containing one sequence per row
#'@param y single integer: represents the occurence that should be counted
#'@param t optional vector that contains the lengths of sequences
#'@param prop if TRUE: proportion will be computed,if FALSE: sum will be computed
#'
#'
#'@return returns a vector containing containing the number of occurences.
#'
#'
#'
#'@examples
#'# Example 1: Small artificial example
#'
#'# Creating data, if sequence ends, rest should be 'NA'
#'seq1<-c(1,0,0,0,1,0,1, NA, NA, NA) # 3 out of 7 Entrys should be round about .43
#'seq2<-c(1,1,1,1, NA, NA, NA, NA, NA, NA) # 4 out of 4 should be 1
#'seq3<-c(1,0,0,0,1,1, NA, NA, NA, NA) # 3 out of 6 should be .50
#'my.data<-rbind(seq1,seq2,seq3)
#'
#'# Determine the proportion of ones in my.data
#'NumbOccur(my.data,1)
#'NumbOccur(my.data,1, prop=FALSE) # compute absolute frequencies
#'
#'
#'# Example 2: Real data dyadic sequences
#'# A researcher is interested in how often was a certain behavior
#'# shown till another one stopped completely
#'
#'my.last<-LastOccur(CouplesCope[,2:49],1) # how long till stress ended?
#'NumbOccur(CouplesCope[,50:97],1, my.last) # how often did dyadic coping occure in this time?
#'
#'@export


NumbOccur<-function(x,y,t=NA, prop=TRUE){
  if(!(is.matrix(x)||is.data.frame(x))) warning("x must be a matrix or dataframe!")
  if(!(is.numeric(t)||is.na(t))) warning("t must be NA or a vector with number elements equal number sequences")
  if(!is.logical(prop)) warning("Argument prob must be logical")

  output<-numeric(length(x[,1]))
  index<-numeric(length(x[,1]))


  if(is.na(t[1])){
    index<-rep(length(x[1,]), length(x[,1]))
    for(k in 1:length(x[,1])){
      output[k] <- sum(x[k,])
    }
  }

  if(is.numeric(t)){index<-t}

  for(k in 1:length(x[,1])){
   output[k]<-sum(as.numeric(x[k,1:index[k]]), na.rm=TRUE)
  }
  if(prop)output<-output/index

  return(output)
}


