#'LastOccur
#'
#'Returns index of last occurence: Each row is scanned for the last column in which a certain integer is shown. 
#'Consider the following example: one row has the values 1-2-2-1-0-4, last occurence of 0 would be the fith column,
#'last occurence of 1 would be the forth column, and so on.  
#'
#'@param x Dataframe or matix containing one sequence per row
#'@param y The value of interest
#'
#'@return returns a vector containing the index of the last 
#'event occurence for every row. 
#'
#'
#'@examples
#'# Example 1: Small artificial data
#'
#'  my.data<-matrix(c(1,0,1,1,
#'                    0,0,1,0,
#'                    1,0,0,0,
#'                    0,0,0,1),4,4, TRUE) # create data
#'
#'  my.data # inspect sampe data
#'  LastOccur(my.data,1) # last Occurence of one
#'  LastOccur(my.data,0) # last Occurence of zero
#'
#'
#'# Example 2: Real data
#'
#'  data(CouplesCope)
#'  LastOccur(CouplesCope[,2:49],1)
#'
#'@export


LastOccur<-function(x,y){
  if(!(is.matrix(x)|is.data.frame(x))) warning("x must be a matrix or dataframe!")

  output<-numeric(length(x[,1]))

  for(k in 1:length(x[,1])){
    for(i in 1:length(x[1,])){
      if(x[k,i]==y) {output[k]<-i}
    }
  }
  return(output)
}


  

