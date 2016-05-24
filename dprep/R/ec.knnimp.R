ec.knnimp <-
function(data,nomatr=0,k = 10) 
{
#xnom: vector containing the indexes of the nominal variables
#data: matrix containing data
  x=data
  N <- dim(x)[1]
  p <- dim(x)[2] 

#Checking if a row has a missing value
#print(sum(nas))
#if(sum(nas)==N) stop("Error: All cases have missing values. Cannot compute neighbors.")
rmiss=which(rowSums(is.na(data))!=0,arr.ind=T)
#print(rmiss)
#submatrix with complete rows
#matrix needed in case xcomplete has only one row
#print(x)
  xcomplete <- x[-rmiss,] 
#print(xcomplete)
  colnames(xcomplete)=seq(p)

#submatrix of rows with at least one missing value
  xbad <- x[rmiss,,drop=FALSE ]
#print(xbad)
#forming logical vector of nominal variables
 xnom=seq(p) %in% nomatr
#print(xnom)
#Locating the missing values in the missing submatrix 
  xnas <- is.na(xbad)
  xbadhat <- xbad
#print(xcomplete)
#print(xbadhat)
  for(i in seq(nrow(xbad))) 
  {
    xinas <- xnas[i,  ]
    xbadhat[i,  ] <- nnmiss(xcomplete, xbad[i,  ],xinas,xnom, K = k)
  }
#print(xbadhat)
  x[rmiss,  ] <- xbadhat
  data2 <-x
  return(data2)
}
