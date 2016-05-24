`gq.bounds` <- function(numrows.pt,numcols.pt,NNtots)
{
# Calculate a matrix of bounds
#  Each row of NNbounds.0 is of the form min of [1, 1], min of [1, 2],
#     . . . min of [1, C], min of [1, 1], min of [1, 2], . . .  . . .
#     min of [R, C]], max of [1, 1], max of [1, 2], . . . max of [1, C],
#     max of [2, 1], max of [2, 2], . . .  . . . max of [R, C]
  num.prec <- nrow(NNtots)
  dim.NNtots <- numcols.pt+numrows.pt
  dim.NNs    <- numcols.pt*numrows.pt
  NNbounds <- matrix(0,nrow=num.prec,ncol=2*dim.NNs)
  for (ii in 1:num.prec){
    tempvec <- NNtots[ii,(numrows.pt+1):dim.NNtots]  #tempvec = precinct ii col sums
    for(jj in 1:numrows.pt){
      for(kk in 1:numcols.pt){
        #  Set the min & max:
        NNbounds[ii,((jj-1)*numcols.pt)+kk] <- max(0,NNtots[ii,jj]-sum(tempvec[-kk]))
        NNbounds[ii,((jj-1)*numcols.pt)+kk+dim.NNs] <- min(NNtots[ii, jj], tempvec[kk])
      }
    }
  }
  return(NNbounds)
}
