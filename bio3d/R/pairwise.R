"pairwise" <-
function(N) {
  
  # deternine the indices for  pairwise
  # comparison of N things (where N = 'number.of.things')
  # Used in function 'identity'
  #
  #    A <- read.fasta(system.file("examples/kinesin.fa",
  #                                 package = "bio3d"))
  #    
  #    N<-nrow(A$ali)
  #    # comparison indices
  #    inds<-pairwise(N)
  #    # store score in 's'
  #    s<-rep(NA,nrow(inds))
  #    # go through all comparisons
  #     for(i in 1:nrow(inds)) {
  #       s[i]<-ide(A$ALI[ inds[i,1], ], A[ inds[i,2], ])
  #     }
  #     # reformat 's' as a NxN matrix
  #     mat<-matrix(NA,ncol=N,nrow=N)
  #     mat[inds]<-s; mat[inds[,c(2,1)]]<-s

  if (!is.numeric(N) || N<0 || length(N)>1) {
    stop("pairwise: N must be positve numeric and of length 1")
  }

  pair <- matrix(NA, ncol=2,
                 nrow= ( ((N*N)/2)-(N/2) )) 
  start<-1
  for(i in 1:(N-1)) {
    end  = N-i
    inds = start:( (start-1)+end )
    pair[inds,1] = rep(i, end)
    pair[inds,2] = (i+1):N
    start = start+end
  }
  return(pair)
}

