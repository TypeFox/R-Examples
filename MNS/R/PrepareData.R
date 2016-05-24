PrepareData <-
function(dat){
  # we prepare data to be put in the format we would expect in the PenalizedLMM function
  # 
  # INPUT:
  #      - dat: a list where each entry is the data for a given subject
  #
  # OUTPUT:
  #      - o: a list where each entry contains 3 elements. The first is the concatenated time series for the ith node across all subjects
  #           and the second in the concatenated design matrix across all subjects (ie time series of remaining nodes)
  #           the final element is the design for the random effects (essentially block diagonals of 2nd element)
  #           UPDATE: o now only contains 2 elements as we improved random effects implementation!
  
  p = ncol(dat[[1]])
  n = nrow(dat[[1]])
  N = length(dat)
  o = lapply(vector("list", p), FUN=function(x){ list(rep(0, N*n), matrix(0, ncol=p-1, nrow=N*n))}) # prepare list to be fill below
  
  for (i in 1:p){
    # fill in response:
    o[[i]][[1]] = unlist(lapply(dat, FUN=function(x){x[,i]}))
    
    # fill in design for fixed and random effects:
    for (j in 1:N){
      o[[i]][[2]][((j-1)*n+1):(j*n), ] = dat[[j]][,-i]
    }    
  }
  return(o)
}
