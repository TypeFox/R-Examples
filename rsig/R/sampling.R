# cv10f = function(n, f = 10L) {
#   groups = (sample(n) - 1L) %% f + 1L
#   split(1:n, groups)
# }

balanced.subset.sample = function(n=10L, group, ratio=.632) {
  ntotal = length(group)
  stopifnot(ntotal > 0L)            
  stopifnot(ratio >= 0L && ratio <= 1L)
  
  f = unique(group)

  trainidx=NULL
  testidx=NULL
  subsets <- lapply(seq_len(n), function(r) {
    for(i in f) {
      idx = which(group==i)
      n = length(idx) 
      train.sub = (sample(n) <= n*ratio)
      trainidx = union(trainidx, idx[train.sub])
      testidx = union(testidx, idx[!train.sub])
    }
    
    return(list(train.idx=trainidx, test.idx=testidx))
  })
  
  return(subsets)
}
