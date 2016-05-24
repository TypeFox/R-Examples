meshgridn = function(L){
  out = list()
  n = sapply(L, length)
  w = 1:length(L)
  for(i in w){
    out[[i]] = rep(rep(L[[i]], rep(prod(n[w<i]), n[i])), prod(n[w>i])) # prod(NULL) == 1, so this works.
  }
  return(out)
}
