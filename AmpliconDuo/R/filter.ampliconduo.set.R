filter.ampliconduo.set <-
  function(x, min.freq = 1, OR = NULL, q = NULL, p = NULL, remove = FALSE){
    data.f <- lapply(X = x, FUN = filter.ampliconduo, min.freq, OR, q, p , remove = remove)
    return(data.f)
  }
