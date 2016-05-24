extract.list <-
function(L, s1, s2=NULL) {
  if(is.null(s2)) {
    sapply(L,function(x) with(x,get(s1)))
  } else {
    sapply(L,function(x) with(with(x,get(s1)),get(s2))  ) 
  }
}
