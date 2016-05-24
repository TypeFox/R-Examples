##= count
count =
# counts the number of elements in x equal to elm; e.g.
# > x = c(1,1,2,2,2,3,3)
# > y = count(x,c(1,2))
# > y
# [1] 2  3
function(x, elm=sort(unique(x))){
  j <- 0;  res <- 0
  if(any(is.na(x)))stop("data contain missing values")
  unique.elm <- unique(elm)
  for (j in 1:length(unique.elm))res[j] = length(x[x==unique.elm[j]])
  return(res)
}
