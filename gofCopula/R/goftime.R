.get.time = function(x){
#   x.year = floor(x/31104000)
#   x.remainder = x%%31104000
#   x.month = floor(x.remainder/2592000)
#   x.remainder = x.remainder%%2592000
  x.day = floor(x/86400)
  x.remainder = x%%86400
  x.hour = floor(x.remainder/3600)
  x.remainder = x.remainder%%3600
  x.min = floor(x.remainder/60)
  x.remainder = x.remainder%%60
  out = list(x.day, x.hour, x.min, x.remainder)
  class(out) = "goftime"
  out
}

print.goftime = function(x, ...){
 # print(paste("The computation will take approximately ", x[[1]], " d ", x[[2]],":",x[[3]],":",x[[4]], ".", sep = ""))
    print(sprintf("The computation will take approximately %d d, %d h, %d min and %d sec.", x[[1]], x[[2]], x[[3]], x[[4]]))
}