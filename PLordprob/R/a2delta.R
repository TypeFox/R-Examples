a2delta <-
function(a){
      res = vector("numeric",length(a)-1)
      count = 2
      while(count > 1 & count <= length(a)){
         res[count-1] = log(a[count] - a[(count-1)])
         count = count + 1
      }
      res
   }
