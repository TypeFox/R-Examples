delta2a <-
function(delta){
      res = vector("numeric",length(delta)+1)
      count = 2
      res[1]  = 0
      while(count > 1 & count <= length(delta)+1){
         res[count] = exp(delta[count-1]) + res[(count-1)]
         count = count + 1
      }
      res
   }
