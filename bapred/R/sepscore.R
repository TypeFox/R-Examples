sepscore <-
function(xba, batch, k=10) {

  allpairs <- combn(length(levels(batch)), 2)

  sum(apply(allpairs, 2, function(y) (sum(batch==y[1]) + sum(batch==y[2]))*sepscoreTwo(xba[batch==y[1],], 
    xba[batch==y[2],], k=k))/(nrow(xba)*(length(levels(batch))-1)))

}
