skewdiv <-
function(xba, batch) {

  xba <- scale(xba)

  allpairs <- combn(length(levels(batch)), 2)

  sum(apply(allpairs, 2, function(y) (sum(batch==y[1]) + sum(batch==y[2]))*skewdivTwo(xba[batch==y[1],], xba[batch==y[2],]))/(nrow(xba)*(length(levels(batch))-1)))

}
