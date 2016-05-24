`computeY` <-
function(g,x) apply(cbind(x), 2, function(z) tapply(z,g,mean))

