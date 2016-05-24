med <-
function(x) {
y <- quantile(x, .5, type=3)
return(y[[1]])
}

