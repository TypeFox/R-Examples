var.correct <-
function (x, incby=1, auto=FALSE) {
if (auto==FALSE) { 
y <- x+incby   
}
else {   
xmin <- abs(min(x))
y <- x+xmin+incby
}
return(y)
}
