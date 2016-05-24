calcDistance <-
function(x, y, type="", method="hamming"){
if(missing(x) || missing(y))
stop("x and/or y is missing.")

if(tolower(method) == "hamming"){
ret <- sum(as.logical(unlist(x)) != as.logical(unlist(y)))

if(tolower(type) == "adjmatrix") #Divide by 2 because of duplicates
ret <- ret / 2
}else{
stop(sprintf("%s is unknown.", as.character(method)))
}

return(ret)
}
