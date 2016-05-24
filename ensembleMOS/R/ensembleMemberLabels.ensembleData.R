`ensembleMemberLabels.ensembleData` <-
function (x) 
{ 
 k <- attr(x, "ensembleSize")
 (dimnames(x)[[2]])[1:k]
}

