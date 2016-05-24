Zstat <-
function(mat, p.null){ (pMME(mat) - p.null) * sqrt(sum(mat))/ sqrt(p.null*(1-p.null))
}
