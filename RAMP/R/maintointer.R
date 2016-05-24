maintointer <-
function(mainIndex){
###expand main effects to possible interaction effects
###output: a list of interaction effects
aa = outer(mainIndex, mainIndex, f <- function(x, y) {
				paste("X", x, "X", y, sep = "")
			})
aa[lower.tri(aa,diag=FALSE)] = NA
bb = as.vector(aa)
bb = bb[!is.na(bb)]
return(bb)
}
