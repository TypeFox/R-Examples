msBP.marginalBeta <- function(y, maxS=10, output="tree")
{
veclen <- 2^(maxS+1) - 1
res <- .C("marginalBeta_C", out= as.double(rep(0, veclen)), as.double(y), as.integer(maxS), PACKAGE = "msBP")
if(output=="tree") vec2tree(res$out)
else res$out
}
