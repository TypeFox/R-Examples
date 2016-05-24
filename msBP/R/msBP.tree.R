msBP.tree <-
function(max.s=10)
{
sizevec <- 2^(max.s+1)-1
S <- vec2tree(rep(0,sizevec))
R <- vec2tree(rep(0,sizevec))
tree <-structure(list(S = S, R = R, a=NULL, b=NULL, max.s=max.s), class  = "msBPTree")
tree
}
