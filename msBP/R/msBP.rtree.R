msBP.rtree <-
function(a,b, max.s=10)
{
if(a<0) stop("a must be a positive real")
if(b<0) stop("b must be a positive real")
sizevec <- 2^(max.s+1)-1
res <- .C("randtree_C", as.double(a), as.double(b), as.integer(max.s), 
	S=as.double(rep(0,sizevec)), R=as.double(rep(0,sizevec)), PACKAGE = "msBP")
S <- vec2tree(res$S)
R <- vec2tree(res$R)
tree <-structure(list(S = S, R = R, a=a, b=b, max.s=max.s), class  = "msBPTree")
tree
}
