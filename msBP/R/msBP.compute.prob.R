msBP.compute.prob <-
function(msBPtree, root=TRUE)
{
if(msBPtree$R$max.s != msBPtree$R$max.s) stop("S tree and R tree must have the same maximum scale")
R <- tree2vec(msBPtree$R)
S <- tree2vec(msBPtree$S)
max.s <- msBPtree$R$max.s
length.ans <- 2^(max.s+1) - 1
res <- .C("computeprob_C", as.double(S), as.double(R), as.integer(max.s), 
	as.double(msBPtree$a), as.double(msBPtree$b), ans=as.double(rep(0, length.ans)), root=as.integer(root), PACKAGE = "msBP")
vec2tree(res$ans)
}
