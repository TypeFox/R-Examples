msBP.rsample <-
function(n, msBPtree)
{
Rvec <- tree2vec(msBPtree$R)
Svec <- tree2vec(msBPtree$S)
res <- .C("rsample_msBP_C", as.integer(n), as.double(Rvec), as.double(Svec), 
		as.double(msBPtree$a), as.double(msBPtree$b), as.integer(msBPtree$max.s), 
		ans=as.double(rep(0, n)), PACKAGE = "msBP")
res$ans
}
