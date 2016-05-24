BOG <-
function(data=NULL, data.type=c("data","pval"), cog.file = NULL, hg.thresh=0.05, 
gsea=FALSE, DIME.K=5, DIME.iter=50, DIME.rep=5)
{
	set.seed(12345)
	data_type=match.arg(data.type)
		
	result=BOGest(data, data.type, cog.file, hg.thresh=hg.thresh, 
	gsea=gsea, DIME.K, DIME.iter, DIME.rep)
	result$call <- match.call()
	class(result) <- "BOG"
	result
}
