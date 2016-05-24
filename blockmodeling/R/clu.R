"clu" <-
function(res,which=1,...){
	res$best[[which]]$clu
}

"partitions" <- 
function(res)lapply(res$best,function(x)x$clu)
