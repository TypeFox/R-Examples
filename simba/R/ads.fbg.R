## function to find the best gradient (fbg())
## whether ads or ads.hot is to be used has to be changed by hand..
"ads.fbg" <- function(nspec, nplots, grad.v, n.iter=100, method="ads", ...){

## internal function for producing a species matrix
## and testing it against the given gradient
## which is derived from a dca with defaults in vegan
"fbg" <- function(nspec, nplots, grad.v, method=method, ...){
	if (!is.na(pmatch(method, "ads"))) 
        method <- "ads"
	METHODS <- c("ads", "ads.hot", "makead")
	method <- pmatch(method, METHODS)
	if (method == 1){
		mat.tmp <- ads(nspec, nplots, grad.v=grad.v, ...)
	}
	else if (method == 2){
		mat.tmp <- ads.hot(nspec, nplots, grad.v=grad.v, ...)
	}
	else if (method == 3){
		mat.tmp <- makead(nspec, nplots, grad.v=grad.v, ...)
	}
	tmp.dca <- decorana(mat.tmp)
	r2.adj <- summary(lm(scores(tmp.dca, display="sites")[,1] ~ I(grad.v/max(grad.v))))	$adj.r.squared
	res <- list(mat.tmp = mat.tmp, r2.adj = r2.adj)
	return(res)
}

## making a result list upon that the function will be applied
best.mat <- vector("list", n.iter)
## applying the function to the list
tmp <- lapply(best.mat, function(x) fbg(nspec, nplots, grad.v, method=method, ...))
## getting the best matrix
mat <- tmp[[which.max(sapply(tmp, function(x) x$r2.adj))]]$mat.tmp
## getting the r2.adj of the best matrix
r2.adj <- tmp[[which.max(sapply(tmp, function(x) x$r2.adj))]]$r2.adj
## putting them together
res <- list(mat = mat, r2.adj = r2.adj)
## return result
return(res)
}