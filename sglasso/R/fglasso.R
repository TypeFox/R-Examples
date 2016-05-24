fglasso <- function(S, model, tp, p, ...){
	this.call <- match.call()
	if(!is.element(class(model), c("list", "matrix")))
	stop("model must be a list or a matrix")
	if(is.list(model)){
		check <- unlist(lapply(model, function(x) length(x) != 2))
		if(any(check))	stop("wrong model specification")
		check <- unlist(lapply(model, function(x) !is.element(x, c("c", "u", "t", "ut", "."))))
		if(any(check))	stop("wrong model specification")
		model <- matrix(unlist(model), length(model), 2, byrow = TRUE)
	}
	if(is.matrix(model)){
		if(dim(model)[1] > tp)
		stop("dim(model)[1] > tp")
		if(dim(model)[2] != 2)
		stop("dim(model)[2] != 2")
	}
	if(model[1, 1] == ".")
	stop("model[1, 1] is equal to zero")
	if(dim(S)[1] != tp * p | dim(S)[2] != tp * p)
	stop("dim(S) != dim(mask)")
	mask <- fglasso_model2mask(model, tp, p)
	out_fglasso <- sglasso(S = S, mask = mask, ...)
	out_fglasso$call <- this.call
	out_fglasso$model <- model
	out_fglasso$p <- p
	out_fglasso$tp <- tp
	class(out_fglasso) <- c("fglasso", class(out_fglasso))
	out_fglasso
}
