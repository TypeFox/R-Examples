ds_eqp_1 <- function(y, lambda = 1)
{
	.Call('dslice_ds_eqp_1', y, lambda, PACKAGE = 'dslice')
}

ds_1 <- function(y, lambda = 1, alpha = 1)
{
	.Call('dslice_ds_1', y, lambda, alpha, PACKAGE = 'dslice')
}

ds_eqp_k <- function(x, xdim, lambda = 1, slice = FALSE)
{
	if(slice) {
		.Call('dslice_dslice_eqp_k', x, xdim, lambda, PACKAGE = 'dslice')
	}else{
		.Call('dslice_ds_eqp_k', x, xdim, lambda, PACKAGE = 'dslice')
	}
}

ds_k <- function(x, xdim, lambda = 1, slice = FALSE)
{
	if(slice){
		.Call('dslice_dslice_k', x, xdim, lambda, PACKAGE = 'dslice')
	}else{
		.Call('dslice_ds_k', x, xdim, lambda, PACKAGE = 'dslice')
	}
}
