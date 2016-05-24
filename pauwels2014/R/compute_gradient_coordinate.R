compute_gradient_coordinate <-
function(coord, theta, fun, step = 1e-6, ...){
	## Finite difference method to compute differential of fun
	## at theta with respect to the coord coordinate
	(fun( add_infinitesimal( theta, coord, c(), step ), ... ) - fun( add_infinitesimal( theta, c(), coord, step ), ... ) )/ (2 * step)
}
