gflux <-
function(ct, c0 = NULL, T, V, A = 1, M = 44, t = 1/60, p = 101325){
	R <- 8.314
	dc <- ifelse(is.null(c0), ct, ct-c0)
	flux <- (dc * V * M * p) / (t * R * (T + 273.15) * A)
	flux
}