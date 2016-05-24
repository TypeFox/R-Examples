avectors <- function(u, v){
	acos(sum(u*v) / (sqrt(sum(u*u)) * sqrt(sum(v*v))))
}