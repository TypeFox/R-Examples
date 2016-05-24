getleafnormal <- function(l){
	if(class(l) == "leaffile")
		m <- l$XYZ
	else
		m <- l
	vec1 <- makevector(m[1,], m[2,])
	vec2 <- makevector(m[1,], m[(nrow(m)-1),])
	N <- xprod(vec1, vec2)
	N
}