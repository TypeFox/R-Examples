nn.ent <-function (th, k = 4){

th <- as.matrix(th)
p <- ncol(th)
n <- nrow(th)

if (p == 1) {
	t <- .C("nnone", as.double(t(th)), as.integer(n), 
        as.integer(k), D = as.double(0), PACKAGE = "abctools")
}
else {
	t <- .C("nnk", as.double(t(th)), as.integer(n), 
	as.integer(p), as.integer(k), D = as.double(0), 
	PACKAGE = "abctools")
}

ent <- t$D

ent
}

