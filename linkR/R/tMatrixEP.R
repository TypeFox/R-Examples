tMatrixEP <- function(v, a){

	v <- uvector(v)

	t0 <- cos(a/2)
	t1 <- v[1]*sin(a/2)
	t2 <- v[2]*sin(a/2)
	t3 <- v[3]*sin(a/2)

	r <- matrix(0, 3, 3)
	r[1, ] <- c(2*(t0^2 + t1^2) - 1, 2*(t1*t2 - t0*t3), 2*(t1*t3 + t0*t2))
	r[2, ] <- c(2*(t1*t2 + t0*t3), 2*(t0^2 + t2^2) - 1, 2*(t2*t3 - t0*t1))
	r[3, ] <- c(2*(t1*t3 - t0*t2), 2*(t2*t3 + t0*t1), 2*(t0^2 + t3^2) - 1)

	r	
}
