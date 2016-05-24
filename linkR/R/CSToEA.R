CSToEA <- function(cs1, cs2){
	
	# FIND ROTATION MATRIX BETWEEN COORDINATE SYSTEMS
	SVD <- svd(t(cs1) %*% cs2)
	L <- diag(SVD$d)
	S <- ifelse(L<0, -1, L)
	S <- ifelse(L>0, 1, L)
	RM <- SVD$v %*% S %*% t(SVD$u)

	# FINDS THE EULER ANGLES FROM ROTATION MATRIX
	if(RM[3, 1] == -1 || RM[3, 1] == -1){
		o1 <- 0
		if(RM[3, 1] == -1){
			t1 <- pi/2
			p1 <- o1 + atan2(RM[1, 2], RM[1, 3])
		}else{
			t1 <- -pi/2
			p1 <- -o1 + atan2(-RM[1, 2], -RM[1, 3])
		}
		
		o2 <- 0
		t2 <- 0
		p2 <- 0
	}else{
		t1 <- -asin(RM[3, 1])
		t2 <- pi - t1
		p1 <- atan2(RM[3, 2]/cos(t1), RM[3, 3]/cos(t1))
		p2 <- atan2(RM[3, 2]/cos(t2), RM[3, 3]/cos(t2))
		o1 <- atan2(RM[2, 1]/cos(t1), RM[1, 1]/cos(t1))
		o2 <- atan2(RM[2, 1]/cos(t2), RM[1, 1]/cos(t2))
	}

	list(c(o1, t1, p1), c(o2, t2, p2))
}