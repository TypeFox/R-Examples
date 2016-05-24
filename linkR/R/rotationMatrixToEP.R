rotationMatrixToEP <- function(R){
	# From: Computing Euler angles from a rotation matrix by Gregory G. Slabaugh

	if(R[3, 1] == -1 || R[3, 1] == -1){
		o1 <- 0
		if(R[3, 1] == -1){
			t1 <- pi/2
			p1 <- o1 + atan2(R[1, 2], R[1, 3])
		}else{
			t1 <- -pi/2
			p1 <- -o1 + atan2(-R[1, 2], -R[1, 3])
		}
		
		o2 <- 0
		t2 <- 0
		p2 <- 0
	}else{
		t1 <- -asin(R[3, 1])
		t2 <- pi - t1
		p1 <- atan2(R[3, 2]/cos(t1), R[3, 3]/cos(t1))
		p2 <- atan2(R[3, 2]/cos(t2), R[3, 3]/cos(t2))
		o1 <- atan2(R[2, 1]/cos(t1), R[1, 1]/cos(t1))
		o2 <- atan2(R[2, 1]/cos(t2), R[1, 1]/cos(t2))
	}

	angles <- list()
	angles[[1]] <- c(o1, t1, p1)
	angles[[2]] <- c(o2, t2, p2)

	angles
}