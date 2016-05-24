
xshewhartrunsrules.matrix <- function(mu, c=1, type="12") {
	# Shewhart chart
	if (type=="1") {
  		p0  <- pnorm(  3*c, mean=mu ) - pnorm( -3*c, mean=mu)
  		Q <- p0
	}
	
	# 2 of 3 beyond +-2 sigma
	if (type=="12") {
		dimQ <- 7
		pl <- pnorm( -2*c, mean=mu ) - pnorm( -3*c, mean=mu)
		p0 <- pnorm(  2*c, mean=mu ) - pnorm( -2*c, mean=mu)
		pr <- pnorm(  3*c, mean=mu ) - pnorm(  2*c, mean=mu)
		
#                      1    2    3    4    5    6    7
#		       0000 1000 0100 0010 0001 1001 0110
#		1 0000 p0   pl   0    pr   0    0    0
#		2 1000 0    0    p0   0    0    0    pr
#               3 0100 p0   0    0    pr   0    0    0
#               4 0010 0    0    0    0    p0   pl   0
#               5 0001 p0   pl   0    0    0    0    0
#               6 1001 0    0    p0   0    0    0    0
#               7 0110 0    0    0    0    p0   0    0

		Q <- diag(0,dimQ)
		
		Q[1,2] <- pl;  Q[1,1] <- p0;  Q[1,4] <- pr
		               Q[2,3] <- p0;  Q[2,7] <- pr
			       Q[3,1] <- p0;  Q[3,4] <- pr
	        Q[4,6] <- pl;  Q[4,5] <- p0
		Q[5,2] <- pl;  Q[5,1] <- p0
		               Q[6,3] <- p0
		               Q[7,5] <- p0
	}	
	
	# 4 of 5 beyond +-1 sigma
	if (type=="13") {
		dimQ <- 29
		pl <- pnorm(  -c, mean=mu ) - pnorm( -3*c, mean=mu)
		p0 <- pnorm(   c, mean=mu ) - pnorm(   -c, mean=mu)
		pr <- pnorm( 3*c, mean=mu ) - pnorm(    c, mean=mu)
		
		Q <- diag(0,dimQ)
		
		Q[ 1, 2] <- pl;  Q[ 1, 1] <- p0;  Q[ 1,11] <- pr
		Q[ 2, 4] <- pl;  Q[ 2, 3] <- p0;  Q[ 2,12] <- pr
		Q[ 3, 5] <- pl;  Q[ 3, 1] <- p0;  Q[ 3,11] <- pr
		Q[ 4, 7] <- pl;  Q[ 4, 6] <- p0;  Q[ 4,13] <- pr
		Q[ 5, 8] <- pl;  Q[ 5, 3] <- p0;  Q[ 5,12] <- pr
		Q[ 6, 9] <- pl;  Q[ 6, 1] <- p0;  Q[ 6,11] <- pr
		                 Q[ 7,10] <- p0;  Q[ 7,14] <- pr
		                 Q[ 8, 6] <- p0;  Q[ 8,13] <- pr
		                 Q[ 9, 3] <- p0;  Q[ 9,12] <- pr
		                 Q[10, 1] <- p0;  Q[10,11] <- pr
		Q[11,16] <- pl;  Q[11,15] <- p0;  Q[11,19] <- pr
		Q[12,17] <- pl;  Q[12,15] <- p0;  Q[12,19] <- pr
		Q[13,18] <- pl;  Q[13,15] <- p0;  Q[13,19] <- pr
		                 Q[14,15] <- p0;  Q[14,19] <- pr
		Q[15, 2] <- pl;  Q[15, 1] <- p0;  Q[15,20] <- pr
		Q[16, 4] <- pl;  Q[16, 3] <- p0;  Q[16,21] <- pr
		Q[17, 8] <- pl;  Q[17, 3] <- p0;  Q[17,21] <- pr
		                 Q[18, 3] <- p0;  Q[18,21] <- pr
		Q[19,23] <- pl;  Q[19,22] <- p0;  Q[19,24] <- pr
		Q[20,16] <- pl;  Q[20,15] <- p0;  Q[20,25] <- pr
		Q[21,17] <- pl;  Q[21,15] <- p0;  Q[21,25] <- pr
		Q[22, 2] <- pl;  Q[22, 1] <- p0;  Q[22,26] <- pr
		Q[23, 4] <- pl;  Q[23, 3] <- p0;  Q[23,27] <- pr
		Q[24,29] <- pl;  Q[24,28] <- p0
		Q[25,23] <- pl;  Q[25,22] <- p0
		Q[26,16] <- pl;  Q[26,15] <- p0
		Q[27,17] <- pl;  Q[27,15] <- p0
		Q[28, 2] <- pl;  Q[28, 1] <- p0
		Q[29, 4] <- pl;  Q[29, 3] <- p0
	}
	
	# 8 on the same side
	if (type=="14") {
		dimQ <- 15
		pl <- pnorm(   0, mean=mu ) - pnorm( -3*c, mean=mu)
		pr <- pnorm( 3*c, mean=mu ) - pnorm(    0, mean=mu)
		
		Q <- diag(0,dimQ)
		
		Q[ 1, 2] <- pl;  Q[ 1, 9] <- pr
		Q[ 2, 3] <- pl;  Q[ 2, 9] <- pr
		Q[ 3, 4] <- pl;  Q[ 3, 9] <- pr
		Q[ 4, 5] <- pl;  Q[ 4, 9] <- pr
		Q[ 5, 6] <- pl;  Q[ 5, 9] <- pr
		Q[ 6, 7] <- pl;  Q[ 6, 9] <- pr
		Q[ 7, 8] <- pl;  Q[ 7, 9] <- pr
		                 Q[ 8, 9] <- pr
		Q[ 9, 2] <- pl;  Q[ 9,10] <- pr
		Q[10, 2] <- pl;  Q[10,11] <- pr
		Q[11, 2] <- pl;  Q[11,12] <- pr
		Q[12, 2] <- pl;  Q[12,13] <- pr
		Q[13, 2] <- pl;  Q[13,14] <- pr
		Q[14, 2] <- pl;  Q[14,15] <- pr
		Q[15, 2] <- pl;
	}

	# ... on the same side (general approach)
	if ( regexpr("SameSide", type)>0 ) {
		anzahl <- as.numeric(gsub("SameSide", "", type))
		dimQ <- 2*anzahl - 1
		hdQ <- anzahl - 1
		pl <- pnorm(   0, mean=mu ) - pnorm( -3*c, mean=mu)
		pr <- pnorm( 3*c, mean=mu ) - pnorm(    0, mean=mu)		
		Q <- diag(0, dimQ)
                for ( i in 1:hdQ ) {
                  Q[i,i+1]     <- pl
                  Q[hdQ+i+1,2] <- pl
                  Q[i,hdQ+2] <- pr
                  Q[hdQ+i,hdQ+i+1] <- pr
                }
	}
	
	# 2 of 2 beyond +-2 sigma
	if (type=="15") {
		dimQ <- 3
		pl <- pnorm( -2*c, mean=mu ) - pnorm( -3*c, mean=mu)
		p0 <- pnorm(  2*c, mean=mu ) - pnorm( -2*c, mean=mu)
		pr <- pnorm(  3*c, mean=mu ) - pnorm(  2*c, mean=mu)
		
#                     1   2   3
#		      00  10  01
#		1 00  p0  pr  pl
#		2 10  p0   0  pl
#               3 01  p0  pr   0

		Q <- diag(0,dimQ)
		
		Q[1,2] <- pl;  Q[1,1] <- p0;  Q[1,3] <- pr
		               Q[2,1] <- p0;  Q[2,3] <- pr
	        Q[3,2] <- pl;  Q[3,1] <- p0;
	}
	
	# 3 of 3 beyond +-3 sigma
	if (type=="19") {
		dimQ <- 5
		pl <- pnorm( -3*c, mean=mu )
		p0 <- pnorm(  3*c, mean=mu ) - pnorm( -3*c, mean=mu)
		pr <- 1 - pnorm(  3*c, mean=mu)
		
#                      1    2    3    4    5
#		       0000 1000 1100 0010 0011
#		1 0000 p0   pr   0    pl   0
#		2 1000 p0   0    pr   pl   0
#               3 1100 p0   0    0    pl   0
#		4 0010 p0   pr   0    0    pl
#		5 0011 p0   pr   0    0    0

		Q <- diag(0,dimQ)
		
		Q[1,4] <- pl;  Q[1,1] <- p0;  Q[1,2] <- pr
		Q[2,4] <- pl;  Q[2,1] <- p0;  Q[2,3] <- pr
	        Q[3,4] <- pl;  Q[3,1] <- p0;
		Q[4,5] <- pl;  Q[4,1] <- p0;  Q[4,2] <- pr
		               Q[5,1] <- p0;  Q[5,2] <- pr
	}
	
	Q
}