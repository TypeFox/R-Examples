numericalDNA <- function(background){

		m <- c(0,0,1,
		   -sqrt(2)/3, sqrt(6)/3, -1/3,
		   -sqrt(2)/3, -sqrt(6)/3, -1/3,
			2*sqrt(2)/3, 0, -1/3,
	(c(0,0,1))*background[1]+
	(c(-sqrt(2)/3, sqrt(6)/3, -1/3))*background[2]+
	(c(-sqrt(2)/3, -sqrt(6)/3, -1/3))*background[3]+
	(c(2*sqrt(2)/3, 0, -1/3))*background[4], 0,0,1,-sqrt(2)/3, sqrt(6)/3, -1/3,
		   -sqrt(2)/3, -sqrt(6)/3, -1/3, 2*sqrt(2)/3, 0, -1/3)
	
	outm <- t(matrix(data=m, nrow=3, ncol=9))
	rownames(outm)<- c("A","C","G","T","-", "a", "c", "g", "t")
	outm	
}



