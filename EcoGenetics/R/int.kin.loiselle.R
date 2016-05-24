#' obtetion of multilocus Loiselle's Fij matrix 
#' @param eco Object of class ecogen. 
#' @keywords internal

# depends of geno and locus
# geno is: eco@A
# locus is: as.numeric(eco@INT@loc.fac)

int.kin.loiselle <- function(geno, locus) {
	
	nloc <- max(locus)
	nal <- ncol(geno)
	N <- nrow(geno)
	alelos <- nrow(geno)
	p.la <- apply(geno, 2, sum, na.rm = TRUE)
	p.locus <- tapply(p.la,locus, sum, na.rm = TRUE)
	p.la <- p.la / p.locus[locus]
	
	p.la.cen <- t(replicate(N, p.la))
	p.center <- geno - p.la.cen
	
	n1 <- 2 * p.locus
	n1 <- n1[locus]

	

	####SUPERIOR COEFFICENT###
	
	#RIGHT TERM	
	T2  <- p.la * (1 - p.la) / (n1 - 2)
	
	#LEFT TERM
	out <- list()
	for(i in 1:nal) {
		out[[i]] <- outer(p.center[,i], p.center[,i])
	}
	
	numer <- matrix(0, N , N)
	for(i in 1:nal){
		numer <- numer + ifelse(is.na(out[[i]]), 0, out[[i]] + matrix(T2[i], N, N))
	}
	
	###INFERIOR COEFFICENT##
	denom<- sum(p.la * (1 - p.la))
	
	#Fij COMPUTATION
	Fij <- (numer / denom)
	
	Fij
  	
}
