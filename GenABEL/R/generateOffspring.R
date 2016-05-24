#' simulates offspring's genotypes 
#' 
#' Simulat offspring's genotypes given genotypes of 
#' the parents. No LD is assumed.
#' 
#' @param g1 vector of genotypes of parent 1, coded 
#' with 0, 1, 2 and NA
#' @param g2 vector of genotypes of parent 2, coded 
#' with 0, 1, 2 and NA
#' @param q vector of the frequencies of effects 
#' alleles (used to simulate missings in parental 
#' genotypes). If NULL, c(0.25,0.5,0.25) is used.
#' 
#' @return a vector containing simulated genotypes 
#' of offspring
#' 
#' @author Yurii Aulchenko
#' 
#' @examples 
#' eaf <- runif(10)
#' g1 <- rbinom(10,2,eaf)
#' g2 <- rbinom(10,2,eaf)
#' generateOffspring(g1,g2)
#' 
generateOffspring <- function(g1,g2,q=NULL) {
	tm <- array(dim=c(3,3,3))
	tm[1,1,] <- c(1,1,1) #c(1,0,0)
	tm[1,2,] <- tm[2,1,] <- c(.5,1,1) #c(.5,.5,0)
	tm[1,3,] <- tm[3,1,] <- c(0,1,1) #c(0,1,0)
	tm[2,2,] <- c(.25,.75,1) #c(.25,.5,.25)
	tm[2,3,] <- tm[3,2,] <- c(0,.5,1) #c(0,.5,.5)
	tm[3,3,] <- c(0,0,1) #c(0,0,1)
	genSOG <- function(localg) {
		vc <- tm[localg[1],localg[2],]
		#print(vc)
		rn <- runif(1)
		#print(rn)
		gt <- which(vc>=rn)[1]-1
		#print(gt)
		gt
	}
# infer missing genotypes, if any
	if (any(is.na(c(g1,g2)))) {
		warning("NA's in g1 and/or g2; inferring")
		infer <- function(q) {
			return(g <- rbinom(length(q),2,q))
		}
		if (!is.null(q)) {
			miss <- which(is.na(g1))
			g1[miss] <- rbinom(length(miss),2,q[miss])
			miss <- which(is.na(g2))
			g2[miss] <- rbinom(length(miss),2,q[miss])
		} else {
			miss <- which(is.na(g1))
			g1[miss] <- rbinom(length(miss),2,rep(0.5,length(miss)))
			miss <- which(is.na(g2))
			g1[miss] <- rbinom(length(miss),2,rep(0.5,length(miss)))
		}
	}
# generate offspring genotypes
	g1g2 <- t(matrix(c(g1+1,g2+1),ncol=2))
	res <- apply(g1g2,MARGIN=2,FUN=genSOG)
	res
}
