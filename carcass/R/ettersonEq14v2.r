ettersonEq14v2 <- function(s,f,J){
#This function calculates the function gfe in Eq. 14 of Etterson, M. 
#Hidden Markov models for estimating animal mortality from anthropogenic 
#hazards.  It has been modified to allow time-varying parameters.
#Parameters are assumed to be indexed to the time since carcass death.
#
#Inputs:
#   s = vector of estimated daily persistence probabilities
#   f = vector of discovery probabilities
#   J = row-vector of intervals between searches
#
#Outputs:
#   gfe = function related to the probability of sampling a carcass 
#   killed between the first and last searches. See ms for details.


#	Note: this function requires expm

n <- length(J)#number of searches
N <- sum(J)#number of days of monitoring

Searches<-cumsum(J)
gfe <- 0
V1<- matrix( c(1,0,0), ncol=3, byrow=TRUE)
V3<- matrix( c(0,0,1), ncol=3, byrow=TRUE)
#Indexing below is for the interval in which an animal is killed.

for (j in 1:N){ #Day the carcass enters the area searched
    A <- diag(3)
	srch_indx <- 0
    for (i in j:N){ #daily fates of the carcass after entering search area
		age <- i-j+1
        qr_i <- s[age]
        pr_i <- 1-qr_i
		Svec <- c(qr_i,pr_i,0,0,1,0,0,0,1)
		S <- matrix( Svec, ncol=3, byrow=TRUE )
        A <- A %*% S
        x <- which(i==Searches)
        if (length(x)!=0){
			srch_indx <- srch_indx + 1
            pd_i <- f[srch_indx]
            qd_i = 1-pd_i
			Dvec <- c(qd_i,0,pd_i,0,1,0,0,0,1)
			D <- matrix( Dvec, ncol=3, byrow=TRUE )
            A <- A %*% D
        }
    }
    gfe <- gfe + V1 %*% A %*% t(V3)
}
p <- gfe/sum(J)
return(p)
}