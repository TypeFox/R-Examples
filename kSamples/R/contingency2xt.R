# Project:  2 * t Contingency Table
# adapted from Angie Zhu's getTable2xtNull.R
# Filename: Contingency2xt.R
# Last modified: 03.23.2012

# This function computes or simulates the null distribution of
# the Kruskal-Wallis statistics in a 2 x t contingency table.  

# Treatment |   1   2   ...    t  | Total
# ---------------------------------------
# Response  |                     |
#    a      | A_1  A_2   ...  A_t | m
#    b      | B_1  B_2   ...  B_t | n
# ---------------------------------------
# Total     | d_1  d_2   ...  d_t | N

# Arguments:  
#		Avec: integer vector of length tnum, containing the 
#                "a" reponses
#		Bvec: integer vector of length tnum, containing the 
#                "b" reponses
#		method: "asymptotic", "simulated", "exact", 
#                  indicating the method of p-value calculation,
#		dist: FALSE (default) means that the exact or
#                estimated null distribution is not requested, 
#                otherwise it is. 
#          tab0: TRUE (default), when method = "simulated",
#			 the null distribution is given in tabular form,
#                when tab0 = FALSE and method = "simulated", the
#                null distribution is given as a single vector 
#                of all simulated values. 
#		Nsim: number of simulations to run, default 10^6
#
#
# Output: a list of class "kSamples" with components
#        test.name  = "2 x t Contingency Table"
#                t  =  number of columns in the table
#          KW.cont  = 2 (or 3 vector) containing the observed KW
#                     value, its asymptotic (and the simulated 
#                     or exact) P-value
#          null.dist  = a 2 x M matrix giving the M unique 
#                       ordered KW values of the null 
#                       distribution in the first column, and 
#                       the corresponding probabilities 
#                       (simulated or exact) in the second 
#                       column. When method = "simulated" and 
# 			        tab0 = FALSE null.dist consists of a
#                       single vector of all simulated values.
#			When dist = NULL or method = "asymptotic", only 
#                           NULL is returned.
#        method = the method used
#        Nsim = number of simulations employed


contingency2xt <- function(Avec, Bvec, method=c("asymptotic","simulated","exact"), 
           		dist = FALSE, tab0 = TRUE, Nsim=1e6) {
    	method <- match.arg(method)
	m <- sum(Avec)
	n <- sum(Bvec)
    	tnum <- length(Avec)
	N <- m + n
	dvec <- Avec + Bvec
	KW.obs <- (N * (N - 1) / (m * n)) * (sum(Avec^2 / dvec) - 
                  	m^2/N )
	pValue.asy <- 1 - pchisq(KW.obs,tnum-1)
	if(dist == TRUE) Nsim <- min(Nsim,1e8) 
	# limits the size of null.dist
	# whether method = "exact" or = "simulated"
	if(method == "asymptotic"){
		KW.cont <- c(KW.obs,pValue.asy)
		null.dist <- NULL
	}else{
		ncomb <- choose(m + tnum - 1, tnum - 1)
		if (ncomb <= Nsim && method == "exact") { 
			if(dist){ 
					ans <- numeric(2+2*ncomb)
				}else{
					ans <- numeric(2)
			}
			out <- .C("contingency2xtExact0",
					Avec = as.integer(Avec), Bvec = as.integer(Bvec),
					tnum = as.integer(tnum), 
					ncomb = as.integer(ncomb),
					getDist = as.integer(dist),
                                	ans = as.double(ans))
				KW.obs <- out$ans[1]
				KW.obs <- (N * (N - 1) / (m * n)) * (KW.obs - 
                                 	m^2 / N) 
				KW.obs <- round(KW.obs, 6)
       				pValue <- out$ans[2]		
				if(dist){	
					out <- out$ans[-(1:2)]
					out <- matrix(out, nrow=ncomb, ncol=2, byrow=FALSE)
					out <- out[out[, 2] > 0, ]
					prob <- out[ , 2]
					KW <- (N * (N - 1) / (m * n)) * (
                                    		out[ , 1] - m^2 / N) 
					KW <- round(KW, 6)
					KW.u <- sort(unique(KW))
					k.u <- length(KW.u)
					prob.u <- sapply(1:k.u, function(j) {
		 					sum( prob[KW == KW.u[j]] ) })
					null.dist <- cbind(KW.u,prob.u)
                    			dimnames(null.dist) <- list(NULL,c("KW","prob"))
				}else{
					null.dist <- NULL
				}
				KW.cont <- c(KW.obs,pValue.asy,pValue)
		}else{
			method <- "simulated"
		}
		if( method == "simulated" ){
			if( Nsim < 100) Nsim <- 100
			if(dist){ 
					ans <- numeric(2+Nsim)
				}else{
					ans <- numeric(2)
			}
			out <- .C("contingency2xtSim0",
					Avec = as.integer(Avec), Bvec = as.integer(Bvec),
					tnum = as.integer(tnum), 
					Nsim = as.integer(Nsim),
					getDist = as.integer(dist),
                                	ans = as.double(ans))	
			KW.obs <- out$ans[1]
			KW.obs <- (N * (N - 1) / (m * n)) * (KW.obs - m^2 / N) 
			KW.obs <- round(KW.obs, 6)
      			pValue <- out$ans[2]		
			if(dist){
				KW <- (N * (N - 1) / (m * n)) * 
                           		(round(out$ans[3:(Nsim+2)],6) - m^2 / N)
				if(tab0==TRUE){
					tab <- table(KW)/Nsim
					null.dist <-
                                   cbind(as.numeric(names(tab)),
                                            as.vector(tab))
                    		dimnames(null.dist) <- 
                             list(NULL,c("KW","prob"))}else{
					null.dist <- KW
				}
			}else{
				null.dist <- NULL
			}
			KW.cont <- c(KW.obs,pValue.asy,pValue)
		}
	}
    if(method=="asymptotic"){
		names(KW.cont) <- c("observed KW",
                               "  asympt. P-value")
    }
    if(method=="exact"){
		names(KW.cont) <- c("observed KW",
                           "  asympt. P-value","exact P-Value")
    }
    if(method=="simulated"){
		names(KW.cont) <- c("observed KW",
                           "  asympt. P-value","sim. P-Value")
    }
    object <- list(test.name =paste("2 x t Contingency Table"),
                     t = tnum, KW.cont = KW.cont, 
                     null.dist = null.dist, 
				method = method, Nsim = Nsim)
    class(object) <- "kSamples"
    object
}


