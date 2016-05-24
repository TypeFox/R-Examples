## Update on Oct 07, 2015

#' Implement physician debias algorithm
#' 
#' This function implements physician debias algorithm proposed in Salter-Townshend 
#' and Murphy (2013).
#'
#' @param data The original data to be used. It is suggested to use similar
#' input as InterVA4, with the first column being death IDs. The only
#' difference in input is InsilicoVA takes three levels: ``present'',
#' ``absent'', and ``missing (no data)''. Similar to InterVA software,
#' ``present'' symptoms takes value ``Y''; ``absent'' symptoms take take value
#' ``NA'' or ``''. For missing symptoms, e.g., questions not asked or answered
#' in the original interview, corrupted data, etc., the input should be coded
#' by ``.'' to distinguish from ``absent'' category. The order of the columns does
#' not matter as long as the column names are correct. Currently it cannot other
#' non-symptom columns such as subpopulation. And the first column should be 
#' the death ID. Everything other than the death ID, physician ID, and physician 
#' codes should be symptoms.
#' @param phy.id vector of column names for physician ID
#' @param phy.code vector of column names for physician code
#' @param phylist vector of physician ID used in physician ID columns
#' @param causelist vector of causes used in physician code columns
#' @param tol tolerance of the EM algorithm
#' @param max.itr maximum iteration to run
#' @param verbose logical indicator for printing out likelihood change

#' @return \item{code.debias}{Individual cause likelihood distribution}
#' @return \item{csmf}{Cause specific distribution in the sample}
#' @return \item{phy.bias}{Bias matrix for each physician}
#' @return \item{cond.prob}{Conditional probability of symptoms given causes}
#' 
#' @examples 
#' 
#' data(RandomPhysician)
#' head(RandomPhysician[, 1:10])
#' \dontrun{
#' causelist <- c("Communicable", "TB/AIDS", "Maternal", 
#'                "NCD", "External", "Unknown")
#' phydebias <- physician_debias(RandomPhysician, phy.id = c("rev1", "rev2"), 
#' phy.code = c("code1", "code2"), phylist = c("doc1", "doc2"), 
#' causelist = causelist, tol = 0.0001, max.itr = 5000)
#' 
#' # see the first physician's bias matrix
#' round(phydebias$phy.bias[[1]], 2)
#' }

#' @references M. Salter-Townshend and T. B. Murphy (2013).\emph{Sentiment 
#' analysis of online media}. \cr \emph{In Algorithms from and for Nature and 
#' Life, pages 137-145, Springer.}
#' 
physician_debias <- function(data, phy.id, phy.code, phylist, causelist, tol = 0.0001, max.itr = 5000, verbose = FALSE){
	# which columns are doctor ID
	docCol <- match(phy.id, colnames(data))
	causeCol <- match(phy.code, colnames(data))

	# delete the data without any physician codes
	noDoc <- apply(data[, docCol], 1, function(x){
		length(which(x %in% phylist == FALSE)) == length(phy.id)
		})
	if(sum(noDoc) == dim(data)[1]){
		stop("No physician coding found for the list of physician IDs")
	}
	if(sum(noDoc) > 0){
		data <- data[-which(noDoc == TRUE), ]
	}
	id <- data[, 1]
	# number of death
	A <- dim(data)[1]
	# number of symptoms
	N <- dim(data)[2] - length(phy.id) - length(phy.code) - 1
	# refine doctor list
	doclist <- intersect(phylist, unique(unlist(data[, docCol])))
	# number of doctors
	K <- length(doclist)
	
	# get numeric doctor matrix
	doc.indiv  <- matrix(0, dim(data)[1], length(docCol))
	for(i in 1:length(docCol)){
		doc.indiv[,i] <- match(data[, docCol[i]], doclist)	
	}
	# get numeric cause matrix
	cause.indiv  <- matrix(0, dim(data)[1], length(causeCol))
	for(i in 1:length(docCol)){
		cause.indiv[,i] <- match(data[, causeCol[i]], causelist)	
	}

	# number of causes
	J <- length(causelist)
	
	#  list of #Doc, each element is a matrix of A * J (#Death * #Cause)
	Y <- vector("list", K) 
	for(i in 1 : K){
		Y[[i]] <- matrix(0, A, J)
		for(j in 1:dim(doc.indiv)[2]){
			which.death <- which(doc.indiv[,j] == i)
			which.diag <- cause.indiv[which.death, j]
			for(k in 1:length(which.death)){
				Y[[i]][which.death[k], which.diag[k]] <- 1
			}
		}
	}

	# Symptom Matrix of A * N (#Death * #Symptom)
	W <-matrix(0, A, N)
	W.raw <- data[, -c(1, docCol, causeCol)]
	W[which(W.raw == "Y")] <- 1
	W[which(W.raw == ".")] <- -1
	noSymp <- which(apply(W, 2, sum) + A == 0)
	N <- N - length(noSymp)
	W <- W[, -noSymp]

	###################################################################
	## 						Initialization                           ## 
	###################################################################
	T.hat<-matrix(0, A, J)          
	for(k in 1 : K){
	  T.hat <- T.hat + Y[[k]]
	}
	K.indiv <- apply(T.hat,1,sum)
	for(i in 1 : A){
		if(K.indiv[i] > 0){
		   T.hat[i,] <- T.hat[i,] / K.indiv[i]
		}		
	}

	T.hat.ori <- T.hat
	p.hat <- apply(T.hat, 2, mean)
	p.hat.original <- p.hat
	
	score1 <- 0
	diff <- 1
	index <- 1
	p.store <- NULL
	Pi.store <- NULL
	theta.store <- NULL
	T.store <- NULL

	###################################################################
	##                     EM Algorithm                              ##
	###################################################################
	loglik<-function(p.hat, Pi.hat, theta.hat, T.hat, A, K, J, N, W, Y){
	  sumall<-0
	  for(i in 1:J){
	    sum1<-rep(log(p.hat[i]),A)
	    sum2<-rep(0,A)
	    for(k in 1:K){
	      for(j in 1:J){
	        if(Pi.hat[[k]][i,j] != 0 && Pi.hat[[k]][i,j] != 1){
	        	sum2<-sum2+Y[[k]][,j]*log(Pi.hat[[k]][i,j]) 
	        }	             
	        if("NaN" %in% sum2)cat(c("i=",i,"j=", j,"k=", k,"\n"))
	      }
	    }
	    sum3<-rep(0,A)
	    for(n in 1:N){
	      if(theta.hat[n,i] != 0 && theta.hat[n, i] != 1 ){
	      	sum3[which(W[,n] >=0)] <- sum3[which(W[,n] >=0)] + 
	      							W[which(W[,n] >=0), n]*log(theta.hat[n, i]) +
	      		  		            (1 - W[which(W[,n]>=0), n])*log(1-theta.hat[n,i])
	      }
	      if("NaN" %in% sum3) cat(c("i=",i,"n=",n ,"\n"))
	    }
	    sum1<-sum1+sum2+sum3 
	    
	    sumall<-sumall+sum(sum1*T.hat[,i])
	    
	  }
	  return(sumall)
	} 


	while(diff >= tol && index < max.itr) {
	  
	  Pi.hat <- list(K)         

	  for(k in 1:K){          
	    Pi.hat[[k]] <- matrix(0, J, J)
	    for(i in 1:J){
	      for(j in 1:J){
	        Pi.hat[[k]][i, j]= T.hat[,i] %*% Y[[k]][,j]            
	      }
	      if(sum(Pi.hat[[k]][i, ]) == 0){
	      	normal.const <- 1
	      }else {
	      	normal.const <- sum(Pi.hat[[k]][i,])
	      }	
	      Pi.hat[[k]][i,] <- Pi.hat[[k]][i,] / normal.const
	    }
	  }
	  
	  theta.hat<-matrix(0, N, J)         # Estimate theta.hat    
	  for(i in 1:N){
	    for(j in 1:J){
	      theta.hat[i, j] <- sum(W[which(W[, i] >= 0),i] * T.hat[which(W[,i] >=0), j])
	      if(sum(T.hat[which(W[, i] >=0), j]) > 0){
	      	theta.hat[i, j] <- theta.hat[i, j] / sum(T.hat[which(W[, i] >=0),j]) 
	      }
	    }
	  }
	   
	  p.hat<-array(0, J)                # Estimate NEW p.hat
	  for(j in 1:J){
	    p.hat[j]<-sum(T.hat[, j])/A
	  }
	 
	  T.hat<-matrix(0, A, J)             # Estimate NEW T.hat   
	  
	  for(i in 1 : J){
	    sum1<-rep(1, A)
	    for(n in 1 : N){
	      notmiss <- which(W[,n] >= 0)
	      sum1[notmiss] <- sum1[notmiss] * theta.hat[n,i]^W[notmiss, n] * (1-theta.hat[n,i])^(1-W[notmiss, n])
	    }
	    
	    sum2<-1
	    for(k in 1:K){
	      for(j in 1:J){
	        sum2<-sum2 * Pi.hat[[k]][i, j] ^ Y[[k]][ ,j]     
	      }
	    }
	    
	    T.hat[,i] <- sum2 * sum1 * p.hat[i]
	  }

	  for(i in 1:A){ 
	  	if(sum(T.hat[i, ]) > 0) T.hat[i,]=T.hat[i,]/sum(T.hat[i,])
	  }
	  
	  p.store[[index]]<-p.hat
	  Pi.store[[index]]<-Pi.hat
	  theta.store[[index]]<-theta.hat
	  T.store[[index]]<-T.hat

	  score2 <- loglik(p.hat,Pi.hat,theta.hat,T.hat, A, K, J, N, W, Y)
	  diff.lik <- score2 - score1
	  score1 <- score2


	  if(index > 1){
		  diff<-max(abs(unlist(p.store[[index-1]])-unlist(p.store[[index]])),
	            abs(unlist(Pi.store[[index-1]])-unlist(Pi.store[[index]])),
	            abs(unlist(theta.store[[index-1]])-unlist(theta.store[[index]])),
	            abs(unlist(T.store[[index-1]])-unlist(T.store[[index]])))	
	  }else{
	  	diff <- diff.lik
	  }
	  
	  if(verbose){
		  cat(paste("itr",index, 
	  			"loglik diff:",round(diff.lik, 4),
	  			"Parameter diff:", round(diff, 4), "\n"))
	  }else{
	  	cat(".")
	  }
	  
	  diff <- max(abs(diff), abs(diff.lik))    
	  index<-index+1
	}
	if(index>=5000 && diff>0.0001) cat("Diverge")
	
	final<-list(4)
	colnames(T.hat) <- colnames(theta.hat) <- names(p.hat) <- causelist
	names(Pi.hat) <- phylist
	
	for(i in 1:length(Pi.hat)){
		colnames(Pi.hat[[i]]) <- rownames(Pi.hat[[i]]) <- causelist
	}

	T.hat <- data.frame(T.hat, check.names=FALSE)
	T.hat <- cbind(ID = id, T.hat)
	final <- list(code.debias = T.hat, 
				  csmf = p.hat, 
				  phy.bias = Pi.hat, 
				  cond.prob  = theta.hat)
	return(final)
}