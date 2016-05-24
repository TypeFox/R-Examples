# 
# Maarten Speekenbrink
# 
# FORWARD-BACKWARD algoritme, 23-3-2008
# 

fb <- function(init,A,B,ntimes=NULL,return.all=FALSE,homogeneous=TRUE,useC=TRUE,na.allow=TRUE) {

	# Forward-Backward algorithm (used in Baum-Welch)
	# Returns alpha, beta, and full data likelihood
	
	# NOTE THE CHANGE IN FROM ROW TO COLUMN SUCH THAT TRANSPOSING A IS NOT NECCESSARY ANYMORE
	# IN COMPUTING ALPHA AND BETA BUT IS NOW NECCESSARY IN COMPUTING XI
	# A = T*K*K array with transition probabilities, from row to column!!!!!!!
	# B = T*D*K matrix with elements ab_{ij} = P(y_i|s_j)
	# init = N*K vector with initial probabilities

    # T = total number of time points
    # K = number of states
    # D = dimension of observations (D>1 is multivariate)
    # N = number of participants
     
	# NOTE: to prevent underflow, alpha and beta are scaled, using sca
	
	# NOTE: xi[t,i,j] = P(S[t] = j & S[t+1] = i) !!!NOTE the order of i and j!!!
	
	#nt <- nrow(B)
	nt <- dim(B)[1]
	ns <- ncol(init)
	
	if(na.allow) B <- replace(B,is.na(B) & !is.nan(B),1)

	B <- apply(B,c(1,3),prod)
	
	if(is.null(ntimes)) ntimes <- nt
	
	lt <- length(ntimes)
	et <- cumsum(ntimes)
	bt <- c(1,et[-lt]+1)
	
 	if(useC) {
		
		alpha <- matrix(0,ncol=ns,nrow=nt)
		sca <- rep(0,nt)
		
		beta <- matrix(0,ncol=ns,nrow=nt)
		xi <- array(0,dim=c(nt,ns,ns))
		
		res <- .C("forwardbackward",
			hom=as.integer(homogeneous),
			ns=as.integer(ns),
			lt=as.integer(lt),
 			nt=as.integer(nt),
 			ntimes=as.integer(ntimes),
			bt=as.integer(bt),
			et=as.integer(et),
 			init=as.double(t(init)),
 			A=as.double(A),
 			B=as.double(t(B)),
 			alpha=as.double(alpha),
 			beta=as.double(beta),
			sca=as.double(sca),
			xi=as.double(xi),
 			PACKAGE="depmixS4")[c("alpha","beta","sca","xi")]
		
		alpha <- matrix(res$alpha,ncol=ns,byrow=TRUE)
		beta <- matrix(res$beta,ncol=ns,byrow=TRUE)
		xi <- array(res$xi,dim=c(nt,ns,ns))
		xi[et,,] <- NA
		sca <- res$sca
				
	} else {
		
		alpha <- matrix(ncol=ns,nrow=nt)
		beta <- matrix(ncol=ns,nrow=nt)
		sca <- vector(length=nt)
		xi <- array(dim=c(nt,ns,ns))
		
		for(case in 1:lt) {
			alpha[bt[case],] <- init[case,]*B[bt[case],] # initialize
			sca[bt[case]] <- 1/sum(alpha[bt[case],])
			alpha[bt[case],] <- alpha[bt[case],]*sca[bt[case]]
						
			if(ntimes[case]>1) {
				for(i in bt[case]:(et[case]-1)) {
					if(homogeneous) alpha[i+1,] <- (A[1,,]%*%alpha[i,])*B[i+1,]
					else alpha[i+1,] <- (A[i,,]%*%alpha[i,])*B[i+1,]
					sca[i+1] <- 1/sum(alpha[i+1,])
					alpha[i+1,] <- sca[i+1]*alpha[i+1,]
				}
			}
			
			beta[et[case],] <- 1*sca[et[case]] # initialize
						
			if(ntimes[case]>1) {
				for(i in (et[case]-1):bt[case]) {
					if(homogeneous) beta[i,] <-(B[i+1,]*beta[i+1,])%*%A[1,,]*sca[i]
					else beta[i,] <-(B[i+1,]*beta[i+1,])%*%A[i,,]*sca[i]
				}
				
				for(i in bt[case]:(et[case]-1)) {
					if(homogeneous) xi[i,,] <- rep(alpha[i,],each=ns)*(B[i+1,]*beta[i+1,]*A[1,,])
					else xi[i,,] <- rep(alpha[i,],each=ns)*(B[i+1,]*beta[i+1,]*A[i,,])
				}
			}
			
		}
	}
	
	gamma <- alpha*beta/sca
	like <- -sum(log(sca))
	
	if(return.all) {
		res <- list(alpha=alpha,beta=beta,gamma=gamma,xi=xi,sca=sca,logLike=like)
	} else {
		res <- list(gamma=gamma,xi=xi,logLike=like)
	}
	
	res
}

