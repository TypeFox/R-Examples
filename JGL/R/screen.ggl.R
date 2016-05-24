screen.ggl <-
function(Y,lambda1,lambda2,weights="equal")
{
	p = dim(Y[[1]])[2]
	K = length(Y)
	n = rep(0,K)
	for(k in 1:K) {n[k] = dim(Y[[k]])[1]}

	# mean-normalize Y:
	for(k in 1:K){
	for(j in 1:p){
		Y[[k]][,j] = Y[[k]][,j]-mean(Y[[k]][,j])
	}}

	# set weights:
	if(length(weights)==1){if(weights == "equal"){
		weights = rep(1,K)
	}}
	if(length(weights)==1){if(weights == "sample.size"){
		weights = n/sum(n)
	}}

	Sabs=list()
	connected = rep(TRUE,p)

	# get p*K matrix of maximum absolute values of off-diagonal elements columns of cov(Y):
	above.threshold = matrix(0,nrow=p,ncol=K)
	for(j in 1:p){
	if((j%%100==0)&(j>0)) {print(paste("screening feature",j))}
	for(k in 1:K){
		# get absolute value of column i of S:
		Sabs[[k]] = 1/(n[k])*Y[[k]][,j]%*%Y[[k]][,setdiff(1:p,j)]
	}
	#now check that Sabs meets the criteria:
	if(length(lambda1)==1) { L1 = 1/weights[k]*lambda1 }
	if(length(lambda1)>1) { L1 = 1/weights[k]*lambda1[j,setdiff(1:p,j)] }
	if(length(lambda2)==1) { L2 = 1/weights[k]*lambda2*(K-1) }
	if(length(lambda2)>1) { L2 = 1/weights[k]*lambda2[j,setdiff(1:p,j)]*(K-1) }

	# criterion vector: evaluates for all of Sabs vector whether each element meets the criteria: 1: |S|<lam1+lam2.  
	criterion = Sabs[[1]]*0
	for(k in 1:K) { criterion = criterion + ( abs(Sabs[[k]]) > (L1)) }  # necessary, but not sufficient for connected
	# is criterion met?
	if(sum(criterion)==0) {connected[j] = FALSE}

	}
#	connected = rowSums(above.threshold)>0
	return(connected=connected)
}

