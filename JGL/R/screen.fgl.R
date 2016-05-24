# status 7/9: not really updated, just a couple notes on possible bugs and places where it should be fixed.
screen.fgl <-
function(Y,lambda1,lambda2,weights="equal")
{
	K = length(Y)
	p = dim(Y[[1]])[2]
	n = rep(0,K)
	for(k in 1:K) {n[k] = dim(Y[[k]])[1]}

	# mean-normalize Y:
	for(k in 1:K){
	#Y[[k]] = scale(Y[[k]],center=TRUE,scale=FALSE)
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

	if(K==2)
	{
	# criterion vector: evaluates for all of Sabs vector whether each element meets the criteria: 1: |S|<lam1+lam2.  2: |S1 + ... SK| < K*lam1.
	criterion = Sabs[[1]]*0
	for(k in 1:K) { criterion = criterion + ( abs(Sabs[[k]]) > (L1+L2)) }
	Sabs.all = Sabs[[1]]
	for(k in 2:K) {Sabs.all = Sabs.all + Sabs[[k]]}
	criterion = criterion + (abs(Sabs.all) > K*L1)
	# is criterion met?
	if(sum(criterion)==0) {connected[j] = FALSE}
	}
	if(K>2)
	{
	# criterion vector: evaluates for all of Sabs vector whether each element meets the sufficient criterion: 1: |S|<lam1
	criterion = Sabs[[1]]*0
	for(k in 1:K) { criterion = criterion + ( abs(Sabs[[k]]) > (L1)) }
	# is criterion met?
	if(sum(criterion)==0) {connected[j] = FALSE}
	}
	}
	return(connected=connected)
}

