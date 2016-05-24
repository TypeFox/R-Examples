weight.seq <-
function(T,C,mu,sigma,n){
	N=length(mu)
	m=length(T)

	# A matrix that contains the values f_j(x_i) in the algorithm.
	A=matrix(0,m,N)
		for(i in 1:m){
		for(j in 1:N){
			A[i,j]=C[i]*dnorm(T[i],mu[j],sigma)
	}}

	# A starting value for the weigth vector
	theta.old=rep(1/N,N)
	abs.diff=1
	epsilon=0.001
	count=0		
	while((abs.diff>epsilon)&(count<500)){

		#Construct the vector of w's
		B=t(t(A)*theta.old)
		
		B=B/(apply(B,1,sum))   ##1 indicates rows, 2 indicates columns

		#The sums over i=1,..,n of the w_ij's
		w=apply(B,2,sum)

		#The system of equations for finding the lambda_hat values to maximize the Q function

		flambda=function(lambda){
			f1<-sum(w/(lambda[1]+mu*lambda[2]+mu^2*(m-lambda[1])/(1-sigma^2)))
			f2<-sum(w*mu/(lambda[1]+mu*lambda[2]+mu^2*(m-lambda[1])/(1-sigma^2)))
			f<-(f1-1)^2+f2^2
			return(f)
		}

		# If function f is equal to zero for some values of lambda, then these lambda solve the system of equations for maximization
		# So find the values of lambda that minimize the non-negative function f
		lambda.hat=nlm(flambda,c(m,1),iterlim=50)$estimate

		#Find the value of theta for the next iteration
		theta.new=w/(lambda.hat[1]+mu*lambda.hat[2]+mu^2*(m-lambda.hat[1])/(1-sigma^2))

		#The max difference of the two consecutive values of theta
		abs.diff=max(abs(theta.new-theta.old))
		if(is.na(abs.diff)){
		abs.diff=0
		theta.new=theta.old}
		
		#Reset the value of theta as the old value to go to the next iteration 
		theta.old=pmax(theta.new,0)
		theta.old=pmin(theta.old,1)
		count=count+1
	}

	############## The EM algorithm stops here giving the maximum likelihood estimate of components theta in theta.old
	#Check whether the value of N is sufficient for the estimation

	theta.hat=pmax(theta.new,0)
	theta.hat=pmin(theta.hat,1)

	return(theta.hat)
}
