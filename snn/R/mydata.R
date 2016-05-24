mydata <-
function(n, d, mu = 0.8, portion = 1/2){
	#data generation function: f1~N(0,I_d); f2~N(mu,I_d) with mu = c(0.8, ..., 0.8) and prior probability pi=1/2.
	
	n1 = floor(n*portion)
	n2 = n - n1	
	X1 = matrix(rnorm(n1*d),n1,d)      				
	X2 = matrix(rnorm(n2*d),n2,d)+ mu  	
	data1 = rbind(X1,X2)
	y = c(rep(1,n1),rep(2,n2))
	DATA = cbind(data1,y)
	DATA
	
}
