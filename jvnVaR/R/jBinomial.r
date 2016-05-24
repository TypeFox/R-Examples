jBinomial <-
function(n,alpha){
if (n <=0 || alpha <=0 || alpha >=1){
return(0)
}
k <- floor(alpha * (n+1))
if (k<=0){
return((1-alpha)^n)
}
if (k>=n){
return(alpha^n)
}
C <- 1
j <- 0
# n!/k!(n-k)! * alpha^k * (1-alpha)^(n-k)
for (i in 1:k){
	 C <- C * (n-i+1)/(k-i+1)*alpha
	 while (C > 1){
		 C <- C * (1-alpha)
		 j <- j + 1
	 }
}
temp <- n*(1-alpha)^(n-k-j)
C <- C * temp
return(C)
}
