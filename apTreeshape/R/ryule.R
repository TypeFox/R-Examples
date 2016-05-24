ryule <-  function(n){

	lst<-(-n):(-1)

	Data <- matrix(NA, nrow = n-1, ncol = 2)
	X <- NULL
	S <- 0
	for (i  in 1:(n-1)){
		Data[i, 1:2] <- xx <- sort(sample(lst, 2,replace= FALSE))
		lst <- c(lst[ (lst !=xx[1]) & (lst != xx[2]) ],i)
	}
	treeshape(Data)
}
