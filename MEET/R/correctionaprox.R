correctionaprox <-function(x,matriu,s) {
	n<-nrow(matriu)
	y<-x-((s-1)/(2*n*log(2,exp(1))))
	y
}

