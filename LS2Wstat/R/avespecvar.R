avespecvar <-function (spectrum) {

m <- spectrum$S
m <- matrix(m,nrow=dim(m)[1])

statistic<-mean(rowVars(m))

return(statistic)

}

