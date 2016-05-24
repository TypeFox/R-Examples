"intervals" <-
function(X,initboundhandl="reflect"){

#X are data points
#intervals is a vector of values denoting the endpoints of the intervals

#the function gives the interval endpoints

n<-length(X)

order<-order(X)   #in case we need it later
X<-sort(X)
X<-as.row(X)

intervals<-matrix(0,1,n+1)

for (i in 2:n){
	intervals[i] <-(X[i-1]+X[i])/2	#does the middle splitpoints
}

if (initboundhandl=="reflect"){
	intervals[1]<-X[1]-(X[2]-X[1])/2
	intervals[n+1]<-X[n]+(X[n]-X[n-1])/2
}
if (initboundhandl=="stop"){
	intervals[1]<-X[1]
	intervals[n+1]<-X[n]
}

return(intervals)

}
