"modjitter" <-function (x,amount){

# fixes jitter to be in correct range so that the test functions can be calculated
# since they take values on [0,1]
# amount now signifies a proportion of grid point distances


d<-1/(length(x)-1)	#finds out how far you can jitter
		  	#since there are length(x)-1 no. intervals

jx<-dojitter(x,amount=amount*d)
jx[1]<-0
jx[length(x)]<-1
if (amount>0){
	jx[2]<-runif(1,0,x[2]+amount*d)
	jx[length(x)-1]<-runif(1,x[length(x)-1]-d*amount,1)
}

jx

}
