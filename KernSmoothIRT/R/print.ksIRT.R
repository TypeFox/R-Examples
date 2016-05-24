print.ksIRT <-
function(x,...){
	toout<-data.frame(Item=1:x$nitem,Correlation=x$itemcor)
	print(toout)
}

