"dojitter" <-
function (x, amount = 0) {
    
#jitters x according to amount

if (amount==0){
	jx<-x
}
else{
	jx<-x + runif(length(x), -amount, amount)
}

jx

}
