ShapleyShubik <-
function(quota,y,Names=NULL){
	
	
	n <- length(y)
	res1 <- permutations(n, n)
	res2 <- apply(res1, 1, 
	function(x) {
    x[sum(cumsum(y[x]) <quota)+1 ] })	

res2<-as.numeric(res2)

Power<-matrix(NA,ncol=1,nrow=n)
for (i in 1:n){
Power[i,1]<-sum(res2==i)
}


SHI<-Power/factorial(n)
TABLE<-rbind(y,y/sum(y),t(SHI))
rownames(TABLE)<-c("Votes","Votes (%)", "Shapley-Shubik")
colnames(TABLE)<-Names

Output<-list(Results=TABLE,Distribution=y,C,Method="PowerIndex",Quota=quota,Names=Names)
class(Output)<-"ShapleyShubik"
return(Output)

}
