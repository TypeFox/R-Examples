CPdimensionalityplot <-
function(A,n,m,p){

lab=paste(A[,1])	
eff_n=n
if (n>m*p){
	eff_n=m*p
}
eff_m=m
if (m>n*p){
	eff_m=n*p
}
eff_p=p
if (p>m*n){
	eff_p=m*n
}
ep=(eff_n+eff_m+eff_p-2)*A[,1]
PLOT1=plot(A[,1],A[,2],xlab="S",ylab="Fit (%)",pch=20)		
text(A[,1],A[,2],lab,pos=1,cex=0.8)
dev.new()
PLOT2=plot(ep,A[,2],xlab="number of effective parameters",ylab="Fit (%)",pch=20)	
text(ep,A[,2],lab,pos=1,cex=0.8)
}
