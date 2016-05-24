T3dimensionalityplot <-
function(A,n,m,p){

lab=paste(A[,1],A[,2],A[,3], sep = "")	
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
ep=eff_n*A[,1]+eff_m*A[,2]+eff_p*A[,3]+A[,1]*A[,2]*A[,3]-A[,1]^2-A[,2]^2-A[,3]^2
PLOT1=plot(A[,5],A[,4],xlab="P+Q+R",ylab="Fit (%)",pch=20)		
text(A[,5],A[,4],lab,pos=1,cex=0.8)
dev.new()
PLOT2=plot(ep,A[,4],xlab="number of effective parameters",ylab="Fit (%)",pch=20)	
text(ep,A[,4],lab,pos=1,cex=0.8)
}
