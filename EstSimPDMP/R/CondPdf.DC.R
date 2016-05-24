CondPdf.DC <-
function( dat , x , t , h=NULL , alpha=1/5 , bound=Inf){
	s<-0
	
	N<-length(dat[1,])-1; m<-length(dat[,1]);	
	A<-dat[,1:N]; B<-as.matrix(dat[ !duplicated(A) ,1:N]);
	
	for (k in 1:length(B[,1])){
		y<-c();
		for (l in 1:N){
			y<-c(y , B[k,l])
		}
		D<-.Tri2( dat , x , y)
		if (length(D)>0){
			s<-s+HR(D , t , h , alpha , bound)*.CondSurv(dat , x , y , t)
			}
		}
	list(time=t,pdf=s)	
}
