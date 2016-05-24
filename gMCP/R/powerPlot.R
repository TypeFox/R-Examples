powerplot <- function() {
	x1<-seq(-10,10,length=41)
	x2<-seq(-10,10,length=41)
	
	f<-function(x1, x2){
		return(x1*x2)
	}
	
	z<-outer(x1, x2, f)
	
	persp(x1, x2, z, 
			main="Power", 
			col="lightgreen", 
			theta=30, phi=20, 
			r=50, 
			d=0.1, 
			expand=0.5, 
			ltheta=90, lphi=180, 
			shade=0.75, 
			ticktype="detailed", 
			nticks=5)
}