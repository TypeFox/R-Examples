RegCircLin<-function(x,y,t,bw,method,tol){
	if (method=="NW"){
		t_x<-outer(t,x,"-")
		if (max(bw)>tol){
			z<-apply(bw*cos(t_x),1,max)-tol
			B<-exp(bw*cos(t_x)-z) 
		}else {
			B<-exp(bw*cos(t_x))
		}
		L<-B/apply(B,1,sum)
	}else{
		x_t<-outer(x,t,"-")
		if (max(bw)>tol){
			z<-apply(bw*cos(x_t),2,max)-tol
			kxt<-exp(t(t(bw*cos(x_t))-z)) 
		}else {
			kxt<-exp(bw*cos(x_t))
		}
		m<-sin(x_t)
		sn1<-apply(kxt*m,2,sum)
		sn2<-apply(kxt*m^2,2,sum)
		bjti<-t(kxt)*(sn2-t(m)*sn1)
		L<-bjti/apply(bjti,1,sum)
	}
	fhat<-as.vector(L%*%y)
	return(fhat)
}
