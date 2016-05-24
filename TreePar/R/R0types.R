R0types<-function(l11,l12,l21,l22,death1,death2){
	L<-l11-l22-death1+death2
	c<-sqrt(L*L+4*l12*l21)
	f1<-(c+L)/(c+L+2*l12)
	R0<-f1*(l11+l12)/death1+(1-f1)*(l22+l21)/death2
	R0
}
