`arrep` <-
function(notation= "arima", phi=c(rep(0,10)), d=0, theta= c(rep(0,10)), Phi=c(rep(0,10)), D=0, Theta= c(rep(0,10)), frequency=1){
	if(!is.null(phi)){ phi<-c(phi,rep(0,(10-length(phi))))}
	if(!is.null(theta)){ theta<-c(theta,rep(0,(10-length(theta))))}
	if(!is.null(Phi)){ Phi<-c(Phi,rep(0,(10-length(Phi))))}
	if(!is.null(Theta)){ Theta<-c(Theta,rep(0,(10-length(Theta))))}

	invert<-function(phi,theta,d,order){
		dif<-c(0,0)
		dif[1]<-ifelse(d>=1,1,0)
		dif[2]<-ifelse(d>=2,1,0)
		if(d!=0 & d!=1 & d!=2) {stop(cat("Differencing value d or D has to be 0,1 or 2 "))}
		#Inverting the non-seasonal component
		ar<-vector(length= order+1)
		ar[1]<- -1
		ar[2]<- ar[1]*-theta[1] + phi[1] + dif[1] + dif[2]
		ar[3]<- ar[2]*-theta[1] + ar[1]*-theta[2] + phi[2] + dif[1]* -phi[1] - dif[2]+ dif[2]* -phi[1]
		ar[4]<- ar[3]*-theta[1] + ar[2]*-theta[2] + ar[1]*-theta[3] + phi[3] + dif[1]* -phi[2] + dif[2]* -phi[2] + dif[2]*phi[1]
		ar[5]<- ar[4]*-theta[1] + ar[3]*-theta[2] + ar[2]*-theta[3] + ar[1]*-theta[4] + phi[4] + dif[1]* -phi[3] + dif[2]* -phi[3] + dif[2]*phi[2]
		ar[6]<- ar[5]*-theta[1] + ar[4]*-theta[2] + ar[3]*-theta[3] + ar[2]*-theta[4] + ar[1]*-theta[5] + phi[5] + dif[1]* -phi[4] + dif[2]* -phi[4] + dif[2]*phi[3]
		ar[7]<- ar[6]*-theta[1] + ar[5]*-theta[2] + ar[4]*-theta[3] + ar[3]*-theta[4] + ar[2]*-theta[5] + ar[1]*-theta[6] + phi[6] + dif[1]* -phi[5] + dif[2]* -phi[5] + dif[2]*phi[4]
		ar[8]<- ar[7]*-theta[1] + ar[6]*-theta[2] + ar[5]*-theta[3] + ar[4]*-theta[4] + ar[3]*-theta[5] + ar[2]*-theta[6] + ar[1]*-theta[7] + phi[7] + dif[1]* -phi[6] + dif[2]* -phi[6] + dif[2]*phi[5]
		ar[9]<- ar[8]*-theta[1] + ar[7]*-theta[2] + ar[6]*-theta[3] + ar[5]*-theta[4] + ar[4]*-theta[5] + ar[3]*-theta[6] + ar[2]*-theta[7] + ar[1]*-theta[8] + phi[8] + dif[1]* -phi[7] + dif[2]* -phi[7] + dif[2]*phi[6]
		ar[10]<- ar[9]*-theta[1] + ar[8]*-theta[2] + ar[7]*-theta[3] + ar[6]*-theta[4] + ar[5]*-theta[5] + ar[4]*-theta[6] + ar[3]*-theta[7] + ar[2]*-theta[8] + ar[1]*-theta[9] + phi[9] + dif[1]* -phi[8] + dif[2]* -phi[8] + dif[2]*phi[7]
		ar[11]<- ar[10]*-theta[1] + ar[9]*-theta[2] + ar[8]*-theta[3] + ar[7]*-theta[4] + ar[6]*-theta[5] + ar[5]*-theta[6] + ar[4]*-theta[7] + ar[3]*-theta[8] + ar[2]*-theta[9] + ar[1]*-theta[10] + phi[10] + dif[1]* -phi[9] + dif[2]* -phi[9] + dif[2]*phi[8]
		ar[12]<- ar[11]*-theta[1] + ar[10]*-theta[2] + ar[9]*-theta[3] + ar[8]*-theta[4] + ar[7]*-theta[5] + ar[6]*-theta[6] + ar[5]*-theta[7] + ar[4]*-theta[8] + ar[3]*-theta[9] + ar[2]*-theta[10] + dif[1]* -phi[10] + dif[2]* -phi[10] + dif[2]*phi[9]
		ar[13]<- ar[12]*-theta[1] + ar[11]*-theta[2] + ar[10]*-theta[3] + ar[9]*-theta[4] + ar[8]*-theta[5] + ar[7]*-theta[6] + ar[6]*-theta[7] + ar[5]*-theta[8] + ar[4]*-theta[9] + ar[3]*-theta[10] + dif[2]*phi[10]
		for(L in 14:(order+1)){
			ar[L]<- ar[L-1]*-theta[1] + ar[L-2]*-theta[2] + ar[L-3]*-theta[3] + ar[L-4]*-theta[4] + ar[L-5]*-theta[5] + ar[L-6]*-theta[6] + ar[L-7]*-theta[7] + ar[L-8]*-theta[8] + ar[L-9]*-theta[9] + ar[L-10]*-theta[10]
		}
		ar
	}
	
	order.sar<-max(trunc(1/(1- abs(Theta[1])^(1/23)))*frequency, 
2*trunc(1/(1- abs(Theta[2])^(1/23)))*frequency, 
3*trunc(1/(1- abs(Theta[3])^(1/23)))*frequency, 
4*trunc(1/(1- abs(Theta[4])^(1/23)))*frequency, 
5*trunc(1/(1- abs(Theta[5])^(1/23)))*frequency, 
6*trunc(1/(1- abs(Theta[6])^(1/23)))*frequency, 
7*trunc(1/(1- abs(Theta[7])^(1/23)))*frequency, 
8*trunc(1/(1- abs(Theta[8])^(1/23)))*frequency, 
9*trunc(1/(1- abs(Theta[9])^(1/23)))*frequency, 
10*trunc(1/(1- abs(Theta[10])^(1/23)))*frequency)
	
	order.ar<-max(trunc(1/(1- abs(theta[1])^(1/23))), 
2*trunc(1/(1- abs(theta[2])^(1/23))), 
3*trunc(1/(1- abs(theta[3])^(1/23))), 
4*trunc(1/(1- abs(theta[4])^(1/23))), 
5*trunc(1/(1- abs(theta[5])^(1/23))), 
6*trunc(1/(1- abs(theta[6])^(1/23))), 
7*trunc(1/(1- abs(theta[7])^(1/23))), 
8*trunc(1/(1- abs(theta[8])^(1/23))), 
9*trunc(1/(1- abs(theta[9])^(1/23))), 
10*trunc(1/(1- abs(theta[10])^(1/23))))
	
	order<-max((max(order.sar, order.ar)/frequency),14)
	
	ar<-invert(phi=phi, theta=theta, d=d, order=trunc(order*frequency))
	sar<-invert(phi=Phi, theta=Theta, d=D, order=trunc(order))
	
	#Combining the non-seasonal and seasonal components into a mixed ar model (mar)
	#First pad the ar series with order*frequency leading 0<U+00ED>s.
	ar.padded<-vector(length= (order*frequency+1 + order*frequency))
	for(i in 1:(order*frequency)){
		ar.padded[i]<-0
	}
	for(i in (order*frequency+1): (order*frequency+1 + order*frequency)){
		ar.padded[i]<-ar[i- (order*frequency)]
	}
	ar.sar<-matrix(nrow= (order*frequency+1), ncol= (order +1))
	mar<-vector(length= (order*frequency+1))
	for(L in 1:(order*frequency+1)){
		for(SL in 1: (order+1)){
			ar.sar[L,SL]<- -ar.padded[(order*frequency) + L-(SL-1)*frequency]*sar[SL]
		}
		mar[L]<-sum(ar.sar[L,])
	}
	seriesna<-vector(length = length(mar))
	for(n in 1: length(mar)){
		seriesna[n]<-ifelse(abs(mar[n])<1.0e-10, NA, n)
	}
	last<-max(seriesna, na.rm=TRUE)
	
	if(notation == "arima") result<- mar[2:last]
	if(notation == "dse1") result<- -mar[1:last]
	return(result)
}

