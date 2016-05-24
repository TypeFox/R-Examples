`mapBreak` <- function(d, ind1, ind2, mr) {
	m = mean(d[ind1:ind2,4])
	if(m>0) {
		if(ind1<length(d[,1]) & ind1>0) { while(d[ind1,4]<mr & ind1) { ind1=ind1+1; if(ind1>length(d[,1]) | ind1<1) { break; } } }
		if(ind1<=length(d[,1]) & ind1>1) { while(d[ind1-1,4]>mr & d[ind1-1,4]>m/2) { ind1=ind1-1; if(ind1>length(d[,1]) | ind1<1) { break; } } }
		if(ind2<=length(d[,1]) & ind2>1) { while(d[ind2,4]<mr) { ind2=ind2-1; if(ind2>length(d[,1]) | ind2<1) { break; } } }
		if(ind2<length(d[,1]) & ind2>0) { while(d[ind2+1,4]>mr & d[ind2+1,4]>m/2) { ind2=ind2+1; if(ind2>length(d[,1]) | ind2<1) { break; } } }
	}
	if(m<0) {
		if(ind1<length(d[,1]) & ind1>0) { while(d[ind1,4]>-mr) { ind1=ind1+1; if(ind1>length(d[,1]) | ind1<1) { break; } } }
		if(ind1<=length(d[,1]) & ind1>1) { while(d[ind1-1,4]< -mr & d[ind1-1,4]< -m/2) { ind1=ind1-1; if(ind1>length(d[,1]) | ind1<1) { break; } } }
		if(ind2<=length(d[,1]) & ind2>1) { while(d[ind2,4]>-mr) { ind2=ind2-1; if(ind2>length(d[,1]) | ind2<1) { break; } } }
		if(ind2<length(d[,1]) & ind2>0) { while(d[ind2+1,4]< -mr & d[ind2+1,4]< -m/2) { ind2=ind2+1; if(ind2>length(d[,1]) | ind2<1) { break; } } }
	}; if(ind2<ind1) { ind2=ind1; }
	return(c(ind1,ind2))
}