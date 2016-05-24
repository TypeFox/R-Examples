mlp<-function(inp,weight,dist,neurons,actfns=c(),layer=NaN, ...){

		sigmoid<-function(x) 1/(1+exp(-x));

		tanhip<-function(x) (1-exp(-2*x))/(1+exp(-2*x));

		gauss<-function(x) exp(-(x^2)/2);

		ident<-function(x) x;

		valuate<-function(fnc){
			value<-list();
			value[1]<-list(inp[watch,]);
			for(i in 2:layer){
				ee<-c();
				for(j in 1:neurons[i]){
					e<-0;
					for(k in 1:neurons[i-1]) e<-e+value[[i-1]][k]*weight[[i-1]][k,j];
					ee<-c(ee,fnc[[i-1]](e+dist[[i]][j]));
				}
				value[i]<-list(ee);
			}
			value[[layer]];
		}

	actfns<-as.list(actfns)
	if ((length(neurons)!=length(actfns)+1)&(length(actfns)>0)) return("The number of activation function must be the same as the number of active layer");
	
	if (length(actfns)!=0){
		talal<-FALSE;
		for(i in 1:length(actfns)){
			if (!is.function(actfns[[i]]))
				if ((actfns[[i]]!=1)&(actfns[[i]]!=2)&(actfns[[i]]!=3)&(actfns[[i]]!=4)) talal<-TRUE;
		}
		if (talal) return("Activation functions is a vector and each element of the vector must be between 1-4 or must be a function.");
	}

	if ((is.na(layer))|(layer>=length(neurons))) layer<-length(neurons)
	if (layer<2) return(inp);

	fnc<-c();

	if (length(actfns)>0){
		fnctype<-actfns;
		for(i in 1:(layer-1)){
				if(is.function(fnctype[[i]])) {fnctype[[i]]<-5;fnc<-c(fnc,actfns[[i]])}
				if(fnctype[[i]]==1) fnc<-c(fnc,sigmoid)
				if(fnctype[[i]]==2) fnc<-c(fnc,tanhip)
				if(fnctype[[i]]==3) fnc<-c(fnc,gauss)
				if(fnctype[[i]]==4) fnc<-c(fnc,ident)
			}
	}
	else{
		for(i in 1:(layer-1)) fnc<-c(fnc,sigmoid)
		fnctype<-as.list(rep(1,times=(layer-1)))
	}

	reslt<-c()
	for(watch in 1:nrow(inp))
		reslt<-rbind(reslt,valuate(fnc))
	reslt;
}
