rbf<-function(inp,weight,dist,neurons,sigma, ...){

		gauss<-function(x,sigma) exp(-(x^2)/(2*sigma^2));

		ident<-function(x) x;


		valuate3<-function(pont){
			value<-vector(2,mode="list");
			value[[1]]<-pont;
			total<-0;
			ee<-c();


			for(j in 1:neurons[2]){
				e<-0;
				for(k in 1:neurons[1]) 
					e<-e+abs(value[[1]][k]-weight[[1]][k,j]);
				total<-total+gauss(e,sigma[k,j]);
			}


			for(j in 1:neurons[2]){
				e<-0;
				for(k in 1:neurons[1]) 
					e<-e+abs(value[[1]][k]-weight[[1]][k,j]);
				ee<-c(ee,gauss(e,sigma[k,j])/total);
			}
			value[[2]]<-ee;

			e<-0;
			for(k in 1:neurons[2]) 
				e<-e+value[[2]][k]*weight[[2]][k,1];
			ident(e+dist[1]);
		}

	value<-c();
	for(i in 1:nrow(inp)) value<-rbind(value,valuate3(inp[i,]));
	value;
}
