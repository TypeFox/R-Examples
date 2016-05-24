GriegSmith <-
function(datapoints,counts=FALSE,env=100){


	if(counts==FALSE){

	if(is.ppp(datapoints)){

		xmin<-datapoints$window$xrange[1]
		xmax<-datapoints$window$xrange[2]

		ymin<-datapoints$window$yrange[1]
		ymax<-datapoints$window$yrange[2]


		datapoints<-cbind(datapoints$x,datapoints$y);


		
	}
	else{

		xmax<-max(datapoints[,1]);
		xmin<-min(datapoints[,1]);

		ymax<-max(datapoints[,2]);
		ymin<-min(datapoints[,2]);





	}


		
	
	numpts<-length(datapoints[,1]);
	startingdim<-ceiling(log(numpts)/(2*log(2)));
	counts<-sums(datapoints,2^startingdim,xmin,xmax,ymin,ymax)
		

	}
	else {
		if (max(datapoints[,1]) != max(datapoints[,2])) stop("Your count data must have equal dimensions")


		datapoints<-datapoints[order(datapoints[,2],datapoints[,1]),]
		numpts<-sum(datapoints[,3])
		startingdim<-ceiling(log(max(datapoints[,1]))/log(2))
		counts<-matrix(datapoints[,3],nrow=2^startingdim,byrow=TRUE);
	
		
	}	


	actual<-iterate(counts,startingdim);
	sims<-envelopes(env,counts,startingdim);
	final<-cbind(actual,sims);
	

	colnames(final)<-c("blocksize","SSr","MSr","MSr.05","MSr.95");

	
	class(final) <- "GriegSmith"


	final;

	


}

