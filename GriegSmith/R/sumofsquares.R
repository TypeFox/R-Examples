sumofsquares <-
function(sizematrix,singlecounts){

#	print(sizematrix);

	xsublength<-sizematrix[1];
	ysublength<-sizematrix[2];


	xsize<-length(singlecounts[1,]);
	ysize<-length(singlecounts[,1]);


	
	xmin<-rep(seq(from=1,to=xsize,by=xsublength),ysize/ysublength);	
	ymin<-sort(rep(seq(from=1,to=ysize,by=ysublength),xsize/xsublength));

	
	submatricies<-cbind(xmin,ymin);
	
	squaredsums<-sum(apply(submatricies,1,addem,data=singlecounts,xmatlen=xsublength,ymatlen=ysublength));


#	print(squaredsums);

	## sum up all the numbers in each matrix, square those numbers and add them;
	


}

