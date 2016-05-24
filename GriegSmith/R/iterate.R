iterate <-
function(counts,startingdim){


	powers<-c(0:startingdim);
	square<-2^powers;


	x_rects<-sort(c(square,square));
	x_rects<-x_rects[c(-1,-length(x_rects))];

	y_rects0<-2^(1:startingdim);
	y_rects1<-2^(0:(startingdim-1));
	y_rects<-c(rbind(y_rects0,y_rects1));
	
	rects<-rbind(cbind(square,square),cbind(x_rects,y_rects));
	
	
	## rects is a 2 column matrix, the first column is the x length;
	## for each iteration of the G-S method, the second column is the y;
	## width for each iteration. We have both vertically and horizontally;
	## oriented blocks, so we will need to average them. 

	 rects<-rects[order(rowSums(rects)),]

	

	checkhere<-apply(rects,1,sumofsquares,singlecounts=counts);
	mid<-cbind(rects[,1]*rects[,2],rects,checkhere);



	ss<-as.matrix(tapply(mid[,4],mid[,1],mean));
	

	
	ss2<-cbind(ss[-1,1],2*ss[-length(ss),1]);
	blocksize<-as.numeric(rownames(ss2))/2
	rownames(ss2)<-blocksize;





	ssrfinal<-cbind(blocksize,ss2[,2]-ss2[,1],(ss2[,2]-ss2[,1])/(2^(2*startingdim)));


	ssrfinal;

}

