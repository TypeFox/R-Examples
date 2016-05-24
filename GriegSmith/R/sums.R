sums <-
function (coords,dim,xmin=min(coords[,1]),xmax=max(coords[,1]),ymin=min(coords[,2]),ymax=max(coords[,2])){


	xints<-((xmax-xmin)/dim);
	yints<-((ymax-ymin)/dim);

	xbins<-seq(from=xmin, to=xmax-xints, by=xints);
	ybins<-seq(from=ymin, to=ymax-yints, by=yints);



	bins<-cbind(c(sapply(xbins,rep,dim)), rep(ybins,dim));
	cnts<-matrix(apply(bins,1,belongtoint,vect=coords,int.x=xints,int.y=yints),nrow=dim,byrow=TRUE);


}

