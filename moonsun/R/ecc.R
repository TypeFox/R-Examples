`ecc` <-
function (lat,long,names=NULL) 
{

	res=data.frame(lat,long);
	row.names(res)=names;
	class(res$lat)="dms";
	class(res$long)="dms";
	
	class(res)=c("ecc","apos","data.frame");
	return(res);

}

