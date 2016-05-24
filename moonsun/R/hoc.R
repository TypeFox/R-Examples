`hoc` <-
function (az,alt,names=NULL) 
{

	res=data.frame(az,alt);
	row.names(res)=names;
	
	class(res$az)="dms";
	class(res$alt)="dms";
	class(res)=c("hoc","apos","data.frame");
	return(res);

}

