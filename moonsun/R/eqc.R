`eqc` <-
function (ra,d,names=NULL) 
{

	res=data.frame(ra,d);
	row.names(res)=names;
	class(res$ra)="time";
	class(res$d)="dms";
	
	class(res)=c("eqc","apos","data.frame");
	return(res);

}

