"plot.tvm" <-
function (x,row=1,...) 
{
	row = x[row,];

	cf = c(row[3],rep(row[5],row[2]),row[4]);
	cf[2]=cf[2]*(1+row[7]);	
	names(cf)=NULL;
	pl = barplot(cf,ylim=1.1*range(c(cf,cumsum(cf))),...);
	lines(pl,cumsum(cf),type="b");

}

