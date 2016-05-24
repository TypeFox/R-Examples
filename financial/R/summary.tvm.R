"summary.tvm" <-
function (object,row=1,...) 
{

cat("\nAmortization Table\n\n");

	x=object;
	row=x[row,];
	n = row[2];  a=row[7];
	i = row[1]/(100*row[8]);
	pv = row[3]; fv = row[4]; pmt=row[5];
	days = row[6]; pyr=row[8];

	bal=pv+a*pmt; 
	res=c()

	for (k in 1:(n-a))
	{
	if (k==1) {
		int = bal*i*(days/(360/pyr));
		prin = pmt+int;
		bal = bal+prin;
		prin = prin+a*pmt;
		res = rbind(res,c(bal,int,prin,pmt*(1+a)));
	}
	else
	{
		int = bal*i;
		prin = pmt+int;
		bal = bal+prin;
		res = rbind(res,c(bal,int,prin,pmt));
	}
	}
	res=rbind(res,c(NA,sum(res[,2]),sum(res[,3]),sum(res[,4])));
	colnames(res)=c("Bal","Int","Prin","PMT");
	rownames(res)=c(1:(n-a),"Total");
	print(round(res,2),na.print="");
	invisible(res);
}

