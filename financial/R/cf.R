"cf" <-
function (x,i=NULL,safe=NULL,rein=safe) 
{
cfnfv = function (x,i) 
{
        return(cfnpv(x,i)*spfv(i,length(x)-1))
}

cfnpv = function (x,i) 
{
        npv = c()       

for (k in 1:length(i)) 
{
        n = length(x)
        j = 0:(n-1)
        pvs = x * sppv(i[k],j)
        npv = c(npv,sum(pvs))
}

return(npv);

}

cfnus = function (x,i) 
{
        cfnpv(x,i)/uspv(i,length(x)-1);
}

cfirr = function (x) 
{
	res = polyroot(x)
	res = Re(res[abs(Im(res))<1e-10]);
	res = (1/res-1)*100;
	return(sort(res));

}

cfext = function (x) 
{
	deriv = x[2:length(x)]
	deriv = deriv*(1:length(deriv));
	res = polyroot(deriv)
	res = Re(res[abs(Im(res))<1e-10]);
	res = (1/res-1)*100;
	return(sort(res));

}


	res = list();
	tab = c();

	if(!is.null(safe)) {

		nn = max(c(length(safe),length(rein)));
		safe=safe+rep(0,nn);
		rein=rein+rep(0,nn);
		modcf=x; modcf[modcf<0]=0;
		nfvp=cfnfv(modcf,rein);
		modcf=x; modcf[modcf>0]=0;
		npvn=-cfnpv(modcf,safe);
		res$mirr=cbind(safe,rein,100*((nfvp/npvn)^(1/(length(x)-1))-1));
		colnames(res$mirr)=c("Safe%","Rein%","MIRR%");
		rownames(res$mirr)=1:nn;
	} else res$mirr=NULL;

	names(x)=1:length(x);
	res$cf=x;
	res$irr=cfirr(x);
	res$ext=cfext(x);

	if (!is.null(i)) {
		for (k in 1:length(i)) {
			tab=rbind(tab,c(i[k],cfnpv(x,i[k]),cfnfv(x,i[k]),cfnus(x,i[k])));
	} 
	rownames(tab)=1:length(i);
	colnames(tab)=c("I%","NPV","NFV","NUS");

	} else tab=NULL;

	res$tab=tab;
	class(res)="cf";
	return(res);

}

