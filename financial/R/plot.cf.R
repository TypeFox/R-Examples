"plot.cf" <-
function (x,type=c("bar","npv"),...) 
{
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

cfnfv = function (x,i) 
{
        return(cfnpv(x,i)*spfv(i,length(x)-1))
}


cfnus = function (x,i) 
{
        cfnpv(x,i)/uspv(i,length(x)-1);
}


s <- match.arg(type)
pt <- switch(s, bar =0, npv =1)

if (pt==0) {
	if (is.null(x$tab)) r = 1.1*range(c(x$cf,cumsum(x$cf)))
	else r = 1.1*range(c(x$cf,cumsum(x$cf),as.vector(x$tab[,2])))

	p = barplot(x$cf,ylim=r,...);
	lines(p,cumsum(x$cf),type="b");

	if (!is.null(x$tab)) {
		abline(h=as.vector(x$tab[,2]),lty=1:nrow(x$tab));
		axis(4,x$tab[,2],x$tab[,1]);
	}

}
else if (pt==1)
	{cf = x$cf; 
	curve(cfnpv(cf,x),1,100,xlab="i%",ylab="NPV");
	  abline(h=0);
	  abline(v=x$irr,lty=2);
	  axis(3,x$irr,round(x$irr,2));
	  abline(v=x$ext,lty=3);
	  axis(3,x$ext,round(x$ext,2));
	}



}

