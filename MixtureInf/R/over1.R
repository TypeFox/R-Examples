over1 <-
function(alpi,mui,sigi,alpj,muj,sigj)
{
	sigi=sqrt(sigi)
	sigj=sqrt(sigj)

	if(sigi==sigj)
	{
		delta=abs(mui-muj)/sigi
		out=pnorm(-delta/2 + log( alpj / alpi)/delta, 0, 1)
	}
	if(sigi>sigj)
	{
		ncp=(mui-muj)*sigi/(sigi^2-sigj^2)
		value=sigj^2*(mui-muj)^2/(sigj^2-sigi^2)^2-sigj^2/(sigi^2-sigj^2)*log(alpi^2*sigj^2/alpj^2/sigi^2 )
		sqrt.value=sqrt(max(value,0))
		out=pnorm(sqrt.value-ncp,0,1)-pnorm(-sqrt.value-ncp,0,1)
	}
	if(sigi<sigj)
	{
		ncp=(mui-muj)*sigi/(sigi^2-sigj^2)
		value=sigj^2*(mui-muj)^2/(sigj^2-sigi^2)^2-sigj^2/(sigi^2-sigj^2)*log(alpi^2*sigj^2/alpj^2/sigi^2 )
		sqrt.value=sqrt(max(value,0))
		out=1-pnorm(sqrt.value-ncp,0,1)+pnorm(-sqrt.value-ncp,0,1)
	}
	out
}
