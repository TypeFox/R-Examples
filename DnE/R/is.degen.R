is.degen <-
function(x)
{
	if(var(x)==0)
	{
		return(data.frame("qchisq"=Inf,"pvalue"=0));
	}
	else
	{
		return(data.frame("state"=-1,"pvalue"=1));
	}
}
