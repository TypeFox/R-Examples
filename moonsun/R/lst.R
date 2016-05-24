`lst` <-
function (...,lambda=getOption("longitude")) 
{
	gstime = gst(...);

	if (is.null(lambda)) 
			{
			lambda = 0;
			warning("Your longitude is not set in the environment, assuming it is equal to 0");
			}

	gstime = gstime + lambda/15;
	gstime = gstime %% 24;

	class(gstime)=c("lst","time");
	return(gstime);
	
}

