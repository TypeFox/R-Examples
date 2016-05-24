`format.dms` <-
function (x,...) 
{
	res = c();

for (ii in 1:length(x)) 
  {

  if (is.finite(x[ii]))
  	{
        h = floor(x[ii]);
        m = (x[ii] %% 1)*60;
        s = round((m %% 1)*60);
        m = floor(m);

        res = c(res,paste(format(h,width=2),"* ",format(m,width=2),"' ",
                        format(s,width=2),"''",sep=""));
      }
	else
	res=c(res,"");
  }
	return(res);
}

