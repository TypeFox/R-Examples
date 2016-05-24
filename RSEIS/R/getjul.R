`getjul` <-
function(year, month, day)
{
	jstart = tojul(year, 1, 1);
   
	jul = tojul(year, month, day)-jstart+1;
   
       return(jul)

}

