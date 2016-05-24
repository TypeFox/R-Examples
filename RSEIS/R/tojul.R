`tojul` <-
function(year, month, day)
{
## /*   given a year a month and day, return the julian day
##   */ 
       
       yy = year;
       mm = month;
       dd = day;
       jul = 0;


	flg = mm>2
	yy[flg] = yy[flg]
	yy[!flg] = yy[!flg]-1

	mm[flg] = mm[flg] -3
	mm[!flg] = mm[!flg] +9
	
     ##  if(mm>2)
##	{
       ##   yy = yy+0;
##	}
       ##   else
	##{
       ##   yy = yy -1;
##	}
          
      ## if(mm>2)
##	{
        ##  mm = mm -3;
##	}
        ##  else
##	{
      ##    mm = mm+9;
##	}
          
       c = trunc(yy/100);
       ya = yy-100*c;
       jul= trunc((146097*c)/4)+trunc((1461*ya)/4)+ trunc((153*mm+2)/5) +dd+1721119;

	return(jul);

}

