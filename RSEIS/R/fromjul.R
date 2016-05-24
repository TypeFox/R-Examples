`fromjul` <-
function(jul, yy)
{
       j=jul-1721119
       yy=trunc((4*j-1)/146097)
       j=4*j-1-146097*yy
       dd=trunc(j/4)
       j=trunc((4*dd+3)/1461)
       dd=4*dd+3-1461*j
       dd=trunc((dd+4)/4)
       mm=trunc((5*dd-3)/153)
       dd=5*dd-3-153*mm
       dd=trunc((dd+5)/5)
       yy=100*yy+j


yy[mm<10]=yy[mm<10]
yy[mm>=10]=yy[mm>=10]+1;
     ##  if(mm<10)
   ## 	{	
       ##       yy=yy+0;
   ## 	}
      ##        else
   ## 	{
      ##        yy=yy+1;
   ## 	}

	flg = mm<10


	mm[flg] = mm[flg]+3;  
	mm[!flg] =  mm[!flg]-9;
     	##  if(mm<10)
	##	{
      	##    mm= mm+3;
	##	}
       	##   else
	##	{
       	##   mm=mm-9;
	##	}



return(list(mo=mm, dom=dd))	
}

