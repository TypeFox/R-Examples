`YRsecdif` <-
function(jd1, hr1, mi1, sec1, jd2, hr2, mi2, sec2, yr1=0, yr2=0)
{
#    subtrace time2 from time1
  if(missing(yr1)) yr1=0
  if(missing(yr2)) yr2=0

  if( all(yr1==yr2) )
    {
      secs1 = jd1*(86400)+hr1*(3600)+mi1*(60)+sec1;
      secs2 = jd2*(86400)+hr2*(3600)+mi2*(60)+sec2;
      return(secs2-secs1)
    }
  else
    {

      EYR = min(c(yr1,yr2)) 

      E1 = EPOCHday(yr1, jd1, origyr=EYR)
      E2 = EPOCHday(yr2, jd2, origyr=EYR)

      secs1 = E1$jday*(86400)+hr1*(3600)+mi1*(60)+sec1;
      secs2 = E2$jday*(86400)+hr2*(3600)+mi2*(60)+sec2;
      
      return(secs2-secs1)
    }

  
 
}

