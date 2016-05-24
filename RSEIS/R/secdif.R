`secdif` <-
function(jd1, hr1, mi1, sec1, jd2, hr2, mi2, sec2)
{
#    subtrace time2 from time1

      secs1 = jd1*(86400)+hr1*(3600)+mi1*(60)+sec1;
      secs2 = jd2*(86400)+hr2*(3600)+mi2*(60)+sec2;
 
  return(secs2-secs1)
}

