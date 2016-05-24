max.day<-function(year,month)
{
  if(!is.numeric(c(year,month)) || !(month%in%1:12) || !(year%in%1984:2100) ) 
     stop("invalid input parameter specification: check year/month")
  
  if(month%in%c(1,3,5,7,8,10,12))
     {max=31}
  else if(month==2)
     {max=ifelse(year%%4==0,29,28)}
  else {max=30}
 
  max

}
