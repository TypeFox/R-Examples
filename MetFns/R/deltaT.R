deltaT<-function(year,month)
{
   if(!is.numeric(c(year,month)) || !(month%in%1:12))
      stop("invalid input parameter specification: check year/month")

   year<-year+(month-0.5)/12
   if(year<1986)
      deltaT<-45.45+1.067*(year-1975)-(year-1975)^2/260-(year-1975)^3/718
      
   else if(year>=1986&year<2005)
           deltaT<-63.86+0.3345*(year-2000)-0.060374*(year-2000)^2+
           0.0017275*(year-2000)^3+0.000651814*(year-2000)^4+0.00002373599*(year-2000)^5

   else if(year>=2005&year<2050)
           deltaT<-62.92+0.32217*(year-2000)+0.005589*(year-2000)^2
           
   deltaT
}
