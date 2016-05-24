solar.long<-function(year, month, day, time,prec=11)
{
   max=max.day(year,month)
   if(!is.numeric(c(year,month,day,time,prec))|| any(c(month,day,prec)<1) || prec>11 || month>12 || day>max || time<0 || time>=24)
     stop("invalid input parameter(s) specification:check year/month/day/time/precision")
   options(digits=14)
   T<-(ymd2jd(year,month,day)+(time+deltaT(year,month)/3600)/24-2451545)/365250
   L0<-4.8950627+6283.0758500*T-0.0000099*T^2

   table0<-matrix(c(334166,4.669257,6283.075850,3489,4.6261,12566.1517,350,2.744,5753.385,
                    342,2.829,3.523,314,3.628,77713.771,268,4.418,7860.419,234,6.135,3930.210,
                    132,0.742,11506.77,127,2.037,529.691,120,1.110,1577.344,99,5.233,5884.927,
                    90,2.045,26.298,86,3.508,398.149,78,1.179,5223.694,75,2.533,5507.553,
                    51,4.58,18849.23,49,4.21,775.52,36,2.92,0.07,32,5.85,11790.63,28,1.90,796.30,
                    27,0.31,10977.08,24,0.34,5486.78,21,4.81,2544.31,21,1.87,5573.14,20,2.46,6069.78,
                    16,0.83,213.30,13,3.41,2942.46,13,1.08,20.78),ncol=3,byrow=T)
   table1<-matrix(c(20606,430,43,2.67823,2.635,1.59,6283.07585,12566.152,3.52),ncol=3)
   table2<-matrix(c(872,29,1.073,0.44,6283.07585,12566.15),ncol=3)
   table3<-c(29,5.84,6283.07585)

   S0<-0
   for(i in 1:28){
       S0<-S0+table0[i,1]*cos(table0[i,2]+table0[i,3]*T)}
   S1<-table1[1,1]*cos(table1[1,2]+table1[1,3]*T)+table1[2,1]*cos(table1[2,2]+table1[2,3]*T)+
       table1[3,1]*cos(table1[3,2]+table1[3,3]*T)
   S2<-table2[1,1]*cos(table2[1,2]+table2[1,3]*T)+table2[2,1]*cos(table2[2,2]+table2[2,3]*T)
   S3<-table3[1]*cos(table3[2]+table3[3]*T)

   C<-(S0+S1*T+S2*T^2+S3*T^3)*10^(-7)
   round(((L0+C)*180/pi)%%360,prec)
  
}
