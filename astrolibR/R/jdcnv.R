#source('intdiv.R')
jdcnv=function( yr, mn, day, hr) {

  if(max(mn)>12 )  
    warning('month number outside of expected range [1-12] ')
  if(max(day)>31 )
    warning('day number outside of expected range [1-31] ')

  L = intdiv(mn-14,12)	#in leap years, -1 for jan, feb, else 0
  julian = day - 32075 + intdiv(1461*(yr+4800+L),4) + 
    intdiv(367*(mn - 2-L*12),12) - intdiv( 3*intdiv(yr+4900+L,100),4)
  julian = julian + (hr/24.0) - 0.5

  return(julian)
}
