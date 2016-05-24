DAYSperYEAR<-function(yr)
{
	n = length(yr)
  YRDAYS = rep(365, n)		
  YRDAYS[yr%%4 == 0 &  yr%%100 != 0 ] =366
  YRDAYS[yr%%4 == 0 &  yr%%100 == 0 & yr%%400 == 0 ] = 366
return(YRDAYS)	
}
