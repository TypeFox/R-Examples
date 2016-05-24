MedianAzimuth <- function (azimuths) 
{
  
  ssize = length(azimuths) 
  his  = Histogram(azimuths, 1)  
  sumle = matrix(0, 360)
  sumri = matrix(0, 360)
  
  for (i in 1:180) {
	for (j1 in i:(179 + i))  {
      sumri[i] <- sumri[i] + his[j1]
    }
	for (k1 in (180 + i):360) {
	  sumle[i] <- sumle[i] + his[k1]
	}
    for (k2 in 1:i)  {
      sumle[i] <- sumle[i] + his[k2]
    }
  }

  for (i in 181:360) {
	for (j1 in i:360)  {
      sumri[i] <- sumri[i] + his[j1]
    }
	for (j2 in 1:(i - 179)) {
	  sumri[i] <- sumri[i] + his[j2]
	}
    for (k3 in (i - 180):i)  {
      sumle[i] <- sumle[i] + his[k3]
    }
  }    
    
	
  dif = abs(sumle - sumri)
  mindif = min(dif)
	
 for (i in 1:360) {
	print(paste(i, sumle[i],"-",sumri[i]))	
	dif = abs(sumle[i] - sumri[i])
	if ( dif <= 1) {
		print(dif)
		median = i
	}
 }
 print(paste ("median =", median, mindif))

  #return(round(median))
}
