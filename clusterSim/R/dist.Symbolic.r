dist.Symbolic<-function(data,type="U_2",gamma=0.5,power=2){
indivIC<-data
individualsNo<-dim(data)[1]
variablesNo<-dim(data)[2]
variableSelection=(1:variablesNo)
distS<-array(0,c(individualsNo,individualsNo))
if (type=="U_2"){
  for (i in 1:(individualsNo-1))
  for (j in (i+1):individualsNo)
  {
    D<-0
    for(k in variableSelection)
    {
       d<-(max(indivIC[i,k,2],indivIC[j,k,2])-min(indivIC[i,k,1],indivIC[j,k,1]))
       if (!((indivIC[i,k,2]<=indivIC[j,k,1]) || (indivIC[j,k,2]<=indivIC[i,k,1])))
      {
         if (indivIC[i,k,1]<=indivIC[j,k,1])
        {
          if (indivIC[i,k,2]<=indivIC[j,k,2])
            d<-d-(1-2*gamma)*abs(indivIC[i,k,2]-indivIC[j,k,1])
          else
            d<-d-(1-2*gamma)*abs(indivIC[j,k,2]-indivIC[j,k,1])
        }
        else
        {
          if (indivIC[j,k,2]<=indivIC[i,k,2])
            d<-d-(1-2*gamma)*abs(indivIC[j,k,2]-indivIC[i,k,1])
          else
            d<-d-(1-2*gamma)*abs(indivIC[i,k,2]-indivIC[i,k,1])
        }
      } 
      d<-d-gamma*abs(indivIC[j,k,2]-indivIC[j,k,1])
      d<-d-gamma*abs(indivIC[i,k,2]-indivIC[i,k,1])
    }
    D<-D+d^power
    distS[i,j]<-D^(1/power)		
    distS[j,i]<-distS[i,j]
  }
  resul<-as.dist(distS)
}
if (type=="M")
{
	xmean<-array(0,c(dim(data)[1],dim(data)[2]))
	for (i in 1:dim(data)[1])
	{
		xmean[i,]<-apply(data[i,,],2,"mean")
		
	}
	resul<-dist(xmean,method="minkowski",p=power)
}
if (type=="H" || type=="S")
{
	zmienne<-0
	zmienne<-dim(data)[2]
	for (i in 1:(individualsNo-1))
	for (j in (i+1):individualsNo)
	{
		if (type=="S")
		{
			distS[j,i]<-2^(zmienne-1)*(dist(rbind(data[i,,1],data[j,,1]))^2+dist(rbind(data[i,,2],data[j,,2]^2)))
			distS[i,j]<-distS[j,i]
		}
		
		if (type=="H")
		{
			s<-0
			for (k in 1:zmienne)
			{
				s<-s+max(c(abs(data[i,k,1]-data[j,k,1]),abs(data[i,k,1]-data[j,k,1])))^2
			}
			distS[i,j]<-distS[j,i]<-s^0.5
		}			
	}
	resul<-as.dist(distS)
}
resul
}


