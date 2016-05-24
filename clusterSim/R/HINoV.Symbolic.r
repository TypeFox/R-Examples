HINoV.Symbolic<-function(x, u=NULL, distance="H", method = "pam", Index = "cRAND")
{
dist.Symbolic<-function(dane,type="U_2",gamma=0.5,power=2){
indivIC<-dane
individualsNo<-dim(dane)[1]
variablesNo<-dim(dane)[2]
variableSelection<-(1:variablesNo)
distS<-array(0,c(individualsNo,individualsNo))
if (type=="U_2")
{
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
	xmean<-array(0,c(dim(dane)[1],dim(dane)[2]))
	for (i in 1:dim(dane)[1])
	{
		xmean[i,]<-apply(dane[i,,],2,"mean")
		
	}
	resul<-dist(xmean,method="minkowski",p=power)
}
if (type=="H" || type=="S")
{
	zmienne<-0
	zmienne<-dim(dane)[2]
	for (i in 1:(individualsNo-1))
	for (j in (i+1):individualsNo)
	{
		if (type=="S")
		{
			distS[j,i]<-2^(zmienne-1)*(dist(rbind(dane[i,,1],dane[j,,1]))^2+dist(rbind(dane[i,,2],dane[j,,2]^2)))
			distS[i,j]<-distS[j,i]
		}
		
		if (type=="H")
		{
			s<-0
			for (k in 1:zmienne)
			{
				s<-s+max(c(abs(dane[i,k,1]-dane[j,k,1]),abs(dane[i,k,1]-dane[j,k,1])))^2
			}
			distS[i,j]<-distS[j,i]<-s^0.5
		}			
	}
	resul=as.dist(distS)
}

resul
}


	z<-x
	if (is.null(u)) stop ("for symbolic data number of classes must be set")
	#if(!require("cluster")) stop ("Please install cluster package")
	#if(!require("e1071")) stop ("Please install e1071 package")
	#if(!require("ade4")) stop ("Please install ade4 package")
	if (is.null(distance)) stop("For hierarchical methods parameter distance cannot be NULL")
	if (Index != "RAND" && Index!="cRAND") stop("Wrong index type, only RAND or cRAND are allowed")
	if (!is.null(distance))
	{
		if (sum(c("U_2","M","H","S")==distance)==0)
			stop("wrong distance")
	}
	liczba_klas=u
	klasyfikacje<-NULL
	cl<-NULL
	for(i in 1:dim(z)[2])
	{
		x<-(z[,i,])
		#print(i)
		#print(x)
		dim(x)<-c(dim(z)[1],1,2)
		d <- dist.Symbolic(x, type=distance)	
		#print("po dist.symbolic")
		if(method=="pam")
		{
			cl<-pam(d,u,diss=TRUE)$clustering
		}
		else
		{
			cl<-cutree(hclust(d,method=method),u)
		}
		klasyfikacje<-cbind(klasyfikacje,cl)
	}
	liczba_zmiennych=ncol(z)
	wynik<-array(0,c(liczba_zmiennych,liczba_zmiennych))
	for (i in 1:liczba_zmiennych)
	for (j in 1:liczba_zmiennych)
	{
		#print(paste(1,j))
		if(i==j)
		{
			wynik[i,j]=1
		}
		else
		{
	#print("DEBUG: A8")
			#print(klasyfikacje[,i])
			#print(klasyfikacje[,j])
			tk=table(klasyfikacje[,i],klasyfikacje[,j])
	#print("DEBUG: A9")
			w<-classAgreement(tk)
			#print(w)
			if (Index=="cRAND")
			wynik[i,j]<-w$crand
			else
			wynik[i,j]<-w$rand
		}
	}
	#print("DEBUG: A10")
	posortowane<-array(0,c(2,liczba_zmiennych))
	for (i in 1:liczba_zmiennych)
	{
		posortowane[1,i]<-i
		posortowane[2,i]<-sum(wynik[i,])-1
	}
	topri<-posortowane
	if(liczba_zmiennych!=1)
	{
    for (i in 1:liczba_zmiennych)
    for (j in 1:(liczba_zmiennych-1))
    {
      if(posortowane[2,j]<posortowane[2,j+1])
      {
        p1<-posortowane[1,j+1]
        p2<-posortowane[2,j+1]
        posortowane[1,j+1]<-posortowane[1,j]
        posortowane[2,j+1]<-posortowane[2,j]
        posortowane[1,j]<-p1
        posortowane[2,j]<-p2
      }
    }
  }
	resul<-list(parim=wynik,topri=t(topri),stopri=t(posortowane))
	resul
}






