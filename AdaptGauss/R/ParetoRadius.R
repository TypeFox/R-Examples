ParetoRadius <- function(Data ,maximumNrSamples = 1000, plotDistancePercentiles = FALSE){
# MT: in Matlab als ParetoRadiusfuerGMM.m benannt
# ParetoRadius <- ParetoRadius(Data)
# function calculates the paretoRadius for passed gauss mixture modell

# INPUT
# Data[1:n,1:d]				data array, cases in rows, variables in columns
# OPTIONAL
# maximumNrSamples			max number for wich the distance calculation can be done, default 1000
# plotDistancePercentiles	should the plot of percentiles of distances be done? TRUE (default) means yes, FALSE means no

# OUTPUT
# PeratoRadius				the paret radius



#require('dbt.general')

	Data <- as.matrix(Data)

	uxPercentiless <- 18 # % Koordinaten des Punktes mit minimalem unrealisiertem Potential (->Ultsch2001)


	nData <- dim(Data)[1]


	if (maximumNrSamples >= nData){ # no sampling necessary
    	sampleData <- Data
    }else{ #  sample with uniform distribution MaximumNrSamples
		sampleInd <- ceiling(runif(maximumNrSamples, min = 0, max = nData)) # floor(nData*c(runif(maximumNrSamples))+1)
		sampleData <- as.matrix(Data[sampleInd,])
 
	}


# calculate distances
	dd <- as.matrix(dist(sampleData,method='euclidean',diag=TRUE,upper=TRUE))
	matrix1=dd
	size <- dim(matrix1)
	 utri <- upper.tri(matrix1,diag=TRUE)# Upper triangular matrix with 1s and 0s.
    ind <- unname(which(utri==1,arr.ind=TRUE)) # Get indizes of upper triangular matrix.
    rm(utri)
    # Adjust the diagonal.
    ind[,2] <- ind[,2]
    ind <- ind[which(!ind[,2]>size[2]),]
    if(length(ind)==2){ # if ind == 2, R will return a vector instead of a matrix --> convert.
      ind <- t(ind)
    }
    # Create result
    Dist=(matrix1[ind])
	rm(ind)
	rm(matrix1)#
	rm(size)
   	#Dist <- triuvec(dd)	# anneinandergereiht als vektor einschliesslich diagonalen (dist=0)
	
	sdist <- sort(Dist)

    sdl <- length(sdist)
   
    nzdist <- sdist[(dim(sampleData)[1]+1):sdl] 


# selction of ParetoRadius


## etwas abweichendes Ergebnis als bei matlab. grund: matlab hat viele NaNs in nzdist
   # pzt <- percentiles(nzdist)
   x=nzdist
y = c(1:100)
ss<-sort(x)
n<-length(x)
i<-n/100
index<-t(seq(i,n,i))
index<-round(index)
index<-index[1:100]
nullInd<-which(index<1)
index[nullInd]<-1
p<-ss[index]
pzt=p[y]
rm(x)
rm(y)
rm(index)
rm(nullInd)
rm(ss)
rm(n)
rm(p)
rm(i)
#MT: Korrektur ist nicht dasselbe, wieso?
# dd=pdist(sampleData,method="euclidean") #Pairwise distance between observations.
# Dist=sort(squareform(dd)) 
# 
# # selction of ParetoRadius
# nzdist <- percentiles(Dist)  
   
	paretoRadius <- pzt[uxPercentiless]
	
   if (paretoRadius == 0) {
       paretoRadius <-  min(pzt[pzt>0]) # take the smallest nonzero
	}


#    plot of distance distribution

	if(plotDistancePercentiles){
		plot(1:100,pzt,type='l',col='blue', main='red = ParetoRatius',xlab='Percentiles',ylab='Distances')
		lines(x=c(uxPercentiless, uxPercentiless),y=c(0,paretoRadius),col='red')	
	}
#MT: 
#ALUs heuristik, in matlab in PDEplot, hier in dieser Funktion, damit AdaptGauss
# die selbe Darstellung benutzt
if (nData>1024){
  paretoRadius = paretoRadius * 4 /(nData^0.2);
}
	return(paretoRadius)
}
