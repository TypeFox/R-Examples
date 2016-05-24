# Provides several methods to binarize a vector consisting of real values.
# The methods are edge-detector-based. With <edge>="firstEdge",
# the first "significant" edge in the sorted data is the threshold for the binarization.
# With the <scaling> factor, the size of the first edge can be adapted.
# The "maxEdge" method searches for the edge with the highest gradient.  
edgeDetector <- function(vector,scaling=1,edge=c("firstEdge","maxEdge"))
{  
  distance<-c()
  binarizeddata<-vector
	
  #sort data
  sortedvector<-sort(vector)
  
  #distance calculation
  for(i in seq_len(length(vector)-1))
  {
    distance[i]<-sortedvector[i+1]-sortedvector[i]  
  }
  
  switch(match.arg(edge),
    #the index of the first edge with distance[i]>threshold is determined
    firstEdge=
      {
        threshold<-scaling*((sortedvector[length(vector)]-sortedvector[1])/((length(vector)-1)))
        index <- 0
        for(i in seq_len(length(vector)-1))
        {  
          if(distance[i]>threshold)
          {
            index<-i
            break  
          }
        }
      },

    #the index of the edge with the maximal gradient is determined
    maxEdge=
      {
       index <- which.max(distance)
      },
    stop("'method' must be one of \"firstEdge\",\"maxEdge\"")
    )
  
  #based on the edge index, the binarization is performed  
  if(index!=0)
  {
    for(i in seq_len(length(vector)))
    {  
      if(vector[i]>=sortedvector[index+1])
        binarizeddata[i]<-1
      else
        binarizeddata[i]<-0
        
    }
    threshold=(sortedvector[index+1]+sortedvector[index])/2
    return(list(bindata=binarizeddata,thresholds=as.numeric(threshold)))
  }
  else
  #if no edge was found, a vector consisting of zeros is returned
  {
  	binarizeddata<-(sapply(vector,function(x) 0))  
    return(list(bindata=binarizeddata,thresholds=NA))
  }  
  
}
