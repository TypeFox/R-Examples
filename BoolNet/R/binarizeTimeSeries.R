# Several binarization methods for time series. <measurements> is a list of matrices with the
# genes in the rows, each specifying one time series.
# If <method> is "kmeans", k-means binarization is used.
# <nstart> and <iter.max> are the corresponding parameters for k-means. This clustering
# is employed for binarization. 
# If <method> is "edgeDetector", an edge detector is used for binarization.
# If <method> is "edgeDetector" and <edge> is "firstEdge",
# the first "significant" edge in the sorted data is the threshold for the binarization.
# If <edge> is "maxEdge", the algorithm searches for the edge with the highest gradient. 
# With the <scaling> factor, the size of the first edge can be adapted.
# If <method> is "scanStatistic", a scan statistic is used to binarize the data
# If <dropInsignificant> is true, insignificant genes (that are not recommended by the statistic) are removed.
# Returns a list containing the binarized matrix list and (for k-means) a vector
# of thresholds and (for scan statistic) a vector specifying genes to remove.
binarizeTimeSeries <- function(measurements, method=c("kmeans","edgeDetector","scanStatistic"), nstart=100, iter.max=1000, edge=c("firstEdge","maxEdge"), scaling=1, windowSize=0.25, sign.level=0.1, dropInsignificant=FALSE)
{
  if (!is.null(dim(measurements)))
    fullData <- measurements
  else
  # in case of list, paste all matrices before clustering
  {
    fullData <- measurements[[1]]
    for (m in measurements[-1])
    {
      fullData <- cbind(fullData,m)
    }
  }
      
  #switch between the different methods
  switch(match.arg(method),
    kmeans={
      cluster <- apply(fullData,1,function(gene)
        # cluster data using k-means
      {
        cl_res <- kmeans(gene, 2, nstart=nstart,iter.max=iter.max)
    
        if (cl_res$centers[1] > cl_res$centers[2])
        # exchange clusters if necessary, so that smaller numbers
        # are binarized to 0, and larger numbers are binarized to 1
          group <- abs(cl_res$cluster-2)
        else
          group <- cl_res$cluster-1
    
           
        # calculate the binarization threshold
        threshold <-  min(cl_res$centers) + dist(cl_res$centers)[1]/2
        list(bin=group, threshold=threshold)
      })
  
  	
      if (is.null(dim(measurements)))
      # split up the collated binarized measurements into a list of matrices of the original size
      {
        startIndex <- 0
        binarizedTimeSeries <- lapply(measurements,function(m)
              {  
                currentSplit <- (startIndex+1):(startIndex+ncol(m))
                startIndex <<- startIndex + ncol(m)
                t(sapply(cluster,function(cl)
                  cl$bin[currentSplit]))          
              })
      }
      else
      {
        binarizedTimeSeries <-   t(sapply(cluster,function(cl)
                cl$bin))
      }
       #colnames(binarizedTimeSeries)<-colnames(fullData)

      return(list(binarizedMeasurements=binarizedTimeSeries,
            thresholds=sapply(cluster,function(cl)cl$threshold)))

    
    },
    edgeDetector={
      #switch between the different edgedetectors
      switch(match.arg(edge),
        firstEdge={
          cluster <- apply(fullData,1,function(gene)
          # cluster data using edgedetector
          {  
            cl_res <- edgeDetector(gene,scaling,edge="firstEdge")
    
            list(bin=cl_res$bindata,thresholds=cl_res$thresholds)
          })
          
          if (is.null(dim(measurements)))
      		# split up the collated binarized measurements into a list of matrices of the original size
      	  {
       			startIndex <- 0
        		binarizedTimeSeries <- lapply(measurements,function(m)
              	{  
                	currentSplit <- (startIndex+1):(startIndex+ncol(m))
                	startIndex <<- startIndex + ncol(m)
                	t(sapply(cluster,function(cl)
                  	cl$bin[currentSplit]))          
              	})
              	threshlist<-sapply(cluster,function(cl) cl$thresholds)
      	  }
      	  else
      	  {
        	 binarizedTimeSeries <-   t(sapply(cluster,function(cl)
             cl$bin))
             threshlist<-sapply(cluster,function(cl) cl$thresholds)
     	  }

      	  
             return(list(binarizedMeasurements=binarizedTimeSeries,thresholds= threshlist))
          },
        maxEdge={
          cluster <- apply(fullData,1,function(gene)
          # cluster data using edgedetector
          {
            cl_res <- edgeDetector(gene,edge="maxEdge")
    
            list(bin=cl_res$bindata,thresholds=cl_res$thresholds)
          })
          
          if (is.null(dim(measurements)))
      		# split up the collated binarized measurements into a list of matrices of the original size
      	   {
       			startIndex <- 0
        		binarizedTimeSeries <- lapply(measurements,function(m)
              	{  
                	currentSplit <- (startIndex+1):(startIndex+ncol(m))
                	startIndex <<- startIndex + ncol(m)
                	t(sapply(cluster,function(cl)
                  	cl$bin[currentSplit]))          
              	})
              	threshlist<-sapply(cluster,function(cl) cl$thresholds)
      	  	}
      		else
      		{
        		binarizedTimeSeries <-   t(sapply(cluster,function(cl)
                cl$bin))
                threshlist<-sapply(cluster,function(cl) cl$thresholds)
     		}

                
          	return(list(binarizedMeasurements=binarizedTimeSeries,thresholds= threshlist))
                  
         },
            
        
        	stop("'method' must be one of \"firstEdge\",\"maxEdge\"")
       	 )
      

    	},
    	scanStatistic={
    		cluster <- apply(fullData,1,function(gene)
          	# cluster data using scanStatistic
          	{  
            	cl_res <- scanStatistic(gene,windowSize,sign.level)
    
            	list(bin= cl_res$bindata,thresholds=cl_res$thresholds,reject=cl_res$reject)
          	})
          	significant <- sapply(cluster,function(cl)(cl$reject==FALSE))
          	# remove not recommended genes    
          	if (dropInsignificant)
          	{	
          		significant <- sapply(cluster,function(cl)(cl$reject==FALSE))
          		cluster <- cluster[significant]
          	
          	}
          
          	if (is.null(dim(measurements)))
      			# split up the collated binarized measurements into a list of matrices of the original size
      	   	{
       			startIndex <- 0
        		binarizedTimeSeries <- lapply(measurements,function(m)
              	{  
                	currentSplit <- (startIndex+1):(startIndex+ncol(m))
                	startIndex <<- startIndex + ncol(m)
                	t(sapply(cluster,function(cl)
                  	cl$bin[currentSplit]))          
              	})
              	rejectlist<-sapply(cluster,function(cl) cl$reject)
              	threshlist<-sapply(cluster,function(cl) cl$thresholds)
      	  	}
      	  	else
      		{
        		binarizedTimeSeries <-   t(sapply(cluster,function(cl)
                cl$bin))
                rejectlist<-sapply(cluster,function(cl) cl$reject)
                threshlist<-sapply(cluster,function(cl) cl$thresholds)
              
     		}
     		if(any(rejectlist))
     		{
     		 	warning("The following genes show a uniform behaviour and should possibly be excluded from binarization:", paste(rownames(fullData)[!significant],collapse=" "))
     		}
     		       		          
                          
        return(list(binarizedMeasurements=binarizedTimeSeries,thresholds=threshlist,reject=rejectlist))
    	
    },
    stop("'method' must be one of \"kmeans\",\"edgeDetector\",\"scanStatistic\"")
    )
    

  }
