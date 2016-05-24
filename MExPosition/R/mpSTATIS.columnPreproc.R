mpSTATIS.columnPreproc <- function(data, column.preprocess = 'None') 
{
  column.processed <- matrix(0,dim(data)[1],dim(data)[2])
  processing <- column.processed
    
  ## Exception
  if(column.preprocess != 'None' && column.preprocess != 'Center' && column.preprocess != '1Norm' && 
     column.preprocess != 'Center_1Norm' && column.preprocess != 'Z_Score')
  {    print(paste('WARNING: Column preprocessing option not recognized. Center was set as default'))
       column.preprocess.error = column.preprocess
       column.processed = scale(data, center=TRUE, scale=FALSE)
  }
	
  ## no column preprocessing
  if(column.preprocess == 'None')
  {    column.processed = data
       print(paste('Preprocessed the Columns of the data matrix using: ', column.preprocess))
  }
			
   ## center only
   if(column.preprocess == 'Center')
   {	column.processed = scale(data, center=TRUE, scale=FALSE)
        print(paste('Preprocessed the Columns of the data matrix using: ', column.preprocess))
   }
		
   ## norm = 1 only 
   if(column.preprocess == '1Norm')	
   {	column.processed = data/sqrt(sum((data-colMeans(data))^2))
        print(paste('Preprocessed the Columns of the data matrix using: ', column.preprocess))
   }
			
   ## Center AND norm = 1
   if(column.preprocess == 'Center_1Norm')
   {	processing = scale(data,center=TRUE,scale=TRUE)
        column.processed = processing/sqrt(sum((processing-colMeans(processing))^2))
        print(paste('Preprocessed the Columns of the data matrix using: ', column.preprocess))
   }

   ## Z Score
   if (column.preprocess == 'Z_Score')
   {	 column.processed = scale(data,center=TRUE,scale=TRUE)
         print(paste('Preprocessed the Columns of the data matrix using: ', column.preprocess))
   }
		
   res.column <- list(column.processed = column.processed)
   return(res.column)
}
