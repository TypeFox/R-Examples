mpSTATIS.rowPreproc <-function(data, row.preprocess = 'None') 

{
    row.processed <- matrix(0,dim(data)[1],dim(data)[2])
    processing <- row.processed
		
		## Exception
		if(row.preprocess != 'None' && row.preprocess != 'Profile' && row.preprocess != 'Hellinger' && row.preprocess != 'Center' && row.preprocess != 'Center_Hellinger')
		{	print(paste('WARNING: Row preprocessing option not recognized. No row preprocessing was set as default'))
			row.preprocess.error = row.preprocess
			row.processed = data
		}
		
		## No row preprocessing
		if (row.preprocess == 'None') 
		{	row.processed = data
			print(paste('Preprocessed the Rows of the data matrix using: ', row.preprocess))
		}
	
		## RowProfile only
		if (row.preprocess == 'Profile')
		{	row.processed <- diag(as.vector((data %*% repmat(1,dim(data)[2],1))^-1)) %*% data
			print(paste('Preprocessed the Rows of the data matrix using: ', row.preprocess))
		}
		
		## Hellinger only
		if (row.preprocess == 'Hellinger')
		{	row.processed <- data/sqrt(rowSums(data %*% t(data)))
			print(paste('Preprocessed the Rows of the data matrix using: ', row.preprocess))
		}
		
		## Row Center only
		if (row.preprocess == 'Center')
		{	row.processed <- t(scale(t(data),center = TRUE, scale=FALSE))
			print(paste('Preprocessed the Rows of the data matrix using: ', row.preprocess))
		}
			
		## Hellinger AND Center
		if (row.preprocess == 'Center_Hellinger')	
		{	processing = data/sqrt(rowSums(data %*% t(data)))
			row.processed = t(scale(t(processing),center = TRUE, scale=FALSE))
			print(paste('Preprocessed the Rows of the data matrix using: ', row.preprocess))
		}
		
res.row <- list(row.processed = row.processed)
return (res.row)
}
