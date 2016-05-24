mpSTATIS.tablePreproc <- function(data, column.design, table.preprocess = 'None') 
{	
   table.processed <- matrix(0,dim(data)[1],dim(data)[2])
   matrixTotal <- sum(data * data)
	
   #table_ids <- matrix(0,dim(column.design)[1],3)
   #colnames(table_ids) <- c('Table_Number', 'From_Column', 'To_Column')
	
   to_total = 0
   from_total = 1
   for(i in 1:dim(column.design)[1])
   {  from = sum(column.design[i-1,]) + from_total
      to = sum(column.design[i,]) + to_total
      to_total = to
      from_total = from
		
      numColumns <- dim(data[,from:to])[2]
	
      ## Table Columns 
     # table_ids[i,1] <- paste('Table',i)
      #table_ids[i,2] <- from
      #table_ids[i,3] <- to
	
      ## Exception
      if (table.preprocess != 'None' && table.preprocess != 'Num_Columns' && table.preprocess != 'Tucker' && table.preprocess != 'Sum_PCA' &&
            table.preprocess != 'RV_Normalization' && table.preprocess != 'MFA_Normalization')
      {   print(paste('WARNING: Table preprocessing option not recognized. Sum PCA was set as default'))
          table.preprocess.error = table.preprocess
          table.processed[,from:to]  <- data[,from:to]/(sqrt(sum(data[,from:to]*data[,from:to])))
      }
	
      ## No whole table preprocessing
      if (table.preprocess == 'None')
      {   table.processed  = data	
      }
	
      ## NumCol Only
      if(table.preprocess == 'Num_Columns')
      {   table.processed[,from:to]  <- data[,from:to]/(numColumns)
      }
	
      ## tucker Only
      if(table.preprocess == 'Tucker')
      {   table.processed[,from:to]  <- data[,from:to]/sqrt(numColumns)
      }
	
      ## SumPCA Only
      if(table.preprocess == 'Sum_PCA')
      {   table.processed[,from:to]  <- data[,from:to]/(sqrt(sum(data[,from:to]*data[,from:to])))
      }
	
      ## RV Only
      if(table.preprocess == 'RV_Normalization')
      { 	table.processed[,from:to] <- data[,from:to]/sqrt(sum(data[,from:to] %*% t(data[,from:to])))
      }
	
      ## MFA Only
      if(table.preprocess == 'MFA_Normalization')
      {   table.processed[,from:to]  <- data[,from:to]/corePCA(data[,from:to])$pdq$Dv[1]
      }
	
   }
	
   ## print
   if (table.preprocess == 'None' || table.preprocess == 'Num_Columns' || table.preprocess == 'Tucker' || table.preprocess == 'Sum_PCA' ||
       table.preprocess == 'RV_Normalization' || table.preprocess == 'MFA_Normalization')
   {	print(paste('Preprocessed the Tables of the data matrix using: ', table.preprocess))
   }

res.table <- list(table.processed = table.processed)
#table.ids = table_ids
return(res.table)
}
