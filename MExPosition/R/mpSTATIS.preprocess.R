mpSTATIS.preprocess <- function(data, column.design = NULL, row.design = NULL, row.preprocess = 'None', column.preprocess = 'None',  
							           table.preprocess = 'None', make.columndesign.nominal = TRUE, make.rowdesign.nominal = TRUE)
							 
{  column.design <- mpTableCheck(data, t(column.design), make.columndesign.nominal)

   if(is.null(row.design) == FALSE)
   {  row.design <- t(mpTableCheck(t(data), row.design, TRUE))
   }
   
   num.obs = dim(data)[1]
   num.groups = dim(column.design)[1]
  
   # row preprocessing
   row.preproc.data <- mpSTATIS.rowPreproc(data, row.preprocess)

   # column preprocessing
   column.preproc.data <- mpSTATIS.columnPreproc(row.preproc.data$row.processed, column.preprocess)

   # table preprocessing
   table.preproc.data  <- mpSTATIS.tablePreproc(column.preproc.data$column.processed, column.design, table.preprocess)
      
   rownames(table.preproc.data$table.processed) <- rownames(data)
   colnames(table.preproc.data$table.processed) <-colnames(data)

   # results
   preprocess <- list(data = table.preproc.data$table.processed, column.design = column.design, row.design = row.design, data.preprocessed = table.preproc.data$table.processed, num.obs = num.obs,
                       num.groups = num.groups, row.preprocess=row.preprocess, column.preprocess=column.preprocess, table.preprocess=table.preprocess)
	
   print('Preprocessing Completed')
   return(preprocess)
}
