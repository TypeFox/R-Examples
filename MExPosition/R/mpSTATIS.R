mpSTATIS <- function(data, column.design, make.columndesign.nominal = TRUE, row.design = NULL, make.rowdesign.nominal = FALSE, statis.prepro.option = 'Plain_STATIS', DESIGN = NULL, make.design.nominal = TRUE, graphs = TRUE)
{	
#########################
# Data preparation
#########################

  main <- deparse(substitute(data))		
  
  if(is.null(data))
  {   stop("You have not provided any data.")
  }
	
  if(!is.matrix(data))
  {   data <- as.matrix(data)
  }
	
  if(sum(is.na(data) > 0))
  {   stop("Missing data not allowed")
  }
  
  DESIGN <- designCheck(data, DESIGN, make.design.nominal)
  
	statis.option.list <- c('Plain_STATIS','MFA','Sum_PCA','Plain_Multitable','ANISOSTATIS_Type1','ANISOSTATIS_Type2','PTA','COVSTATIS','Power_STATIS','Customization')
	optimization.option.list <- c('None', 'STATIS', 'RV_Matrix', 'STATIS_Power1', 'ANISOSTATIS_Type1', 'ANISOSTATIS_Type2')
	table.option.list <- c('None', 'Num_Columns', 'Tucker', 'Sum_PCA', 'RV_Normalization', 'MFA_Normalization')
	row.option.list <- c('None', 'Center', 'Hellinger', 'Center_Hellinger','Profile')
	col.option.list <- c('None', 'Center', '1Norm', 'Center_1Norm', 'Z_Score')  

#########################
# Preset STATIS Options
#########################
	 if(!statis.prepro.option %in% statis.option.list || statis.prepro.option=='Plain_STATIS')
	 {	if(!statis.prepro.option %in% statis.option.list)
	 	{	print('WARNING: Preset STATIS Option not recognized. Plain STATIS was set as default')
		}
		row.selection <- 'None'
		column.selection <- 'Center_1Norm'
		table.selection <- 'Sum_PCA'
		optimization.selection <- 'STATIS'
		
		##Customization selected		
		}
		else if(statis.prepro.option=='Customization')
		{	print('You have selected Customization. Please choose from the follow options by selecting a number.')
			keep.row.going <- TRUE
			while(keep.row.going)
			{	cat('Row Preprocessing Selection\n')
				cat(paste(paste(1:length(row.option.list),row.option.list,sep=": "),"\n",sep=""))
				valueCaptured <- as.numeric(readline())
				if(valueCaptured%in%1:length(row.option.list))
				{	row.selection <- row.option.list[valueCaptured]
					keep.row.going <- FALSE
				}
			}
			keep.col.going <- TRUE
			while(keep.col.going)
			{	cat('Column Preprocessing Selection\n')
				cat(paste(paste(1:length(col.option.list),col.option.list,sep=": "),"\n",sep=""))
				valueCaptured <- as.numeric(readline())
				if(valueCaptured%in%1:length(col.option.list))
				{	column.selection <- col.option.list[valueCaptured]
					keep.col.going <- FALSE
				}
			}	
			keep.table.going <- TRUE
			while(keep.table.going)
			{	cat('Table Preprocessing Selection\n')
				cat(paste(paste(1:length(table.option.list),table.option.list,sep=": "),"\n",sep=""))
				valueCaptured <- as.numeric(readline())
				if(valueCaptured%in%1:length(table.option.list))
				{	table.selection <- table.option.list[valueCaptured]
					keep.table.going <- FALSE
				}
			}	
			keep.optimization.going <- TRUE
			while(keep.optimization.going)
			{	cat('Optimization Selection\n')
				cat(paste(paste(1:length(optimization.option.list),optimization.option.list,sep=": "),"\n",sep=""))
				valueCaptured <- as.numeric(readline())
				if(valueCaptured%in%1:length(optimization.option.list))
				{	optimization.selection <- optimization.option.list[valueCaptured]
					keep.optimization.going <- FALSE
				}
			}	
		}
	#this will get replaced soon with functions for each of these...
	#Select STATIS type by statis.prepro
	else{
		if(statis.prepro.option=='MFA'){
			row.selection <- 'None'
			column.selection <- 'Center_1Norm'
			table.selection <- 'MFA_Normalization'
			optimization.selection <- 'None'	
		}
		if(statis.prepro.option=='Sum_PCA'){
			row.selection <- 'None'
			column.selection <- 'Center_1Norm'
			table.selection <- 'Sum_PCA'
			optimization.selection <- 'None'		
		}	
		if(statis.prepro.option=='Plain_Multitable'){
			row.selection <- 'None'
			column.selection <- 'None'
			table.selection <- 'None'
			optimization.selection <- 'None'	
		}
		if(statis.prepro.option=='ANISOSTATIS_Type1'){
			row.selection <- 'None'
			column.selection <- 'Center_1Norm'
			table.selection <- 'Sum_PCA'
			optimization.selection <- 'ANISOSTATIS_Type1'			
		}
		if(statis.prepro.option=='ANISOSTATIS_Type2'){
			row.selection <- 'None'
			column.selection <- 'Center_1Norm'
			table.selection <- 'Sum_PCA'
			optimization.selection <- 'ANISOSTATIS_Type2'			
		}
		if(statis.prepro.option=='PTA'){
			row.selection <- 'None'
			column.selection <- 'Center_1Norm'
			table.selection <- 'Sum_PCA'
			optimization.selection <- 'STATIS'			
		}

	}
	res.preproc <- mpSTATIS.preprocess(data, column.design, row.design = NULL, row.preprocess=row.selection, column.preprocess=column.selection, table.preprocess=table.selection, make.columndesign.nominal)
	res.proc <- mpSTATIS.optimize(res.preproc$data.preprocessed, res.preproc$num.obs, res.preproc$column.design, res.preproc$num.groups,optimization.option=optimization.selection)
	
######################### 
# Results
#########################

statis.overview <- list(data = res.preproc$data, groupmatrix = res.preproc$column.design, preprocess.data = res.preproc$data.preprocessed, 
						num.groups = res.preproc$num.groups, num.obs = res.preproc$num.obs, row.preprocess = res.preproc$row.preprocess, 
						column.preprocess = res.preproc$column.preprocess, table.preprocess = res.preproc$table.preprocess)

statis.innerproduct <- list(S=res.proc$S, C = res.proc$C, RVMatrix = res.proc$RVMatrix, ci = res.proc$ci, cj = res.proc$cj, eigs.vector = res.proc$eigs.vector, 
						eigs = res.proc$eigs, fi = res.proc$fi, t = res.proc$tau, alphaWeights = res.proc$alphaWeights)

statis.compromise <- list(compromise = res.proc$compromise, compromise.ci = res.proc$compromise.ci, compromise.cj = res.proc$compromise.cj, 
						compromise.eigs.vector = res.proc$compromise.eigs.vector, compromise.eigs = res.proc$eigs, compromise.fi = res.proc$compromise.fi,
						compromise.t = res.proc$compromise.tau)

statis.table <- list(m = res.proc$masses, eigs = res.proc$table.eigs, eigs.vector = res.proc$table.eigs.vector, 
            		Q = res.proc$table.loadings, fi = res.proc$table.fi, partial.fi = res.proc$table.partial.fi,
         		    cj = res.proc$table.cj, ci = res.proc$table.ci, t =res.proc$table.tau, partial.fi.array = res.proc$table.partial.fi.array)  
       
class(statis.overview) <- c("statis.overview", "list")
class(statis.innerproduct) <- c("statis.innerproduct", "list")
class(statis.compromise) <- c("statis.compromise","list")
class(statis.table) <- c("statis.table","list")

res <- list(Overview = statis.overview, InnerProduct = statis.innerproduct, Compromise = statis.compromise, Table = statis.table)
 
class(res) <- c("mpSTATIS","list")

print('Processing Complete')

mpPlotInfo <- mpGraphs(res = res, table = res$Overview$groupmatrix, DESIGN = DESIGN, main = main, graphs = graphs)
return(mpOutputHandler(res=res, mpPlotInfo=mpPlotInfo))
}
