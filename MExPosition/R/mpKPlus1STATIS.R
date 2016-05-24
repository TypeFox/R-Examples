mpKPlus1STATIS <- function(data, plus1data, column.design, make.columndesign.nominal = TRUE, row.preprocess = 'None', column.preprocess = 'Center', 
                          table.preprocess = 'Sum_PCA', optimization.option = 'STATIS', DESIGN = NULL, make.design.nominal = TRUE, graphs = TRUE)
{	
#########################
# Data preparation
#########################

  main <- deparse(substitute(data))		
  
  if(is.null(data) || is.null(plus1data))
  {   stop("You have not provided the required datasets.")
  }
	
  if(!is.matrix(data))
  {   data1 <- as.matrix(data)
  }

  if(!is.matrix(plus1data))
  {   data2 <- as.matrix(plus1data)
  }
	
  if(sum(is.na(data) > 0))
  {   stop("Missing data in 1st the data matrix; missing data not allowed")
  }

  if(sum(is.na(plus1data) > 0))
  {   stop("Missing data in external matrix; missing data not allowed")
  }
  
  DESIGN <- designCheck(data, DESIGN, make.design.nominal)
  
  res.preproc <- mpSTATIS.preprocess(data = data, column.design = column.design, row.design = NULL, row.preprocess=row.preprocess, column.preprocess=column.preprocess, table.preprocess=table.preprocess, make.columndesign.nominal)
  
  plus1.design <- makeNominalData(matrix(1,1,dim(plus1data)[2]))
  plus1.row.preprocess <- mpSTATIS.rowPreproc(plus1data, row.preprocess)
  plus1.column.preprocess <- mpSTATIS.columnPreproc(plus1.row.preprocess$row.processed, column.preprocess)

  plus1.table.preprocess <- mpSTATIS.tablePreproc(plus1.column.preprocess$column.processed, (plus1.design), table.preprocess)

  res.proc <- mpKPlus1STATIS.core(data = res.preproc$data.preprocessed, plus1data = plus1.table.preprocess$table.processed, num.obs = res.preproc$num.obs, res.preproc$column.design, res.preproc$num.groups, optimization.option = optimization.option)

######################### 
# Results
#########################

KPlus1.statis.overview <- list(data = res.preproc$data.processed, plus1data = plus1.table.preprocess$table.processed, column.design = res.preproc$column.design,
						row.preprocess= res.preproc$row.preprocess, column.preprocess = res.preproc$column.preprocess,
						table.preprocess = res.preproc$table.preprocess, num.groups = res.preproc$num.groups)

KPlus1.statis.innerproduct <- list(S=res.proc$S, S.star=res.proc$S.star, C = res.proc$C, rvMatrix = res.proc$rvMatrix, ci = res.proc$ci, cj = res.proc$cj, 
						eigs.vector = res.proc$eigs.vector, eigs = res.proc$eigs, fi = res.proc$fi, t = res.proc$tau.star, 
            alphaWeights = res.proc$alphaWeights.star)

KPlus1.statis.compromise <- list(compromise = res.proc$compromise, compromise.ci = res.proc$compromise.ci, compromise.cj= res.proc$compromise.cj, 
						compromise.P = res.proc$compromise.P, compromise.eigs = res.proc$compromise.eigs.value, compromise.eigs.vector = res.proc$compromise.eigs.vector,
            compromise.fi = res.proc$compromise.fi, compromise.t = res.proc$compromise.tau)

KPlus1.statis.table <- list(m = res.proc$masses, partial.fi.array = res.proc$table.partial.fi.array, cj = res.proc$table.cj, ci = res.proc$table.ci, 
                        eigs = res.proc$table.eigs,  t = res.proc$table.tau, eigs.vector = res.proc$table.eigs.vector, 
                        Q = res.proc$table.loadings,  fi = res.proc$table.fi,
                        partial.fi = res.proc$table.partial.fi)
      
class(KPlus1.statis.overview) <- c("KPlus1.statis.overview", "list")
class(KPlus1.statis.innerproduct) <- c("KPlus1.statis.innerproduct", "list")
class(KPlus1.statis.compromise) <- c("KPlus1.statis.compromise","list")
class(KPlus1.statis.table) <- c("KPlus1.statis.table","list")

res <- list(Overview = KPlus1.statis.overview, InnerProduct = KPlus1.statis.innerproduct, Compromise = KPlus1.statis.compromise, Table = KPlus1.statis.table)
 
class(res) <- c("mpKPlus1STATIS","list")

print('Processing Complete')

mpPlotInfo <- mpGraphs(res = res, table = res$Overview$column.design, DESIGN = DESIGN, main = main, graphs = graphs)
return(mpOutputHandler(res=res, mpPlotInfo=mpPlotInfo))
}
