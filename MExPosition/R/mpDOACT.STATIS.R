mpDOACT.STATIS <- function(data1, column.design.1, make.columndesign.1.nominal = TRUE, data2, column.design.2, make.columndesign.2.nominal = TRUE, 
						   row.preprocess.data1 = 'None', column.preprocess.data1 = 'Center', table.preprocess.data1 = 'Sum_PCA',  
						   row.preprocess.data2 = 'None', column.preprocess.data2 = 'Center', table.preprocess.data2 = 'Sum_PCA', DESIGN = NULL, 
						   make.design.nominal = TRUE, graphs = TRUE)
{	
#########################
# Data preparation
#########################


  main <- deparse(substitute(data1))		
  
  if(is.null(data1) || is.null(data2))
  {   stop("You have not provided two datasets.")
  }
	
  if(!is.matrix(data1))
  {   data1 <- as.matrix(data1)
  }

  if(!is.matrix(data2))
  {   data2 <- as.matrix(data2)
  }
	
  if(sum(is.na(data1) > 0))
  {   stop("Missing data in 1st data matrix; missing data not allowed")
  }

  if(sum(is.na(data2) > 0))
  {   stop("Missing data in 2nd data matrix; missing data not allowed")
  }
  
  DESIGN <- designCheck(data1, DESIGN, make.design.nominal)
  
  res.preproc.data.1 <- mpSTATIS.preprocess(data = data1, column.design = column.design.1, row.design = NULL, row.preprocess=row.preprocess.data1, column.preprocess=column.preprocess.data1, table.preprocess=table.preprocess.data1, make.columndesign.1.nominal)
  res.preproc.data.2 <- mpSTATIS.preprocess(data = data2, column.design = column.design.2, row.design = NULL, row.preprocess=row.preprocess.data2, column.preprocess=column.preprocess.data2, table.preprocess=table.preprocess.data2, make.columndesign.2.nominal)
  res.proc <- mpDOACT.STATIS.core(res.preproc.data.1$data.preprocessed, res.preproc.data.1$column.design, res.preproc.data.2$data.preprocessed, res.preproc.data.2$column.design)
	
######################### 
# Results
#########################

doact.statis.overview <- list(data1 = res.preproc.data.1$data, data2 = res.preproc.data.2$data, column.design.1 = res.preproc.data.1$column.design, column.design.2 = res.preproc.data.2$column.design,
						row.preprocess.data1 = res.preproc.data.1$row.preprocess, column.preprocess.data1 = res.preproc.data.1$column.preprocess,
						table.preprocess.data1 = res.preproc.data.1$table.preprocess, row.preprocess.data2 = res.preproc.data.1$row.preprocess, 
						column.preprocess.data2 = res.preproc.data.2$column.preprocess, table.preprocess.data2 = res.preproc.data.2$table.preprocess, num.groups.1 = res.preproc.data.1$num.groups,
						num.groups.2 = res.preproc.data.2$num.groups)

doact.statis.innerproduct <- list(S.1=res.proc$S.1, S.2=res.proc$S.2, C = res.proc$C, RVMatrix = res.proc$RVMatrix, ci = res.proc$ci, cj = res.proc$cj, 
						eigs.vector = res.proc$eigs.vector, eigs = res.proc$eigs, fi = res.proc$fi, t = res.proc$tau, alphaWeights = res.proc$alphaWeights,
						betaWeights = res.proc$betaWeights)

doact.statis.compromise <- list(compromiseMatrix.1 = res.proc$compromiseMatrix.1, compromise.ci.1 = res.proc$compromise.ci.1, compromise.cj.1 = res.proc$compromise.cj.1, 
						compromise.eigs.vector.1 = res.proc$compromise.P.1, compromise.eigs.value.1 = res.proc$compromise.eigs.value.1, compromise.fi.1 = res.proc$compromise.fi.1, 
						compromise.t.1 = res.proc$compromise.tau.1, compromiseMatrix.2 = res.proc$compromiseMatrix.2, compromise.ci.2 = res.proc$compromise.ci.2, 
						compromise.cj.2 = res.proc$compromise.cj.2, compromise.eigs.vector.2 = res.proc$compromise.P.2,  compromise.eigs.value.2 = res.proc$compromise.dd.2, 
						compromise.fi.2 = res.proc$compromise.fi.2, compromise.t.2 = res.proc$compromise.tau.2)

doact.statis.table <- list(m.1 = res.proc$masses.1, partial.fi.array.1 = res.proc$table.partial.fi.array.1, cj.1 = res.proc$table.cj.1, ci.1 = res.proc$table.ci.1, 
                        eigs.1 = res.proc$table.eigs.1,  t.1 = res.proc$table.tau.1, eigs.vector.1 = res.proc$table.eigs.vector.1, 
                        Q.1 = res.proc$table.loadings.1,  fi.1 = res.proc$table.fi.1, partial.fi.1 = res.proc$table.partial.fi.1, tau.1 = res.proc$table.tau.1, 

                        m.2 = res.proc$masses.2, partial.fi.array.2 = res.proc$table.partial.fi.array.2, cj.2 = res.proc$table.cj.2, ci.2 = res.proc$table.ci.2, 
                        eigs.2 = res.proc$table.eigs.2,  t.2 = res.proc$table.tau.2, eigs.vector.2 = res.proc$table.eigs.vector.2, 
                        Q.2 = res.proc$table.loadings.2,  table.fi.1 = res.proc$table.fi.1, partial.fi.2 = res.proc$table.partial.fi.2) 

      
class(doact.statis.overview) <- c("doact.statis.overview", "list")
class(doact.statis.innerproduct) <- c("doact.statis.innerproduct", "list")
class(doact.statis.compromise) <- c("doact.statis.compromise","list")
class(doact.statis.table) <- c("doact.statis.table","list")

res <- list(Overview = doact.statis.overview, InnerProduct = doact.statis.innerproduct, Compromise = doact.statis.compromise, Table = doact.statis.table)
 
class(res) <- c("mpDOACT.STATIS","list")

print('Processing Complete')


mpPlotInfo <- mpGraphs(res = res, table = res$Overview$column.design.1, DESIGN = DESIGN, main = main, graphs = graphs)
return(mpOutputHandler(res=res, mpPlotInfo=mpPlotInfo))
}
