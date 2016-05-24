## Canonical STATIS
mpCANOSTATIS <- function(data, column.design, row.design, normalization = 'MFA', row.preprocess = 'None', column.preprocess = 'Center_1Norm', table.preprocess = 'Sum_PCA', make.columndesign.nominal = TRUE, make.rowdesign.nominal = TRUE, DESIGN = NULL, make.design.nominal = TRUE, graphs = TRUE)
{	
	main <- deparse(substitute(data))
	DESIGN <- designCheck(data, DESIGN, make.design.nominal)

	res.preproc <- mpSTATIS.preprocess(data, column.design = column.design, row.design = row.design , row.preprocess=row.preprocess, column.preprocess=column.preprocess, table.preprocess=table.preprocess, make.columndesign.nominal=make.columndesign.nominal, make.rowdesign.nominal=make.rowdesign.nominal)

	res.proc <- mpCANOSTATIS.core(res.preproc$data, num.obs = res.preproc$num.obs, res.preproc$column.design, res.preproc$row.design, num.groups = res.preproc$num.groups, normalization = normalization, masses = NULL)
	
	print('Processing Completed')
	#########################
	# Results
	#########################

	canostatis.overview <- list(data = res.preproc$data, groupmatrix = res.preproc$column.design, row.design = res.preproc$row.design,
						preprocess.data = res.preproc$data.preprocessed, num.groups = res.preproc$num.groups, num.obs = res.preproc$num.obs, 
						row.preprocess = res.preproc$row.preprocess, column.preprocess = res.preproc$column.preprocess, table.preprocess = res.preproc$table.preprocess)

	canostatis.innerproduct <- list(mahalanobis = res.proc$mahalanobis, S=res.proc$scalarProductMatrices, C = res.proc$C, RVMatrix = res.proc$rvMatrix, 
						ci = res.proc$ci, cj = res.proc$cj, eigs.vector = res.proc$eigs.vector, eigs = res.proc$eigs, fi = res.proc$fi, t = res.proc$tau, 
						alphaWeights = res.proc$alphaWeights)

	canostatis.compromise <- list(compromise = res.proc$compromise, compromise.ci = res.proc$compromise.ci, compromise.cj = res.proc$compromise.cj, 
						compromise.eigs.vector = res.proc$compromise.eigs.vector, compromise.eigs = res.proc$eigs, compromise.fi = res.proc$compromise.fi,
						compromise.t = res.proc$compromise.tau)

	canostatis.table <- list(m = res.proc$masses, eigs = res.proc$table.eigs, eigs.vector = res.proc$table.eigs.vector, 
            			Q = res.proc$table.Q, fi = res.proc$table.fi, partial.fi = res.proc$table.partial.fi,
         		    	cj = res.proc$table.cj, ci = res.proc$table.ci, t =res.proc$table.tau, partial.fi.array = res.proc$table.partial.fi.array)  

	class(canostatis.overview) <- c("canostatis.overview","list")
	class(canostatis.innerproduct) <- c("canostatis.innerproduct","list")
	class(canostatis.compromise) <- c("canostatis.compromise","list")
	class(canostatis.table) <- c("canostatis.table","list")

	res <- list(Overview = canostatis.overview, InnerProduct = canostatis.innerproduct, Compromise = canostatis.compromise, Table = canostatis.table)

	class(res) <- c("mpSTATIS", "list")
	
	mpPlotInfo <- mpGraphs(res = res, table = res$Overview$groupmatrix, DESIGN = DESIGN, main = main, graphs = graphs)
	
	return(mpOutputHandler(res = res, mpPlotInfo = mpPlotInfo))	
}


