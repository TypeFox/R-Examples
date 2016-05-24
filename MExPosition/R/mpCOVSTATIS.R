mpCOVSTATIS <- function(data, normalization = 'None', masses= NULL, table = NULL, make.table.nominal = TRUE,DESIGN = NULL, make.design.nominal = TRUE, graphs=TRUE)
{ 
#######################
# Data preparation
#######################
  main <- deparse(substitute(data))		
  
  if(is.null(data))
  {   stop('You have not provided any data')
  }
	
  if(!is.matrix(data))
  {   data <- as.matrix(data)
  }
	
  if(sum(is.na(data)>0))
  { stop('Missing data not allowed')
  }
  
  DESIGN <- designCheck(data, DESIGN, make.design.nominal)
	
#########################
# Running DISTATIS
#########################

# Core processing of COVSTATIS
   res.proc <- mpCOVSTATIS.core(data, normalization = normalization, masses = masses, table = table, make.table.nominal= TRUE)


  print('Processing Complete.')
 		
#########################
# Results
#########################

covstatis.overview <- list(data = res.proc$data, normalization = res.proc$normalization, table = res.proc$table,
					           num.groups=dim(res.proc$table.partial.fi.array)[3])

covstatis.innerproduct <- list(S = res.proc$S, C=res.proc$C, rvMatrix = res.proc$rvMatrix, ci = res.proc$ci, cj = res.proc$cj, eigs.vector= res.proc$eigs.vector, 
                     eigs = res.proc$eigs, fi = res.proc$fi, t = res.proc$tau, alphaWeights = res.proc$alphaWeights)

covstatis.compromise <- list(compromise = res.proc$compromise, compromise.ci = res.proc$compromise.ci, compromise.cj=res.proc$compromise.cj, 
                      compromise.eigs.vector = res.proc$compromise.eigs.vector, compromise.eigs = res.proc$compromise.eigs, compromise.fi = res.proc$compromise.fi,
                      compromise.t = res.proc$compromise.tau)

covstatis.table <- list(m = res.proc$masses, eigs = res.proc$table.eigs, eigs.vector = res.proc$table.eigs.vector, 
                      fi = res.proc$table.fi, partial.fi = res.proc$table.partial.fi, ci = res.proc$table.ci, t =res.proc$table.tau, 
                      partial.fi.array = res.proc$table.partial.fi.array, Q = res.proc$table.Q, cj = res.proc$table.cj)   					

class(covstatis.overview) <- c("covstatis.overview","list")
class(covstatis.innerproduct) <- c("covstatis.innerproduct","list")
class(covstatis.compromise) <- c("covstatis.compromise","list")
class(covstatis.table) <- c("covstatis.table","list")
  
res<- list(Overview = covstatis.overview, InnerProduct = covstatis.innerproduct, Compromise = covstatis.compromise, Table=covstatis.table)  

class(res) <- c("mpCOVSTATIS","list")

mpPlotInfo <- mpGraphs(res, table = res$Overview$table, DESIGN = DESIGN, main = main, graphs = graphs)
  
return(mpOutputHandler(res=res, mpPlotInfo = mpPlotInfo))
  
}
