mpDISTATIS <- function(data, sorting = 'No', normalization = 'None', masses= NULL, table = NULL, make.table.nominal = TRUE,DESIGN = NULL, make.design.nominal = TRUE, graphs=TRUE)
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

# Core processing of DISTATIS
   res.proc <- mpDISTATIS.core(data, sorting = sorting, normalization = normalization, masses = masses, table = table, make.table.nominal= TRUE)
 		
#########################
# Results
#########################

distatis.overview <- list(data = res.proc$data, sorting = res.proc$sorting, normalization = res.proc$normalization, table = res.proc$table,
					           num.groups=dim(res.proc$table.partial.fi.array)[3])

distatis.innerproduct <- list(S = res.proc$S, norm.S = res.proc$norm.S, C=res.proc$C, ci = res.proc$ci, cj = res.proc$cj, eigs.vector= res.proc$eigs.vector, 
                     eigs = res.proc$eigs, fi = res.proc$fi, t = res.proc$tau, alphaWeights = res.proc$alphaWeights) 

distatis.compromise <- list(compromise = res.proc$compromise, compromise.ci = res.proc$compromise.ci, compromise.cj=res.proc$compromise.cj, 
                      compromise.eigs.vector = res.proc$compromise.eigs.vector, compromise.eigs = res.proc$compromise.eigs, compromise.fi = res.proc$compromise.fi,
                      compromise.t = res.proc$compromise.tau)

distatis.table <- list(m = res.proc$masses, eigs = res.proc$table.eigs, eigs.vector = res.proc$table.eigs.vector, 
            		      fi = res.proc$table.fi, partial.fi = res.proc$table.partial.fi, ci = res.proc$table.ci, t =res.proc$table.tau, 
                      partial.fi.array = res.proc$table.partial.fi.array, Q = res.proc$table.Q, cj = res.proc$table.cj)  		

class(distatis.overview) <- c("distatis.overview","list")
class(distatis.innerproduct) <- c("distatis.innerproduct","list")
class(distatis.compromise) <- c("distatis.compromise","list")
class(distatis.table) <- c("distatis.table","list")

print('Processing Complete.')
	
res<- list(Overview = distatis.overview, InnerProduct = distatis.innerproduct, Compromise = distatis.compromise, Table=distatis.table)	

class(res) <- c("mpDISTATIS","list")

mpPlotInfo <- mpGraphs(res, table = res$Overview$table, DESIGN = DESIGN, main = main, graphs = graphs)
  
return(mpOutputHandler(res=res, mpPlotInfo = mpPlotInfo))
  
}
