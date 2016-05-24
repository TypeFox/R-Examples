#' @include coclusterStrategy.R
#' @include optionclasses.R
#' @include Rcoclust.R
#' 
NULL

#' Co-Clustering function for categorical data-sets.
#' 
#' This function performs Co-Clustering (simultaneous clustering of rows and columns )
#' Categorical data-sets using latent block models. It can also be used to perform
#' semi-supervised co-clustering.  
#' 
#' @param data Input data as matrix (or list containing data matrix.)
#' @param semisupervised Boolean value specifying whether to perform semi-supervised co-clustering or not. Make sure to provide row and/or
#' column labels if specified value is true. The default value is false.
#' @param rowlabels Vector specifying the class of rows. The class number starts from zero.
#' Provide -1 for unknown row class. 
#' @param collabels Vector specifying the class of columns.
#' The class number starts from zero. Provide -1 for unknown column class.
#' @param model This is the name of model. The following models exists for various types of data:
#' \tabular{rlll}{
#'     pik_rhol_multi \tab categorical \tab unequal \tab unequal \cr
#'     pi_rho_multi \tab categorical \tab equal \tab unequal \cr
#' }
#' 
#' @param nbcocluster Integer vector specifying the number of row and column clusters respectively.
#' @param strategy Object of class \code{\linkS4class{strategy}}.
#' @param a First hyper-parameter in case of Bayesian settings. Default is 1 (no prior).
#' @param b Second hyper-parameter in case of Bayesian settings. Default is 1 (no prior).
#' 
#' @return Return an object of \code{\linkS4class{BinaryOptions}} or \code{\linkS4class{ContingencyOptions}}
#' or \code{\linkS4class{ContinuousOptions}} depending on whether the data-type is Binary, Contingency or Continuous
#' respectively.
#' 
#' @export
#' 
#' @exportPattern "^[[:alpha:]]+"
#' @useDynLib RCocluster
#' 
#' @examples
#' 
#' ## Simple example with simulated categorical data
#' ## load data
#' data(categoricaldata)
#' ## usage of coclusterCategorical function in its most simplest form
#' out<-coclusterCategorical(categoricaldata,nbcocluster=c(3,2))
#' ## Summarize the output results
#' summary(out)
#' ## Plot the original and Co-clustered data 
#' plot(out)
#' 
coclusterCategorical<-function(data, semisupervised = FALSE, rowlabels = numeric(0), collabels = numeric(0),
                               model = character(0), nbcocluster, strategy = coclusterStrategy(), a=1, b=1)
{
	#Check for data
	if(missing(data)){ stop("Data is missing.")}

		if(!is.list(data))
    {
			if(!is.matrix(data))
      { stop("Data should be matrix.")}
		}
    else
		{
			if(!is.matrix(data[[1]]))
				stop("Data should be matrix.")
			if(!is.numeric(data[[2]])||!is.numeric(data[[3]]))
				stop("Row/Column effects should be numeric vectors.")
			if(length(data[[2]])!=dim(data[[1]])[1]||length(data[[3]])!=dim(data[[1]])[2])
				stop("Dimension mismatch in Row/column effects  and Data.")
		}
  
  #check for row and column labels
  if (semisupervised)
  {
    if(missing(rowlabels)&&missing(collabels))
      stop("Missing row and column labels. At-least one should be provided to perform semi-supervised Co-clustering.")
    if(!missing(rowlabels)&&!is.numeric(rowlabels))
      stop("Row labels should be a numeric vector.")
    if(!missing(collabels)&&!is.numeric(collabels))
      stop("Column labels should be a numeric vector.")
    
    if(!is.list(data)) dimData = dim(data)
    else               dimData = dim(data[[1]])
    
    if(missing(rowlabels))      rowlabels = rep(-1,dimData[1])
    else if(missing(collabels)) collabels = rep(-1,dimData[2])
    
    if(dimData[1]!=length(rowlabels))
      stop("rowlabels length does not match number of rows in data (also ensure to put -1 in unknown labels)")
    
    if(dimData[2]!=length(collabels))
      stop("collabels length does not match number of columns in data (also  ensure to put -1 in unknown labels)")
  }
	#check for number of coclusters
	if(missing(nbcocluster))
	{ stop("Mention number of CoClusters.")}
  
 	if(!is.list(data)) dimData = dim(data)
	else               dimData = dim(data[[1]])
		
	if(dimData[1]<nbcocluster[1]) stop("Number of Row cluters exceeds numbers of rows.")
	if(dimData[2]<nbcocluster[2])	stop("Number of Column cluters exceeds numbers of columns.")
	#check for Algorithm name (and make it compatible with version 1)
	if(strategy@algo=="XEMStrategy")
  {
    warning("The algorithm 'XEMStrategy' is renamed as BEM!")
    strategy@algo == "BEM"
  }
  else if(strategy@algo == "XCEMStrategy")
  {
    warning("The algorithm 'XCEMStrategy' is renamed as BCEM!")
    strategy@algo = "BCEM"
  }
  else if(strategy@algo!="BEM" && strategy@algo!="BCEM" && strategy@algo!="BSEM" )
  { stop("Incorrect Algorithm, Valide algorithms are: BEM, BCEM, BSEM") }
	#check for stopping criteria
	if(strategy@stopcriteria!="Parameter" && strategy@stopcriteria!="Likelihood")
		stop("Incorrect stopping criteria, Valid stopping criterians are: Parameter, Likelihood")
  # check model
  if(missing(model)){ model = "pik_rhol_multi"}
  else 
  {
    if(model!="pi_rho_multi" && model!="pik_rhol_multi")
    {
      stop("Incorrect Model, Valid categorical models are:pik_rhol_multi, pi_rho_multi")
		}
  }
  #check init method
  if(length(strategy@initmethod)==0){ strategy@initmethod = "randomInit"}
  # checking for hyperparameters
  if ((a <=0)||(b<=0)) { stop("hyper-parameters must be positive")}
  strategy@hyperparam = c(a,b);
  # create object
  inpobj<-new("CategoricalOptions",data = data
             , rowlabels = rowlabels, collabels = collabels
             , semisupervised = semisupervised
             , datatype = "categorical"
             , model = model
             , nbcocluster = nbcocluster
             , strategy = strategy)

  .Call("CoClustmain",inpobj,PACKAGE = "blockcluster")
  cat(inpobj@message,"\n")
  return(inpobj)
}

