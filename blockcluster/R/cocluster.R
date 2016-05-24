#' @include coclusterStrategy.R
#' @include optionclasses.R
#' @include Rcoclust.R
#' 
NULL

#' Co-Clustering function.
#' 
#' This function performs Co-Clustering (simultaneous clustering of rows and columns ) for Binary, Contingency
#' and Continuous data-sets using latent block models.It can also be used to perform semi-supervised co-clustering.  
#' 
#' @param data Input data as matrix (or list containing data matrix, numeric vector for row effects and numeric 
#'        vector column effects in case of contingency data with known row and column effects.)
#' @param datatype This is the type of data which can be "binary" , "contingency", "continuous" or "categorical".
#' @param semisupervised Boolean value specifying whether to perform semi-supervised co-clustering or not. Make sure to provide row and/or
#' column labels if specified value is true. The default value is false.
#' @param rowlabels Vector specifying the class of rows. The class number starts from zero. Provide -1 for unknown row class. 
#' @param collabels Vector specifying the class of columns. The class number starts from zero. Provide -1 for unknown column class.
#' @param model This is the name of model. The following models exists for various types of data:
#' \tabular{rlll}{
#'     Model  \tab Data-type \tab Proportions \tab Dispersion/Variance \cr
#'     pik_rhol_epsilonkl(Default) \tab binary \tab unequal \tab unequal \cr
#'     pik_rhol_epsilon \tab binary \tab unequal \tab equal \cr
#'     pi_rho_epsilonkl \tab binary \tab equal \tab unequal \cr
#'     pi_rho_epsilon \tab binary \tab equal \tab equal \cr
#'     pik_rhol_sigma2kl(Default) \tab continuous \tab unequal \tab unequal \cr
#'     pik_rhol_sigma \tab continuous \tab unequal \tab equal \cr
#'     pi_rho_sigma2kl \tab continuous \tab equal \tab unequal \cr
#'     pi_rho_sigma2 \tab continuous \tab equal \tab equal \cr
#'     pik_rhol_unknown(default) \tab contingency \tab unequal \tab N.A \cr
#'     pi_rho_unknown \tab contingency \tab equal \tab N.A \cr
#'     pik_rhol_known \tab contingency \tab unequal \tab N.A \cr
#'     pi_rho_known \tab contingency \tab equal \tab N.A \cr
#'     pik_rhol_multi \tab categorical \tab unequal \tab unequal \cr
#'     pi_rho_multi \tab categorical \tab equal \tab unequal \cr
#' }
#' 
#' @param nbcocluster Integer vector specifying the number of row and column clusters respectively.
#' @param strategy Object of class \code{\linkS4class{strategy}}.
#' @return Return an object of \code{\linkS4class{BinaryOptions}} or \code{\linkS4class{ContingencyOptions}}
#' or \code{\linkS4class{ContinuousOptions}} depending on whether the data-type is Binary, Contingency or Continuous
#' respectively.
#' 
#' @export
#' 
#' @exportPattern "^[[:alpha:]]+"
#' @useDynLib blockcluster
#' 
#' @examples
#' 
#' # Simple example with simulated binary data
#' #load data
#' data(binarydata)
#' #usage of cocluster function in its most simplest form
#' out<-cocluster(binarydata,datatype="binary",nbcocluster=c(2,3))
#' #Summarize the output results
#' summary(out)
#' #Plot the original and Co-clustered data 
#' plot(out)
#' 
#' 


cocluster<-function( data, datatype
                   , semisupervised = FALSE, rowlabels = numeric(0), collabels = numeric(0)
                   , model = character(0), nbcocluster, strategy = coclusterStrategy()) 
{
  #Check for data
  if(missing(data)){
    stop("Data is missing.")
  }
  else
  {
    if(!is.list(data)){
      if(!is.matrix(data))
        stop("Data should be matrix.")
    }else
    {
      if(!is.matrix(data[[1]]))
        stop("Data should be matrix.")
      if(!is.numeric(data[[2]])||!is.numeric(data[[3]]))
        stop("Row/Column effects should be numeric vectors.")
      if(length(data[[2]])!=dim(data[[1]])[1]||length(data[[3]])!=dim(data[[1]])[2])
        stop("Dimension mismatch in Row/column effects  and Data.")
    }
  }
  
  #check for row and column labels
  if(semisupervised)
  {
    if(missing(rowlabels)&&missing(collabels))
      stop("Missing row and column labels. At-least one should be provided to perform semi-supervised Co-clustering.")
    if(!missing(rowlabels)&&!is.numeric(rowlabels))
      stop("Row labels should be a numeric vector.")
    if(!missing(collabels)&&!is.numeric(collabels))
      stop("Column labels should be a numeric vector.")
    
    if(!is.list(data))
      dimention = dim(data)
    else
      dimention = dim(data[[1]])
    
    if(missing(rowlabels))
      rowlabels = rep(-1,dimention[1])
    else if(missing(collabels))
      collabels = rep(-1,dimention[2])
    
    
    if(dimention[1]!=length(rowlabels))
      stop("rowlabels length does not match number of rows in data (also ensure to put -1 in unknown labels)")
    
    if(dimention[2]!=length(collabels))
      stop("collabels length does not match number of columns in data (also  ensure to put -1 in unknown labels)")
  }
  
  #check for number of coclusters
  if(missing(nbcocluster))
  {
    stop("Mention number of CoClusters.")
  }
  else
  {
    if(!is.list(data))
      dimention = dim(data)
    else
      dimention = dim(data[[1]])
    
    if(dimention[1]<nbcocluster[1])
      stop("Number of Row cluters exceeds numbers of rows.")
    
    if(dimention[2]<nbcocluster[2])
      stop("Number of Column cluters exceeds numbers of columns.")
  }
  
  
  #check for Algorithm name (and make it compatible with version 1)
  if(strategy@algo=="XEMStrategy"){
    warning("The algorithm 'XEMStrategy' is renamed as BEM!")
    strategy@algo == "BEM"
  }else if(strategy@algo == "XCEMStrategy"){
    warning("The algorithm 'XCEMStrategy' is renamed as BCEM!")
    strategy@algo = "BCEM"
  }else if(strategy@algo!="BEM" && strategy@algo!="BCEM" && strategy@algo!="BSEM" )
    stop("Incorrect Algorithm, Valide algorithms are: BEM, BCEM, BSEM") 
  
  #check for stopping criteria
  
  if(strategy@stopcriteria!="Parameter" && strategy@stopcriteria!="Likelihood")
    stop("Incorrect stopping criteria, Valid stopping criterians are: Parameter, Likelihood")
  
  #check for datatype and models and create input object to be passed in .Call function.
  if (missing(datatype)) {
    stop("Mention datatype.")
  } 
  else{
    if(datatype == "binary"){
      if(missing(model)){
        model = "pik_rhol_epsilonkl"
      }
      else 
      {
        if(model!="pik_rhol_epsilonkl" && model!="pik_rhol_epsilon" && 
            model!="pi_rho_epsilonkl" && model!="pi_rho_epsilon")
        {
          stop("Incorrect Model, Valid Binary models are:pik_rhol_epsilonkl, pik_rhol_epsilon
                  pi_rho_epsilonkl, pi_rho_epsilon")
        }
      }
      
      if(length(strategy@initmethod)==0){
        strategy@initmethod = "CEMInit"
      }
      
      inpobj<-new("BinaryOptions",data = data, rowlabels = rowlabels, collabels = collabels, semisupervised = semisupervised,
          datatype = datatype, model = model,nbcocluster = nbcocluster, strategy = strategy)
    }
    
    else if(datatype == "continuous"){
      if(missing(model)){
        model = "pik_rhol_sigma2kl"
      }
      else if(model!="pik_rhol_sigma2kl" && model!="pik_rhol_sigma2" && 
          model!="pi_rho_sigma2kl" && model!="pi_rho_sigma2"){
        stop("Incorrect Model, Valid Continuous models are: pik_rhol_sigma2kl, pik_rhol_sigma2, pi_rho_sigma2kl, pi_rho_sigma2")
      }
      
      if(length(strategy@initmethod)==0){
        strategy@initmethod = "CEMInit"
      }
      inpobj<-new("ContinuousOptions",data = data,rowlabels = rowlabels, collabels = collabels, semisupervised = semisupervised, 
          datatype = datatype, model = model,nbcocluster = nbcocluster,strategy = strategy)
    }
    
    else if(datatype == "contingency"){
      if(missing(model)&& !is.list(data)){
        model = "pik_rhol_unknown"
      }
      else if(missing(model) && is.list(data))
      {
        model = "pik_rhol_known"
      }
      else if(model!="pik_rhol_unknown" && model!="pik_rhol_known" && 
          model!="pi_rho_unknown" && model!="pi_rho_known"){
        stop("Incorrect Model, Valid Contingency models are:pik_rhol_unknown, pik_rhol_known, pi_rho_unknown, pi_rho_known")
      }
      
      if((model=="pi_rho_known"||model=="pik_rhol_known")&& (length(data)!=3))
      {
        stop("Missing Row/Column effects.") 
      }
      
      if(length(strategy@initmethod)==0){
        if((model=="pi_rho_known"||model=="pik_rhol_known"))
        {
          strategy@initmethod = "RandomInit"
        }
        else{
          strategy@initmethod = "CEMInit"
        }
      }
      else
      {
        if(strategy@initmethod!="RandomInit"&&(model=="pi_rho_known"||model=="pik_rhol_known"))
        {stop("Incorrect initialization method, valid method(s) are: RandomInit")}
        else if(strategy@initmethod!="CEMInit"&&(model=="pi_rho_unknown"||model=="pik_rhol_unknown"))
          stop("Incorrect initialization method, valid method(s) are: CEMInit")
      }
      if(!is.list(data))
        inpobj<-new("ContingencyOptions",data = data,semisupervised = semisupervised,rowlabels = rowlabels, collabels = collabels,
            datatype = datatype, model = model, nbcocluster = nbcocluster, strategy = strategy)
      else
        inpobj<-new("ContingencyOptions",data = data[[1]], semisupervised = semisupervised,rowlabels = rowlabels, collabels = collabels,
            datatype = datatype, model = model, nbcocluster = nbcocluster, strategy = strategy,datamui=data[[2]],datanuj=data[[3]])
      
    }
    else if(datatype == "categorical"){
      if(missing(model)){
        model = "pik_rhol_multi"
      }
      else 
      {
        if(model!="pi_rho_multi" && model!="pik_rhol_multi")
        {
          stop("Incorrect Model, Valid categorical models are:pik_rhol_multi, pi_rho_multi")
        }
      }
      
      if(length(strategy@initmethod)==0){
        strategy@initmethod = "RandomInit"
      }
      
      inpobj<-new("CategoricalOptions",data = data, rowlabels = rowlabels, collabels = collabels, semisupervised = semisupervised,
          datatype = datatype, model = model,nbcocluster = nbcocluster, strategy = strategy)
    }
    else
    {stop("Invalid datatype, Valid types are: binary , contingency, continuous and categorical")}
  }
  
  
  .Call("CoClustmain",inpobj,PACKAGE = "blockcluster")
  
  
  cat(inpobj@message,"\n")
  
  return(inpobj)
}


