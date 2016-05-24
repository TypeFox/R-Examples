#'Checks for preconditions
#' Checks to see if input satisfies preconditions 
#'checks if inputs are in an expected format
#'and stops right away if something is wrong
#'@noRd
#'@param formula formula for models
#'@param data input dataframe
#'@keywords internal
precondition_checks <-function(formula, data)
{
  if( nrow(data)==0)
    stop("data contains 0 rows after removing NA in dependent variable", call.=F)
  
  depvar_name <- as.character(formula)[2]
  if(nlevels(data[,depvar_name]) !=2)
    stop("Dependent variables must have two levels", call.=F) 
}


#'fast implementation of adaboost using Rcpp
#'
#'main internal function where both M1 and real 
#'adaboost calls are redirected
#'
#'@noRd
#'@import rpart
#'@param formula Formula for models
#'@param data Input dataframe
#'@param nIter no. of classifiers 
#'@param method a string, either "M1" or "SAMME.R"
#'@return adaboost_object An object of class adaboost/real_adaboost containing each tree, and its corresponding weight
#'@keywords internal
adaboost_fast<-function(formula, data, nIter, method)
{
  if(!all(method %in% c("M1","SAMME.R") ))
    stop(paste("method must be M1 or SAMME.R. It is",method), call.=F)
  #initial cleaning
  #(1) remove NAs in the dependent variable
  #(2) make dependent variable a factor if not already
  depvar_name <- as.character(formula)[2]
  data <-data[ !(is.na(data[,depvar_name])) , ]
  
  if( class(data[,depvar_name]) !="factor" )
    data[,depvar_name] <- factor(data[,depvar_name])
    
  #check if preconditions are satisfied
  precondition_checks(formula, data)
  
  #convert the factor to 0/1 and remember names
  classnames_map <- levels(data[,depvar_name])
  
  names(classnames_map) <- c("A","B")
  vardep = ifelse(data[,depvar_name]==classnames_map["A"],0,1)  
  
  #print(depvar_name)
  #print(levels(data[,depvar_name]))
  #print(classnames_map)
  #print(vardep)
  
  rcpp_ada_obj <- adaboost_main_loop_(formula, data, nIter, wrap_rpart,
                                      vardep, classnames_map, method)
  trees_list <- rcpp_ada_obj$trees
  coeff_vector<- rcpp_ada_obj$weights
  
  
  
  ceoff_vector <- rcpp_ada_obj$coeff_vector
  adaboost_object <- list(formula = formula, trees=trees_list, weights = coeff_vector, 
                          classnames = classnames_map, dependent_variable = depvar_name)
  if(method=="M1")
    class(adaboost_object) <- "adaboost"
  else
    class(adaboost_object) <-"real_adaboost"
  return(adaboost_object)
  
}

