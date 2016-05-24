#' A harvested classification tree 
#' 
#' The main function of the package, aiming at develop the harvest classification tree. Training data input and
#' @param training original data where 'y' stores classmembership 0 and 1,in the first column, with explanatory variable stores in the second to the last column.
#' @param num.var number of explanatory variables
#' @param numeric.info the vector stores the number of which variable is continuous
#' @param sig significance level (default 0.95)
#' @details The function will return the harvested tree model. Missing values are allowed, and they will be treated accordingly. To use the trained tree model to predict, you can use \code{predict} function in this package.
#' @return An object of class "harvest", which is the result of  algorithm with the following elements for each nodes(nodes are ordered in sequence of harvesting):
#' @return \code{rule} constriants of the node
#' @return \code{total} total number of data points in the node
#' @return \code{'1'} the number of data points belonging to class 1 in the node
#' @return \code{'logchange'} the improvement of log likelihood of deleting the redundent rules by the algorithm for the node
#' @examples
#' data(training)
#' harvest(training,4,3)
#' @export
#' @import rpart stats

harvest=function(training,num.var,numeric.info,sig=0.95){
  names(training)=c('y',paste("x", 1:num.var, sep = ""))
  if(!is.na(numeric.info[1])){
    training.m=ranking(training,numeric.info)$data
  }
  else{
    training.m=training
  }
  formula.rank='y~'
  for(i in 1:num.var){
    if(i!=1){
      formula.rank=paste(formula.rank,paste("x",i,sep=""),sep="+")
    }
    else{
      formula.rank=paste(formula.rank,paste("x",i,sep=""),sep="")
    }
  }
  
  rank.rp <- rpart(as.formula(formula.rank), data=training.m, method="class", parms=list(split="information"), cp=0, minsplit=10, minbucket=5, maxsurrogate=0)
  
  
  
  harnode <- harfunc(rank.rp, data=training.m, varname=paste("x", 1:num.var, sep = ""),sig=0.95)
  
  ori.rule=harnode$original_set
  har.rule=harnode$har.rule
  for(i in 1:length(har.rule)){
    for(j in 1:dim(har.rule[[i]]$rule)[1]){
      for(p in 1:2){
        
        if(!is.na(har.rule[[i]]$rule[j,p])&abs(har.rule[[i]]$rule[j,p])!=Inf){
          variable = dimnames(har.rule[[i]]$rule)[[1]][[j]]
          number = as.numeric(substr(variable,2,nchar(variable)))
          lower=max(training.m[training.m[,number+1]<=har.rule[[i]]$rule[j,p],number+1])
          upper=min(training.m[training.m[,number+1]>=har.rule[[i]]$rule[j,p],number+1])
          harnode$har.rule[[i]]$rule[j,p] = 0.5*(unique(training[training.m[,number+1]==lower,number+1])+unique(training[training.m[,number+1]==upper,number+1]))
        }
      }
    }
  }
  
  for(i in 1:length(ori.rule)){
    for(j in 1:dim(ori.rule[[i]]$bounds)[1]){
      for(p in 1:2){
        
        if(!is.na(ori.rule[[i]]$bounds[j,p])&abs(ori.rule[[i]]$bounds[j,p])!=Inf){
          lower=max(training.m[training.m[,j+1]<=ori.rule[[i]]$bounds[j,p],j+1])
          upper=min(training.m[training.m[,j+1]>=ori.rule[[i]]$bounds[j,p],j+1])
          harnode$original_set[[i]]$bounds[j,p] = 0.5*(unique(training[training.m[,j+1]==lower,j+1])+unique(training[training.m[,j+1]==upper,j+1]))
        }
      }
    }
  }
  
  
  
  return(harnode$har.rule)
}