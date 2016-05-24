#' Learning Hybric Bayesian Networks
#' 
#' Get a directed acyclic graph using the method \bold{hill climbing}.
#' It is a mathematical optimization technique which belongs to the 
#' family of local search. It is an iterative algorithm which deals with
#' discrete and continuous variables. 
#' 
#' @param dataset A dataset with discrete and continuous variables. Discrete variables 
#' must be of class \code{"factor"}, if not they are transformed into factores. 
#' @param numIntervals A \code{"numeric"} value containing the number of intervals 
#' to create a discrete dataset. The method used to split the values is by equal width.
#' By default it is \code{NULL} and means the variables are not discretized to get the network.
#' @return The output is a \code{"bn"} object containing the learned graph.
#' @seealso \link{hc}
#' @importFrom bnlearn hc
#' @export
#' @examples
#' 
#' ## Data
#' data(ecoli)
#' ecoli <- ecoli[,-1] ## Sequence Name
#' 
#' ## DAG1
#' dag1 <- LearningHC(ecoli)
#' dag1
#' plot(dag1)
#' 
#' ## DAG2
#' dag2 <- LearningHC(ecoli, numIntervals = 10)
#' dag2
#' plot(dag2)
#' 
#' 
LearningHC <- function(dataset, numIntervals=NULL)
{
  if(is.numeric(numIntervals)){
    ## Get the discrete dataset
    dataset <- discretizeVariablesEWdis(dataset, numIntervals, factor=TRUE)
    
    ## Estimate DAG
    dag <- hc(dataset, score="loglik")
  } else{
    
    ## Get discrete variables as factor
    pos <- which(!sapply(dataset, is.numeric))
    if(length(pos)!=0) for(i in pos) dataset[,i] <- as.factor(dataset[,i])
    
    ## Estimate DAG
    dag <- hc(dataset)
  }
  return(dag)
}

#' Get Relationships in a Network
#' 
#' Extract the relationship between the variables of a dataset using
#' the obtained network
#' 
#' @param graph A structural network of the class \code{"graphNEL"},
#' \code{"network"} or \code{"bn"}.
#' @param nameVars A character array giving the names of the variables in the graph. By default it's NULL,
#' only put it when a graph of class \code{"network"} is used.
#' @return A list of elements. Each element contains a vector with the name of a child 
#' and their parents.
#' @export
#' @examples
#' 
#' ## Data
#' data(ecoli)
#' ecoli <- ecoli[,-1] ## Sequence Name
#' 
#' ## DAG1
#' dag1 <- LearningHC(ecoli)
#' dag1
#' plot(dag1)
#' getChildParentsFromGraph(dag1)
#' 
#' ## DAG2
#' dag2 <- LearningHC(ecoli, numIntervals = 10)
#' dag2
#' plot(dag2)
#' getChildParentsFromGraph(dag2)
#' 
getChildParentsFromGraph <- function(graph, nameVars=NULL)
{
  if(class(graph)=="graphNEL") type <- 1
  if(class(graph)=="network")  type <- 2
  if(class(graph)=="bn")       type <- 3
  switch(type, 
         
## type <- 1 
{namenodes <- graph@nodes;
 edgenodes <- graph@edgeL; 
 childrenAndParents <- list();
 for(i in 1:length(edgenodes)) childrenAndParents[[length(childrenAndParents)+1]] <- 
   c(namenodes[i],namenodes[edgenodes[[i]]$edges]);
},

## type <- 2
{childrenAndParents = list();
 for(j in 1:length(graph$nodes)){
   Child <- nameVars[graph$nodes[[j]]$idx]
   Parents <- nameVars[graph$nodes[[j]]$parents]
   childrenAndParents[[length(childrenAndParents)+1]] <- c(Child, Parents)
 }
},

## type <- 3
{if(is.null(nameVars)) nameVars <- names(graph$nodes);
 childrenAndParents = list();
 for(i in 1:length(graph$nodes)){
   Child <- nameVars[i]
   Parents <- graph$nodes[[i]]$parents
   childrenAndParents[[length(childrenAndParents)+1]] <- c(Child, Parents)
 }
}
 )
return(childrenAndParents)
}

