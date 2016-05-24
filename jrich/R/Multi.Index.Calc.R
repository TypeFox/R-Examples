#'
#' @title Jack-knife indices in n topologies one time.
#'
#' @description The function calculates the indices values for a MultiData list 
#' one time. 
#'
#' @param MultiData is the list of Trees and distributions to evaluate, a list object.
#'
#' @param  jtip is the proportion of terminals to delete, real (range 0-1).
#'
#' @param jtopol is the proportion of topologies to delete, real (range 0-1).
#' 
#' @return Returns the indices values.
#'

#' @examples
#'  ## get the library
#'  library(jrich)
#'  
#'  ## load the data
#'  data(Multitaxon1) 
#' 
#'  Multi.Index.Calc(Multitaxon1, jtip = 0, jtopol = 0)
#'
#'@author Miranda-Esquivel Daniel R.
#'
#'


Multi.Index.Calc <-
function (MultiData = MultiData, jtip = 0, jtopol = 0) {
  
  for (i in 1:length(MultiData)){    
    if (jtopol > runif(1)){
      temp.Index.Value <- Calculate.Index(tree = MultiData[[i]][[1]],
                                     distribution = MultiData[[i]][[2]], jtip)
    }else{
      temp.Index.Value <- Calculate.Index(tree = MultiData[[i]][[1]],
                                     distribution = MultiData[[i]][[2]], 0)
      temp.Index.Value$jtip <- abs(jtip)
    }
    
    if (i == 1){
      def.Index.Value <- temp.Index.Value 
    }else{
      def.Index.Value <-  Sum.Indices.2.Topologies(temp.Index.Value, def.Index.Value)
    }
  }
  
  def.Index.Value$jtopol <-  abs(jtopol)
  
  return(def.Index.Value)
}
