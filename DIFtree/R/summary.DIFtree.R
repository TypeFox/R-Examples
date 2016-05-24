#' Summary for fitted Item focussed Trees 
#' 
#' @description
#' The function takes an object of class \code{"DIFtree"} and returns an useful summary 
#' with an overiew of all executed splits during the estimation procedure.
#' 
#' @param object Object of class \code{\link[DIFtree]{DIFtree}}
#' @param x Object of class \code{\link[DIFtree]{summary.DIFtree}}
#' @param ... Further arguments passed to or from other methods 
#' 
#' @return Object of class \code{"summary.DIFtree"}. 
#' An object of class \code{"summary.DIFtree"} is a list containing the following components:
#' 
#' \item{stats}{Useful overview of detected DIF items, responsible variables and executed splits}
#' \item{nosplits}{Total number of executed splits during the estimation procedure}
#' 
#' @author Moritz Berger <moritz.berger@stat.uni-muenchen.de> \cr \url{http://www.statistik.lmu.de/~mberger/}
#' 
#' @references 
#' Berger, Moritz and Tutz, Gerhard (2015): Detection of Uniform and Non-Uniform Differential Item Functioning 
#' by Item Focussed Trees, Cornell University Library, arXiv:1511.07178
#' 
#' Tutz, Gerhard and Berger, Moritz (2015): Item Focused Trees for the Identification of Items
#' in Differential Item Functioning, Psychometrika, published online, DOI: 10.1007/s11336-015-9488-3 
#' 
#' @seealso \code{\link[DIFtree]{DIFtree}}, \code{\link[DIFtree]{plot.DIFtree}}, \code{\link[DIFtree]{predict.DIFtree}}
#' 
#' @examples 
#' data(data_sim)
#'  
#' Y <- data_sim[,1]
#' X <- data_sim[,-1]
#'  
#' \dontrun{
#'  
#' mod <- DIFtree(Y=Y,X=X,model="Logistic",type="udif",alpha=0.05,nperm=1000,trace=TRUE)
#'  
#' summary(mod)
#' }
#' 
#
#' @method summary DIFtree
#' @exportClass summary.DIFtree 
#' @export

summary.DIFtree <-
function(object, # object of class DIFtree
                            ...){
  
  to_return <- list(call=object$call)
  
  model <- which(c("Rasch","Logistic") %in% paste(object$call))
  if(model==2){
    type <- which(c("udif","dif","nudif") %in% paste(object$call))
  } else{
    type <- 1
  }
  nitems <- object$items
  overview <- infos_summary(object$splits,1:nitems,model,type)
  
  if(model==2 & type==2){
    nos <- nrow(rbind(object$splits[[1]],object$splits[[2]]))
  } else{
    nos <- nrow(object$splits)
  }
  to_return$stats  <- overview 
  to_return$nosplits <- nos  

  class(to_return) <- "summary.DIFtree"
  to_return
  
}

#' @rdname summary.DIFtree 
#' @method print summary.DIFtree
#' @export

print.summary.DIFtree <-
  function(x, # object of class summary.DIFtree 
           ...){
    
    model <- which(c("Rasch","Logistic") %in% paste(x$call))
    if(model==2){
      type <- which(c("udif","dif","nudif") %in% paste(x$call))
    } else{
      type <- 1
    }
    cat("\n")
    if(model==1){
      cat("Item focussed Trees based on the Rasch Model:\n")
    } else{
      cat("Item focussed Trees based on the Logistic Regression Approach",c("(uniform","(uniform and non-uniform","(non-uniform")[type],"DIF):\n")
    }
    cat("\n")
    cat("Call:\n",paste(deparse(x$call)),"\n")
    cat("\n")
    cat("----------\n")
    cat("\n")
    cat("Overview:\n")
    cat("\n")
    print(x$stats)
    cat("\n")
    cat("Total number of Splits:", x$nosplits)
  }
