#' Item Focused Trees for the Identification of Items in Differential Item Functioning 
#' 
#' @description
#' A function to estimate item focused trees for simultaneous selection of items and variables 
#' that induce DIF (Differential Item Functioning) based on the Rasch Model or the 
#' Logistic Regression Approach for DIF detection.
#' The basic method of item focussed recursive partitioning in Rasch Models is described in Tutz and Berger (2015).
#' 
#' @param Y Matrix or Data.frame of binary 0/1 response (rows correspond to persons, columns correspond to items)
#' @param X Data.frame of (not scaled) covariates (rows correspond to persons, columns correspond to covariates)
#' @param model Type of model to be fitted; can be \code{"Rasch"} or \code{"Logistic"}.
#' @param type  Type of DIF to be modelled; one out of \code{"udif"}, \code{"dif"} and \code{"nudif"}. 
#' For \code{"Rasch"} model only uniform DIF can be modelled and therefore \code{type} will be ignored.
#' @param alpha Global significance level for the permutation tests
#' @param nperm Number of permutations used for the permutation tests
#' @param trace If true, information about the estimation progress is printed
#' @param penalize If true, a small ridge penalty is added to ensure existence of model parameters; only for \code{"Rasch"} model.
#' @param x Object of class \code{"DIFtree"}
#' @param ... Further arguments passed to or from other methods
#' 
#' @details 
#' The methods require 0/1 coded answers on binary items. 
#' Items with DIF are gradually identified by recursive partitioning.
#' 
#' For \code{"Rasch"} model one yields a model with linear predictors 
#' \deqn{eta_{pi}=theta_p-tr_i(x_p),}
#' where \eqn{theta_p} correspond to the ability and \eqn{x_p} correspond to the covariate vector of person p. 
#' 
#' For \code{"Logistic"} model one yields a model with linear predictors 
#' \itemize{
#' \item Uniform DIF, \code{type="udif"}
#' \deqn{eta_{pi}=S_p beta_i+tr_i(x_p),}
#' where \eqn{S_p} corresponds to the test score and \eqn{x_p} corresponds to the covariate vector of person p.
#' \item DIF and Non-Uniform DIF, \code{type="dif", "nudif"}
#' \deqn{eta_{pi}=tr_i(x_p)+tr_i(S_p,x_p),}
#' where \eqn{S_p} corresponds to the test score and \eqn{x_p} corresponds to the covariate vector of person p. 
#' }
#'  
#' 
#' Significance of each split is verified by permutation tests. The result of the permutation tests 
#' can strongly depend on the number of permutations \code{nperm}.
#' In the case of pure terminal nodes estimates of the model do not exist. If \code{penalize=TRUE} 
#' a small ridge penalty is added during estimation to ensure existence of all parameters. 
#'
#' @return Object of class \code{"DIFtree"}. 
#' An object of class \code{"DIFtree"} is a list containing the following components:
#' 
#' \item{splits}{Matrix with detailed information about all executed splits during the estimation process}
#' \item{coefficients}{List of estimated coefficients for items with and without DIF. 
#' Structure of \code{coefficients} depends on \code{model} and \code{type}.}
#' \item{pvalues}{P-values of each permutation test during the estimation process}
#' \item{devs}{Maximal value statistics \eqn{T_j} of selected variable of each iteration during the estimation process}
#' \item{crit}{Critical values of each permutation test during the estimation process}
#' \item{Y}{Response matrix used in the estimation process}
#' \item{X}{Model matrix used in the estimation process}
#' \item{persons}{Number of persons} 
#' \item{items}{Number of items} 
#' 
#' @author Moritz Berger <moritz.berger@stat.uni-muenchen.de> \cr \url{http://www.statistik.lmu.de/~mberger/}
#' 
#' @references 
#' Berger, Moritz and Tutz, Gerhard (2015): Detection of Uniform and Non-Uniform Differential Item Functioning 
#' by Item Focussed Trees, Cornell University Library, arXiv:1511.07178
#' 
#' Swaminathan, Hariharan and Rogers, H Jane (1990): Detecting differential item functioning 
#' using logistic regression procedures, Journal of Educational measurements 27(4), 361-370
#' 
#' Tutz, Gerhard and Berger, Moritz (2015): Item Focused Trees for the Identification of Items
#' in Differential Item Functioning, Psychometrika, published online, DOI: 10.1007/s11336-015-9488-3 
#' 
#' @seealso \code{\link[DIFtree]{plot.DIFtree}}, \code{\link[DIFtree]{predict.DIFtree}}, \code{\link[DIFtree]{summary.DIFtree}}
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
#' print(mod)
#' }
#' 
#
#' @exportClass DIFtree
#' @export
#' @importFrom penalized penalized 
#' @importFrom stats binomial coef deviance formula glm predict quantile var


DIFtree <-
function(Y,
                    X,
                    model=c("Rasch","Logistic"),
                    type=c("udif","dif","nudif"),
                    alpha=0.05,
                    nperm=1000,
                    trace=FALSE,
                    penalize=FALSE,
                    ...){
  UseMethod("DIFtree")
}

#' @rdname DIFtree 
#' @method print DIFtree
#' @export

print.DIFtree <-
  function(x, # object of class DIFtree 
           ...){
    
    model <- which(c("Rasch","Logistic") %in% paste(x$call))
    if(model==2){
      type <- which(c("udif","dif","nudif") %in% paste(x$call))
    } else{
      type <- 1
    }
    
    npersons  <- x$persons
    nitems    <- x$items
    
    if(model==2 & type==2){
      dif_items <- unique(c(x$splits[[1]][,"item"],x$splits[[2]][,"item"]))
    } else{
      dif_items <- unique(x$splits[,"item"])
    }
    
    if(model==2 & type==2){
      splits <- lapply(1:2,function(j) x$splits[[j]][,c("item","variable","threshold")])
    } else{
      splits    <- x$splits[,c("item","variable","threshold")]
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
    cat("Number of persons:", npersons, "\n")
    cat("Number of items:", nitems, "\n")
    if(!is.null(dif_items)){
      cat("DIF items:", dif_items, "\n")
    } else{
      cat("DIF items: no DIF item\n")
    }
    cat("\n")
    cat("Overview of executed splits:\n")
    if(!is.null(dif_items)){
      cat("\n")
      if(model==2 & type==2){
        cat("Intercept:\n")
        if(!is.null(splits[[1]])){
          print(splits[[1]])
        } else{
          cat("no split performed")
        }
        cat("\n")
        cat("Slope:\n")
        if(!is.null(splits[[2]])){
          print(splits[[2]])
        } else{
          cat("no split performed")
        }
      } else{
        print(splits)
      }
    } else{
      cat("no split performed")
    }
    invisible(x)
  }  
