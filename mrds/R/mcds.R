#' MCDS function definition
#'
#' Creates model formula list for multiple covariate distance sampling using
#' values supplied in call to \code{\link{ddf}}
#'
#'
#' @param formula formula for scale function
#' @param key string identifying key function (currently either "hn"
#'   (half-normal),"hr" (hazard-rate), "unif" (uniform) or "gamma" (gamma
#'   distribution)
#' @param adj.series string identifying adjustment functions cos (Cosine), herm
#'   (Hermite polynomials), poly (simple polynomials) or NULL
#' @param adj.order vector of order of adjustment terms to include
#' @param adj.scale whether to scale the adjustment terms by "width" or "scale"
#' @param adj.exp if TRUE uses exp(adj) for adjustment to keep f(x)>0
#' @param shape.formula formula for shape function
#' @return A formula list used to define the detection function model
#'   \item{fct}{string "mcds"} \item{key}{key function string}
#'   \item{adj.series}{adjustment function string} \item{adj.order}{adjustment
#'   function orders} \item{adj.scale}{adjustment function scale type}
#'   \item{formula}{formula for scale function} \item{shape.formula}{formula
#'   for shape function}
#' @author Jeff Laake; Dave Miller
#' @keywords utility
mcds <- function(formula=NULL, key=NULL, adj.series=NULL, adj.order=c(NULL),
                 adj.scale="width", adj.exp=FALSE, shape.formula=~1){

  if(is.null(formula)&&key!="unif"){
    stop("Missing formula needed for scale")
  }

  if(class(formula)!="formula"){
    if(class(try(as.formula(formula)))=="formula"){
      formula <- as.formula(formula)
    }else{
      stop("Invalid formula")
    }
  }

  if(!is.null(shape.formula)){
    if(class(shape.formula)!="formula"){
      if(class(try(as.formula(shape.formula)))=="formula"){
        shape.formula <- as.formula(shape.formula)
      }else{
        stop("Invalid shape.formula")
      }
    }
  }

  key <- match.arg(key,c("hn","hr","unif","gamma","th1","th2"))
  if(key%in%c("hn","unif")){
    shape.formula <- NULL
  }

  # What to do if we have adjustment terms
  if(!is.null(adj.series)){
    adj.series <- match.arg(adj.series,c("cos","herm","poly"))
    adj.scale <- match.arg(adj.scale,c("width","scale"))
    if(key=="unif" ){
      if(adj.scale=="scale"){
        message("Setting adj.scale to width for uniform key\n")
      }
      adj.scale <- "width"
    }
    adj.check.order(adj.series,adj.order,key)
  }

  return(list(fct           = "mcds",
              formula       = formula,
              shape.formula = shape.formula,
              key           = key,
              adj.series    = adj.series,
              adj.order     = adj.order,
              adj.scale     = adj.scale,
              adj.exp       = adj.exp))

}
