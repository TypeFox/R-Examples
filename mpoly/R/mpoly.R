#' Define an mpoly object.
#' 
#' mpoly is the most basic function used to create objects of class
#' mpoly. However, it is not a general purpose function; for that
#' see mp.
#' 
#' @param list a list from which to construct an mpoly object
#' @param varorder (optional) a character vector setting the
#'   intrinsic variable order of the polynomial
#' @return Object of class mpoly.
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @seealso \code{\link{mp}}
#' @export
#' @examples
#' list <- list(
#'   c(x = 1, coef = 1, y = 0),
#'   c(x = 0, y = 1, coef = 2),  
#'   c(y = 1, coef = -6),  
#'   c(z = 1, coef = -3, x = 2),  
#'   c(x = 1, coef = 0, x = 3),
#'   c(t = 1, coef = 4, t = 2, y = 4),
#'   c(x = 1),
#'   c(x = 1),
#'   c(coef = 5),
#'   c(coef = 5),
#'   c(coef = -5)
#' )
#' 
#' mpoly(list) # 3 x  -  4 y  -  3 x^2 z  +  4 y^4 t^3  +  5
#' mpoly(list, varorder = c("y", "z", "t", "x"))
#' 
#' list <- list(  c(x = 5, x = 2, coef = 5, coef = 6, y = 0) )
#' mpoly(list)
#' 
#' 
mpoly <- function(list, varorder){

  ## argument checking
  if(!is.list(list)) stop("input to mpoly must be a list.", call. = FALSE)
  
  if(!all(vapply(list, is.numeric, logical(1)))){
    stop("each element of list must be of type numeric.", call. = FALSE)  
  }
  
  if(any( unlist(lapply(list, function(v) nchar(names(v))), use.names = FALSE) == 0 )){
    stop("each element of list must be named for every element.", call. = FALSE)  
  }    
  
  flatList <- unlist(list)
  flatList <- flatList[names(flatList) != "coef"]
  if(!all(is.wholenumber(flatList)) || any(flatList < 0)){
    stop("degrees must be nonnegative integers.", call. = FALSE)
  }  
  
  
  # give terms without a coefficient a coef of 1
  addCoefIfMissing <- function(v){
    if(any("coef" == names(v))) return(v)
    c(v, coef = 1)
  }
  list <- lapply(list, addCoefIfMissing)
  
  
  
  ## organize 
    
  # remove terms with coef 0
  list <- filterOutZeroTerms(list)
  
  
  # remove 0 degrees, combine like degrees, single coef as rightmost element
  list <- lapply(list, function(v){  	
    # separate vardegs from coefs
    coef_ndx <- which(names(v) == "coef")
    coefs <- v[coef_ndx]    
    v     <- v[-coef_ndx]
  	
  	# combine like degrees (sum)
  	if(length(names(v)) != length(unique(names(v)))) v <- fastNamedVecTapply(v, sum)   
  	
  	# combine like coefficients (product)
  	coefs <- c(coef = prod(coefs))
  	
  	# remove zero degree elements, combine and return  
    c(v[v != 0], coefs)
  })
  
  
  
  ## set intrinsic varorder - done again after 0 degrees are removed
  vars <- unique(names(flatList))
  
  
  # deal with varorder argument
  if(!missing(varorder)){  	
    if( !setequal(vars, varorder) ){
      stop(paste(
        "if specified varorder must be a permutation of",
        paste(vars, collapse = ", ")        
      ), call. = FALSE)	
    }    
    vars <- varorder
  }

  
  # sort variables in terms
  list <- lapply(list, function(v){
    p <- length(v) - 1L
    if(p == 0L) return(v)
    c( (v[1:p])[intersect(vars, names(v[1:p]))], v["coef"] )
  })   
  
   
   
  ## prepare to check if like terms are present
  monomials <- vapply(list, function(v){
  	p <- length(v) - 1 # remove coef on monomials
    paste(names(v[1:p]), v[1:p],  sep = "", collapse = "")
  }, character(1))  
  # e.g. c("x1", "y1", "y1", "x2z1", "y4t3", "x1", "x1", "coef5")
  
  unique_monomials <- unique(monomials)
  
  
  
  ## check if like terms are present and, if so, correct  
  if(length(monomials) != length(unique_monomials)){

    matchedMonomials <- match(monomials, unique_monomials)
    matchedMonomials <- factor(matchedMonomials, levels = 1:max(matchedMonomials))
    ndcs2combine     <- split.default(1:length(list), matchedMonomials)
    
    list <- lapply(ndcs2combine, function(v){
      if(length(v) == 1) return(list[[v]])      
      coef <- sum(vapply(list[v], `[`, double(1), length(list[[v[1]]])))
      v <- list[[v[1]]]
      v["coef"] <- coef
      v
    })
    
    names(list) <- NULL # i.e. list <- unname(list)
  }
  
  
  
  ## combine constant terms
  ## mpoly(list(c(x = 1, coef = 1), c(coef = 1), c(coef = 2)))    
  nonConstantTerms <- fastFilter(isNotLengthOne, list)
  constantTerms    <- fastFilter(isLengthOne, list)
  
  if(length(constantTerms) > 0){
    list <- c(
      nonConstantTerms, 
      list(Reduce(`+`, constantTerms))
    )  
  }
  
  
  
  ## re-organize after like-terms combined
  list <- filterOutZeroTerms(list)

  
  
  ## return classed list
  class(list) <- "mpoly"
  list  
}











filterOutZeroTerms <- function(list){
  
  # a term is zero if any of its coefficients are 0
  hasNonZeroCoef <- function(v) all(v["coef"] != 0)
  
  # remove terms with coef 0
  list <- fastFilter(hasNonZeroCoef, list)
  
  # if all terms have been eliminated, recreate it with a 0 coef
  if(length(list) == 0) list <- list(c(coef = 0))  
  
  # return
  list  
}







isLengthOne    <- function(x) length(x) == 1L
isNotLengthOne <- function(x) length(x) > 1L
fastFilter     <- function(f, x) x[vapply(x, f, logical(1))]







fastNamedVecTapply <- function(x, f, type = double(1)){
  uniqueNames  <- unique(names(x))
  matchedNames <- match(names(x), uniqueNames) # indices
  matchedNames <- factor(matchedNames, levels = 1:max(matchedNames))
  groupIndices <- split.default(1:length(x), matchedNames)
  out <- vapply(groupIndices, function(ndcs) f(x[ndcs]), type)
  names(out) <- uniqueNames
  out
}
# x <- 1:10
# names(x) <- sample(letters[1:3], 10, replace = TRUE)
# fastNamedVecTapply(x, sum)
# tapply(x, names(x), sum)



