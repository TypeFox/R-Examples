#' Reorder a multivariate polynomial.
#'
#' This function is used to set the intrinsic order of a multivariate
#' polynomial. It is used for both the in-term variables and the
#' terms.
#'
#' @param x an object of class mpoly
#' @param varorder the order of the variables
#' @param order a total order used to order the terms
#' @param ... additional arguments
#' @usage \method{reorder}{mpoly}(x, varorder = vars(x), order, ...)
#' @return An object of class mpoly.
#' @export
#' @examples
#'	
#' list <- list(
#'   c(x = 1, y = 2, z = 1, coef = 1),
#'   c(x = 2, y = 0, z = 2, coef = 1),
#'   c(x = 0, y = 0, z = 2, coef = 1),
#'   c(x = 3, y = 0, z = 0, coef = 1)
#' )
#' (p <- mpoly(list)) # -> x y^2 z  +  x^2 z^2  +  z^2  +  x^3
#' reorder(p) # -> x y^2 z  +  x^2 z^2  +  z^2  +  x^3
#' reorder(p, varorder = c("x","y","z"), order = "lex")
#'     # -> x^3  +  x^2 z^2  +  x y^2 z  +  z^2
#' reorder(p, varorder = c("x","y","z"), order = "glex")
#'     # -> x^2 z^2  +  x y^2 z  +  x^3  +  z^2
#' reorder(p, varorder = c("x","y","z"), order = "grlex")
#'     # -> x y^2 z  +  x^2 z^2  +  x^3  +  z^2
#' 
#' reorder(mp("x + 1"), varorder = c("y","x","z"), order = "lex")
#' reorder(mp("x + y"), varorder = c("y","x","z"), order = "lex")
#' reorder(mp("x y + y + 2 x y z^2"), varorder = c("y","x","z"))
#' reorder(mp("x^2 + y x + y"), order = "lex")
#'
#'
reorder.mpoly <- function(x, varorder = vars(x), order = "lex", ...){
  
  ## set basic quantities
  vars <- vars(x)
  p <- length(vars)
  n <- length(x)
  
  
  
  
  ## argument checking  
  # if(missing(varorder) && missing(order)) return(x)
  
  if(!missing(varorder)){
    if(!all(vars %in% varorder)){
      error <- paste(
        "if specified, varorder must contain all computed vars - ",
        paste(vars, collapse = ", "),
        sep = ""
      )
      stop(error, call. = FALSE)
    }
    
    # reduce varorder to size of vars
    varorder <- varorder[varorder %in% vars]    
    
    ## reorder in-term variables
    x <- lapply(x, function(v){
      if(length(v) == 1) return(v) # constant case
      v <- v[intersect(c(varorder, "coef"), names(v))]
    })
    class(x) <- "mpoly"    
  }  
  
  if(missing(varorder) && !missing(order)){
    message <- paste(
      "using variable ordering - ",
      paste(vars, collapse = ", "),
      sep = ""
    )
    message(message)
  }
  

  
  
  ## argument check - order
  match.arg(order, c("lex","glex","grlex"))
  
  
  
  
  
  ## order == "lex"  
  if(order == "lex"){
  	
  	if(n == 1) return(x)  # single term polynomial
  	
  	# add zeros
    l <- lapply(x, function(v){
      z <- rep(0, p + 1)
      names(z) <- c(varorder, "coef")
      z[names(v)] <- v
      z
    })
    
    # construct matrix and sort appropriately
    m <- matrix(unname(unlist(l)), nrow = length(l), ncol = p + 1, byrow = TRUE)
    dimnames(m) <- list(1:nrow(m), c(varorder, "coef"))
    for(k in p:1) m <- m[order(m[,k], decreasing = TRUE),]	
    
    
    # split into list and add names
    list4mpoly <- unname(lapply(split(m, 1:n), function(v){
      names(v) <- c(varorder, "coef")
      v
    }))
    
    # mpoly it and return
    return( mpoly(list4mpoly, varorder = varorder) )  
  }
  
  
  
  
  
  ## order == "glex"  
  if(order == "glex"){
  	
  	if(n == 1) return(x)  # single term polynomial
  	
  	# add zeros
    l <- lapply(x, function(v){
      z <- rep(0, p + 1)
      names(z) <- c(varorder, "coef")
      z[names(v)] <- v
      z
    })
    
    # construct matrix and sort appropriately
    m <- matrix(unname(unlist(l)), nrow = length(l), ncol = p + 1, byrow = TRUE)
    dimnames(m) <- list(1:nrow(m), c(varorder, "coef"))
    for(k in p:1) m <- m[order(m[,k], decreasing = TRUE),]	
    
    m <- m[order(apply(m[, 1:p, drop = FALSE],1,sum), decreasing = TRUE),]    
    
    # split into list and add names
    list4mpoly <- unname(lapply(split(m, 1:n), function(v){
      names(v) <- c(varorder, "coef")
      v
    }))
    
    # mpoly it and return
    return( mpoly(list4mpoly, varorder = varorder) )
  }  
  
  
  
  
  
  ## order == "grlex"  
  if(order == "grlex"){
  	
  	if(n == 1) return(x)  # single term polynomial
  	
  	# add zeros
    l <- lapply(x, function(v){
      z <- rep(0, p + 1)
      names(z) <- c(varorder, "coef")
      z[names(v)] <- v
      z
    })
    
    # construct matrix and sort appropriately
    m <- matrix(unname(unlist(l)), nrow = length(l), ncol = p + 1, byrow = TRUE)
    dimnames(m) <- list(1:nrow(m), c(varorder, "coef"))
    for(k in 1:p) m <- m[order(m[,k]),]	
    m <- m[order(apply(m[, 1:p, drop = FALSE],1,sum), decreasing = TRUE),]    
    
    # split into list and add names
    list4mpoly <- unname(lapply(split(m, 1:n), function(v){
      names(v) <- c(varorder, "coef")
      v
    }))
    
    # mpoly it and return
    return( mpoly(list4mpoly, varorder = varorder) )    
  }    

}










