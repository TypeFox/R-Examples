#' Solve a System of Polynomial Equations
#'
#' \code{polySolve} solves a system of polynomial equations, specifiable in any of several ways.
#' 
#' @param lhs a mpolyList or character vector of left hand sides
#' @param rhs a mpolyList or character vector of right hand sides
#' @param varOrder variable order (see examples)
#' @param ... stuff to pass to bertini
#' @return an object of class bertini
#' @export polySolve
#' @seealso \code{\link{variety}}, \code{\link{bertini}}
#' @examples
#' \dontrun{
#' 
#' # it can solve linear systems 
#' # (here where the line y = x intersects y = 2 - x)
#' polySolve(c("y", "y"), c("x", "2 - x"), c("x", "y"))
#'
#' # or nonlinear systems
#' polySolve(c("y", "y"), c("x^2", "2 - x^2"), c("x", "y"))
#' 
#' # perhaps an easier specification is equations themselves
#' # with either the " = " or " == " specifications
#' # varOrder is used to order the solutions returned
#' polySolve(c("y = x^2", "y = 2 - x^2"), varOrder = c("x", "y"))
#' polySolve(c("y == x^2", "y == 2 - x^2"), varOrder = c("x", "y"))
#'
#'
#' # mpoly objects can be given instead of character strings
#' lhs <- mp(c("y - (2 - x)", "x y"))
#' rhs <- mp(c("0","0"))
#' polySolve(lhs, rhs, varOrder = c("x", "y"))
#'
#' # if no default right hand side is given, and no "=" or "==" is found,
#' # rhs is taken to be 0's.
#' # below is where the lines y = x and y = -x intersect the unit circle
#' polySolve(c("(y - x) (y + x)", "x^2 + y^2 - 1"))
#'
#' # the output object is a bertini object
#' out <- polySolve(c("(y - x) (y + x)", "x^2 + y^2 - 1"))
#' str(out,1)
#' 
#' # here is the code that was run :
#' cat(out$bertiniCode)
#'
#' # the finite and real solutions:
#' out$finite_solutions
#' out$real_finite_solutions
#'
#'
#'
#' 
#' # example from Riccomagno (2008), p. 399
#' polySolve(c(
#'   "x (x - 2) (x - 4) (x - 3)",
#'   "(y - 4) (y - 2) y",
#'   "(y - 2) (x + y - 4)",
#'   "(x - 3) (x + y - 4)"
#' ))
#'
#' }
#'
polySolve <- function(lhs, rhs, varOrder, ...){

  if(missing(rhs)){
  	
  	if(is.character(lhs)){
      if(all(str_detect(lhs, "=="))){
        split <- strsplit(lhs, " == ")
        lhs <- sapply(split, function(x) x[1])
        rhs <- sapply(split, function(x) x[2])          
      } else if(all(str_detect(lhs, "="))){
        split <- strsplit(lhs, " = ")      
        lhs <- sapply(split, function(x) x[1])
        rhs <- sapply(split, function(x) x[2])                
      } else {
        rhs <- mp(as.character(rep(0, length(lhs))))
        if(length(rhs) == 1){
          rhs <- list(rhs)
          class(rhs) <- "mpolyList"
        }
      }
    } else {
      if(is.mpoly(lhs)){
        rhs <- list(mp("0"))
        class(rhs) <- "mpolyList"
      } else if(is.mpolyList(lhs)){
        rhs <- mp(as.character(rep(0, length(lhs))))    
      }
    }
  }

  if(is.character(lhs)) lhs <- mp(lhs)

  if(is.mpoly(lhs)){
    lhs <- list(lhs)
    class(lhs) <- "mpolyList"
  }

  if(!missing(rhs)){
    if(is.character(rhs)) rhs <- mp(rhs)  
    if(is.mpoly(rhs)){
      rhs <- list(rhs)
      class(rhs) <- "mpolyList"
    }
  }  
  stopifnot(is.mpolyList(lhs) && is.mpolyList(rhs)) 
  gens <- lhs - rhs


  # sort out variables
  vars <- vars(gens)
  
  if(!missing(varOrder) && !all(sort(vars) == sort(varOrder))){
    stop("if varOrder is provided, it must contain all of the variables.", call.=F)
  }
  if(!missing(varOrder)) vars <- varOrder
  
  
  
  # format code
  code <- "\nINPUT\n\nvariable_group vars;\nfunction funNames;\n\nfuns\n\nEND;\n"
  code <- str_replace(code, "vars", paste(vars, collapse = ", "))
  code <- str_replace(code, "funNames", 
    paste(paste("f", 1:length(gens), sep = ""), collapse = ", ")
  )

  funs <- suppressMessages(print(gens, stars = TRUE))
  funs <- str_replace_all(funs, "  ", " ")
  funs <- str_replace_all(funs, "\\*\\*", "^")  
  code <- str_replace(code, "funs",
    paste(
      paste(paste("f", 1:length(gens), sep = ""), " = ", funs, ";", sep = ""),
      collapse = "\n"
    )  
  )

  bertini(code, ...)
}
