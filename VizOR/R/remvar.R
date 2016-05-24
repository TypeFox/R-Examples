##' Removes variables from a regression model formula
##'
##' This is a utility function, intended to support \code{fastbw}
##' @title remvar
##' @param f A model formula
##' @param vars A character vector giving names of variables to remove from \code{f}
##' @return A modified formula, omitting the specified variables
##' @author David C. Norris
##' @export remvar
remvar <- function(f, vars){
  ## TODO: Consider making this a 'private' function of VizOR
  if(length(f)==1){ # Base case of atomic f
    if(as.character(f) %in% vars)
      NULL
    else
      f
  } else { # f is compound
    if(f[[1]]=="~"){ # remove vars from only the RHS of a '~' expression
      rhs <- remvar(f[[3]], vars)
      f[[3]] <- remvar(f[[3]], vars)
      if(is.null(rhs))
        f[[3]] <- 1
      else
        f[[3]] <- rhs
      f
    } else if(f[[1]]=="(") { # simply return the contents of a grouped expression
      f <- remvar(f[[2]], vars)
      f
    } else if(f[[1]]=="+") { # pluck any vars from a sum, retaining others
      ## If either of the terms is one of the 'vars', then we return the other term
      lhs <- remvar(f[[2]], vars)
      rhs <- remvar(f[[3]], vars)
      if(is.null(lhs))
        rhs
      else if(is.null(rhs))
        lhs
      else
        substitute(lhs + rhs, list(lhs=lhs, rhs=rhs))
    } else { # I treat other operators (including functions of any arity) like '*'
      for(i in 2:length(f)){
        f.i <- remvar(f[[i]], vars)
        if(is.null(f.i)) # If any operand/arg is NULL..
          return(NULL)   # ..then return NULL.
        f[[i]] <- f.i    # Otherwise, transform elementwise.
      }
      f
    }
  }
}
