#' Correct typos in restricted numeric data
#'
#' Attempt to fix violations of linear (in)equality restrictions imposed on a
#' record by replacing values with values that differ from the original values
#' by typographical errors.
#' 
#' @section Details:
#' 
#' The algorithm works by proposing candidate replacement values and checking
#' whether they are likely to be the result of a typographical error. A value is
#' accepted as a solution when it resolves at least one equality violation. An
#' equality restriction \code{a.x=b} is considered satisfied when
#' \code{abs(a.x-b)<eps}. Setting \code{eps} to one or two units of measurement
#' allows for robust typographical error detection in the presence of 
#' roundoff-errors.
#' 
#' The algorithm is meant to be used on numeric data representing integers.
#' 
#'
#' @param dat An R object holding numeric (integer) data.
#' @param x An R object holding linear data validation rules
#' @param ... Options to be passed to \code{\link[stringdist]{stringdist}}
#'   which is used to determine the typographic distance between the original
#'   value and candidate solutions. By default, the optimal string alignment distance
#'   is used, with all weights equal to one.
#'
#' @return \code{dat}, with values corrected.
#'
#' @references 
#' \itemize{
#' \item{The first version of the algorithm was described by S. Scholtus (2009). Automatic
#' correction of simple typing errors in numerical data with balance edits. Statistics 
#' Netherlands, Discussion Paper \href{http://www.cbs.nl/en-GB/menu/methoden/onderzoek-methoden/discussionpapers/archief/2009/2009-46-x10-pubpdf.htm}{09046}
#' }
#' \item{The generalized version of this algorithm that is implemented for this package is
#' described in M. van der Loo, E. de Jonge and S. Scholtus (2011). Correction of rounding, typing and sign errors with the deducorrect package.
#' Statistics Netherlands, Discussion Paper \href{http://www.cbs.nl/nl-NL/menu/methoden/onderzoek-methoden/discussionpapers/archief/2011/correction-of-rounding-typing-and-sign-errors-with-the-deducorrect-package.htm}{2011019}
#' }
#' }
#'
#' @example ../examples/correct_typos.R
#'
#' @export
setGeneric("correct_typos", function(dat, x,...) standardGeneric("correct_typos"))

#' @rdname correct_typos
#' 
#' @param eps \code{[numeric]} maximum roundoff error
#' @param maxdist \code{[numeric]} maximum allowd typographical distance
#' @param fixate \code{[character]} vector of variable names that may not be changed
#' 
#' 
setMethod("correct_typos", c("data.frame","validator")
  , function(dat,x,fixate=NULL, eps=1e-8,maxdist=1, ...){

  lc <- x$linear_coefficients()
  a <- lc$b
  # separate equalities and inequalities
  eq <- lc$operators == "=="
  F <- x[!eq,]  # inequalities
  E <- x[eq,]   # equalities
  a <- a[eq]
  vars <- validate::variables(x)
  
  fixate <- if( is.null(fixate) ) {
      rep(FALSE, length(vars))
    } else { 
      vars %in% fixate
    }
  
  # align names of x and dat, beware m contains only constrained, numeric
  # variables at this point
  m <- as.matrix(dat[vars])
  n <- nrow(m)
  
  # only loop over complete records
  cc <- which(complete.cases(m))
  for (i in cc){
    r <- m[i,]
    chk <- getTypoCorrection(E, r, fixate=fixate, eps=eps, maxdist=maxdist, ...)
    
    if (chk$status %in% c("valid", "invalid")){ #nothing we can do...
      next
    }
  
    cor <- chk$cor
    
    sol <- tree(chk$B, cor[,5])
    if (nrow(sol) > 1){
       # if a correction is valid for all found solutions, then it can be applied
       partialsol <- colSums(sol) == nrow(sol)
       if (any(partialsol)){
          sol[1,] <- partialsol
     #     status[i] <- "partial"
       }
       else {
     #     status[i] <- "invalid"
          next
       }
    }
    cor <- cor[sol[1,],,drop=FALSE]
    
    
    r[cor[,1]] <- cor[,3]
    
    # only accept solutions that do not violate any new inequality restrictions 
    v1 <- values(confront(vec2df(r),F))
    i1 <- if (length(v1)>0) which(!v1) else integer(0)
    v2 <- values(confront(vec2df(m[i,]),F))
    i2 <- if (length(v2)>0) which(!v2) else integer(0)
    if (all(  i1 %in% i2 ) ){
       m[i,] <- r
    }
    else {
       next
    }
    # check if record is now valid with the corrections applied
    cor <- cbind(row=rep(i, nrow(cor)), cor)
  }

  # recreate data.frame dat in original column order, but with the corrections applied
  dat[vars] <- as.data.frame(m)[]
  dat
})


vec2df <- function(x){
  setNames(as.data.frame(matrix(x,nrow=1)),names(x))
}


getTypoCorrection <- function( E, x, fixate=FALSE, eps=1e-8, maxdist=1,...){
  eps2 <- 1e-8
  ret <- list(status=NA)
  L <- E$linear_coefficients()   
  a <- L$b
  M <- L$A
  
  # violated edits (ignoring rounding errors)
  E1 <- (abs(a-M%*%x) > eps)
  
  #non violated edits
  E2 <- !E1
  
  if (all(E2)){
    #record is valid ignoring rounding errors
    ret$status <- "valid"
    return(ret)
  }
   
  B <- abs(M) > eps2
  # set of variables that are involved in the violated edits
  V1 <- if (any(E1)) colSums(B[E1,,drop=FALSE]) != 0 else FALSE
             
  # set of variables that are not involved in the non-violated edits and therefore can be edited
  I0 <- if (any(E2)) colSums(B[E2,,drop=FALSE]) == 0 else TRUE
  
  # restrict I0 to the set of variables involved in violated edits that can be changed
  I0 <- V1 & I0 & !fixate
   
  if (sum(I0) == 0){
    # cannot correct this error
    ret$status <- "invalid"
    return(ret)
  }
   
  names(I0) <- variables(E)
  names(x) <- NULL
  # retrieve correction canditates for variables that can be changed
  cor <- lapply( which(I0)
              , function(i){
                   # edits valid for current variable v_i
                   eqs <- E1 & (B[,i])
                   # correction candidates
                   #TODO check if solution has to be rounded!!!)
                   x_i_c <- ( (a[eqs]-(M[eqs,-i, drop=FALSE] %*% x[-i])) / (M[eqs,i]))
                   # count their numbers
                   kap <- table(x_i_c)
                   x_i_c <- as.integer(rownames(kap))
                   kap <- as.integer(kap)
                   # and retrieve their distance from the current x[i]
                   sapply( seq_along(kap)
                         , function(j){
                              c( var = i
                               , old = x[i]
                               , new = x_i_c[j]
                               , dist = stringdist::stringdist(x_i_c[j], x[i],...)
                               , kappa = kap[j]
                               )
                           }
                         )
                }
              )
  cor <- t(do.call(cbind,cor))
  # filter out the corrections that have dist > maxdist
  valid <- cor[,4] <= maxdist
  
  if (sum(valid) == 0){
    # cannot correct this error
    ret$status <- "invalid"
    return(ret)
  }
   
  cor <- cor[valid,,drop=FALSE]
  # optimization matrix
  B <- B[E1,cor[,1], drop=FALSE] != 0
  ret$cor <- cor
  ret$B <- B
  ret$status <- "partial"
  ret
}

#' Solve an optimization problem using a tree algorithm as described in Scholtus (2009)
#' @keywords internal
#' @param B binary matrix with suggested corrections per violated edit
#' @param kappa frequency of suggested corrections
#' @param delta \code{logical} vector with partial solution (starts with NA)
#' @param sol current best solution. (starts with null)
#' 
#' @return sol
tree <- function( B
                , kappa
                , delta=as.logical(rep(NA, ncol(B)))
                , sol = NULL
                ) {
   if (any(is.na(delta))){
      i_t <- match(NA,delta); # first element of partial solution that is not determined
      
      # leftnode delta_i_t == FALSE
      delta[i_t] <- FALSE
      sol <- tree(B, kappa, delta, sol)
      
      # rightnode  delta_i_t == TRUE
      # set other corrections involved in this edit to FALSE
      # edits involved in i_t
      E2 <- B[,i_t]
      delta[colSums(B[E2,,drop=FALSE]) > 0] <- FALSE
      delta[i_t] <- TRUE
      sol <- tree(B, kappa, delta, sol)
   }
   else {
      value = kappa%*%delta
      delta  <- matrix(delta, nrow=1)
      if (is.null(sol)){
         sol <- delta
      }
      else {
         vals <- kappa %*% sol[1,]
         if (vals < value){
            sol <- delta
         }
         else if (vals == value){
            sol <- rbind(sol, delta)
         }
      }
   }
   sol
}


