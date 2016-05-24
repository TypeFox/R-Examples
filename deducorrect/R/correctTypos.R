#' Correct records under linear restrictions using typographical error suggestions 
#'
#' This algorithm tries to detect and repair records that violate linear equality constraints by correcting simple typo's as described in Scholtus (2009).
#' The implemention of the detection of typing errors differs in that it uses the restricted Damerau-Levensthein distance. Furthermore it solves a broader class of 
#' problems: the original paper describes the class of equalities: \eqn{Ex=0} (balance edits) and this implementation allows for  \eqn{Ex=a}.
#' 
#' For each row in \code{dat} the correction algorithm first detects if row \code{x} violates the equality constraints of \code{E} taking possible rounding errors into account.
#' Mathematically:
#' \eqn{|\sum_{i=1}^nE_{ji}x_i - a_j| \leq \varepsilon,\quad \forall j }
#'
#' It then generates correction suggestions by deriving alternative values for variables only involved in the violated edits. The correction suggestions must be within a typographical
#' edit distance (default = 1) to be selected. If there are more then 1 solutions possible the algorithm tries to derive a partial solution, otherwise the solution is applied to the data.
#'
#' \code{correctTypos} returns an object of class \code{\link[=deducorrect-object]{deducorrect}} object describing the status of the record and the corrections that have been applied.
#'
#' Inequalities in editmatrix \code{E} will be ignored in this algorithm, so if this is the case, the corrected records
#' are valid according to the equality restrictions, but may be incorrect for the given inequalities.
#'
#' Please note that if the returned status of a record is "partial" the corrected record still is not valid.
#' The partially corrected record will contain less errors and will violate less constraints. 
#' Also note that the status "valid" and "corrected" have to be interpreted in combination with \code{eps}.
#' A common scenario is first to correct for typo's and then correct for rounding errors. This means that in the first
#' step the algorithm should allow for typo's (e.g. \code{eps==2}). The returned "valid"  record therefore may still contain 
#' rounding errors.
#'
#' @export
#' @example ../examples/correctTypos.R
#' @seealso \code{\link{damerauLevenshteinDistance}}
#'
#' @param E \code{editmatrix} or \code{editset}
#' @param dat \code{data.frame} with data to be corrected.
#' @param ... arguments to be passed to other methods.
#' @param fixate \code{character} with variable names that should not be changed.
#' @param cost for a deletion, insertion, substition or transposition.
#' @param eps \code{numeric}, tolerance on edit check. Default value is \code{sqrt(.Machine$double.eps)}. Set to 2 
#' to allow for rounding errors. Set this parameter to 0 for exact checking.
#' @param maxdist \code{numeric}, tolerance used in finding typographical corrections. Default value 1 allows for one error. Used in combination with \code{cost}.
#'
#' @return \code{\link[=deducorrect-object]{deducorrect}} object with corrected data.frame, applied corrections and status of the records.
#' 
#' @references
#' 
#' Scholtus S (2009). Automatic correction of simple typing errors in numerical data with balance edits.
#' Discussion paper 09046, Statistics Netherlands, The Hague/Heerlen.
#' 
#' Damerau F (1964). A technique for computer detection and correction of
#' spelling errors. Communications of the ACM, 7,issue 3
#'
#' Levenshtein VI (1966). Binary codes capable of correcting deletions, insertions, 
#' and reversals. Soviet Physics Doklady 10: 707-10
#'
#' A good description of the restricted DL-distance can be found on wikipedia: http://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance 
#'
correctTypos <- function(E, dat, ...){
    UseMethod('correctTypos')
}

#' @method correctTypos editset
#' @rdname correctTypos
#' @export
correctTypos.editset <- function(E, dat,...){
    correctAndRevert(correctTypos.editmatrix, E, dat, ...)
}


#' @method correctTypos editmatrix
#' @rdname correctTypos
#' @export
correctTypos.editmatrix <- function( E
                        , dat
                        , fixate = NULL
                        , cost = c(1,1,1,1)
                        , eps = sqrt(.Machine$double.eps)
                        , maxdist = 1
                        , ...
                        ){
                        
   stopifnot(is.editmatrix(E), is.data.frame(dat))
   
   # separate equalities and inequalities
   a <- getb(E)
   eq <- getOps(E) == "=="
   F <- E[!eq,]
   E <- E[eq,]
   a <- a[eq]
   vars <- getVars(E)
   
   fixate <- if(is.null(fixate)) {rep(FALSE, length(vars))}
             else vars %in% fixate
   
   #align names of E and dat, beware m contains only constrained, numeric variables at this point
   m <- as.matrix(dat[vars])
   n <- nrow(m)
   
   status <- status(n)
   corrections <- NULL
   
   # only loop over complete records
   cc <- which(complete.cases(m))
	for (i in cc){
      x <- m[i,]
      chk <- getTypoCorrection(E, x, fixate=fixate, eps=eps, maxdist=maxdist)
      
      status[i] <- chk$status
      
      if (chk$status %in% c("valid", "invalid")){
         #nothing we can do...
         next
      }

      cor <- chk$cor
      
      sol <- tree(chk$B, cor[,5])
      if (nrow(sol) > 1){
         # if a correction is valid for all found solutions, then it can be applied
         partialsol <- colSums(sol) == nrow(sol)
         if (any(partialsol)){
            sol[1,] <- partialsol
            status[i] <- "partial"
         }
         else {
            status[i] <- "invalid"
            next
         }
      }
      cor <- cor[sol[1,],,drop=FALSE]
      
      #m[i, cor[,"var"]]  <- cor[,"new"]
      x[cor[,1]] <- cor[,3]
      #TODO if any violatedEdits then solution is always partial
#     print(violatedEdits(F,x) )
      if (all(which(violatedEdits(F, x)) %in% which(violatedEdits(F,m[i,])))){
         m[i,] <- x
      }
      else {
         status[i] <- "invalid"
         next
      }
      # check if record is now valid with the corrections applied
      status[i] <- if (sum(abs(a-getA(E)%*%m[i,]) > eps) == 0) "corrected"
                   else "partial"
                   
      cor <- cbind(row=rep(i, nrow(cor)), cor)
      corrections <- rbind(corrections, cor)      
	}
      
   # recreate data.frame dat in original column order, but with the corrections applied
   corrected <- dat   
   corrected[vars] <- as.data.frame(m)[]
   
   cdf <- data.frame( row=corrections[,1]
                    , variable=vars[corrections[,2]]
                    , old=corrections[,3]
                    , new=corrections[,4]
                    )
    if ( nrow(cdf) == 0 ){
        cdf <- data.frame(
            row = integer(0),
            variable = factor(levels=vars),
            old = numeric(0),
            new = numeric(0)
        )
    }

    return(
        newdeducorrect(
            status = data.frame(status=status),
            corrected = corrected,
            corrections = cdf,
            ...
        )
    )
}

#' Check record validity and suggest typo corrections
#'
#' This function is the working horse for \code{\link{correctTypos}}
#' @keywords internal
#' @param E editmatrix
#' @param x numerical record to be checked
#' @param eps tolerance for an edit to be valid
#' @param maxdist maximum edit distance to be valid as a correction suggestion
#' @return list with members
#' \tabular{ll}{
#' status \tab \cr
#' cor    \tab suggested corrections \cr
#' B      \tab reduced binary editmatrix with violated edits, needed for choosing the suggested corrections\cr
#'}
getTypoCorrection <- function( E, x, fixate=FALSE, eps=sqrt(.Machine$double.eps), maxdist=1){
   ret <- list(status=NA)
   
   a <- getb(E)
   M <- getA(E)
   
   #violated edits (ignoring rounding errors)
   E1 <- (abs(a-M%*%x) > eps)
   
   #non violated edits
   E2 <- !E1
   
   if (all(E2)){
      #record is valid ignoring rounding errors
      ret$status <- "valid"
      return(ret)
   }
   
   B <- M != 0
   # set of variables that are involved in the violated edits
   V1 <- if (any(E1)) colSums(B[E1,,drop=FALSE]) != 0
         else FALSE
               
   # set of variables that are not involved in the non-violated edits and therefore can be edited
   I0 <- if (any(E2)) colSums(B[E2,,drop=FALSE]) == 0
         else TRUE

   # restrict I0 to the set of variables involved in violated edits that can be changed
   I0 <- V1 & I0 & !fixate
   
   if (sum(I0) == 0){
		# cannot correct this error
      ret$status <- "invalid"
		return(ret)
   }
   
   names(I0) <- getVars(E)
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
                                 , dist = damerauLevenshteinDistance(x_i_c[j], x[i])
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
