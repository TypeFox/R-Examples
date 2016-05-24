
#' Localize errors using a MIP approach.
#' 
#' @section Details:
#' \code{errorLocalizer_mip} uses \code{E} and \code{x} to define a mixed integer problem
#' and solves this problem using \code{lpSolveApi}. This function can be much faster then \code{\link{errorLocalizer}} 
#' but does not return the degeneracy of a solution. However it does return an bonus: 
#' \code{x_feasible}, a feasible solution.
#'
#' 
#'
#' @param E an \code{\link{editset}}, \code{\link{editmatrix}}, or \code{\link{editarray}}
#' @param x named \code{numeric} with data
#' @param weight  \code{numeric} with weights
#' @param maxduration number of seconds that is spent on finding a solution
#' @param verbose verbosity argument that will be passed on to \code{solve} lpSolveAPI
#' @param lpcontrol named \code{list}  of arguments that will be passed on to \code{\link[lpSolveAPI]{lp.control}}. \code{maxduration} will override
#'  \code{lpSolve}'s \code{timeout} argument.
#' @param ... other arguments that will be passed on to \code{solve}.
#' @return list with solution weight \code{w}, \code{logical} \code{adapt} stating what to adapt,  
#'  \code{x_feasible} and the lp problem (an \code{lpExtPtr} object)
#'
#' @seealso \code{\link{localizeErrors}}, \code{\link{errorLocalizer}}, \code{\link{errorLocation}}
#'
#' @references
#'  E. De Jonge and Van der Loo, M. (2012) Error localization as a mixed-integer program in 
#'  editrules (included with the package)
#'
#'  lp_solve and Kjell Konis. (2011). lpSolveAPI: R Interface for
#'  lp_solve version 5.5.2.0. R package version 5.5.2.0-5.
#'  http://CRAN.R-project.org/package=lpSolveAPI
#' @export
errorLocalizer_mip <- function( E
                            , x
                            , weight=rep(1, length(x))
                            , maxduration=600L
                            , verbose="neutral"
                            , lpcontrol = getOption("er.lpcontrol")
                            , ...
                            ){

   vars <- getVars(E)
   E <- as.editset(E)
   DUMP <- FALSE
   
   wp <- perturbWeights(as.vector(weight))
   
   t.start <- proc.time()
   
   #elm <- writeELAsMip(E=E, x=x, weight=wp, ...)
   elm <- as.mip(E=E, x=x, weight=wp, prefix="adapt.",...)
   
   lps <- as.lp.mip(elm)
   # end TODO

   lpcontrol$timeout <- maxduration
   lpcontrol$lprec <- lps
   do.call(lp.control,lpcontrol)

   if (DUMP) write.lp(lps, "test3.lp")
   
   statuscode <- solve(lps)
   degeneracy <- get.solutioncount(lps)
   
   sol <- get.variables(lps)
   # lps may have optimized and removed redundant adapt.variables, so retrieve names of variable...
   names(sol) <- colnames(lps)
   
   #print(sol)
   ## weights are perturbed so objective function is recalculated as sum(weights[adapt])
   # w <- get.objective(lps)
   # get the positions of the adapt.variables
   aidx <- grepl("^adapt\\.", names(sol))
   # split solution in a value and a adapt part
   sol.adapt <- sol[aidx]
   sol.values <- sol[!aidx]
   
   cat.idx <- names(sol.values) %in% names(elm$binvars)
   sol.cat <- asLevels(sol.values[cat.idx])
   
   delta.idx <- grepl("^delta\\.", names(sol.values))
   sol.num <- sol.values[!cat.idx & !delta.idx]
   
   #print(list(sol.cat=sol.cat, sol.num=sol.num))
   
   names(sol.adapt) <- sub("^adapt\\.","",names(sol.adapt))
   
   if (DUMP) write.lp(lps, "test4.lp")
   #print(list(idx=idx, sol=sol))
   adapt <- sapply(x, function(i) FALSE)
   adapt[names(sol.adapt)] <- (sol.adapt > 0)
   #browser()
   
   x_feasible <- x
   idx <- match(names(sol.values), names(x), nomatch=0)
   
   x_feasible[names(sol.num)] <- sol.num 
    
   if (length(sol.cat))
      x_feasible[names(sol.cat)] <- sol.cat
#    if (is.cateditmatrix(E)){
#      x_feasible[idx] <- asLevels(sol.values)
#    } else {
#      x_feasible[idx] <- sol.values
#    }
   

   t.stop <- proc.time()
   duration <- t.stop - t.start
   list( w=sum(weight[adapt])
       , adapt = adapt
       , x_feasible = x_feasible
       , duration = duration
       , maxdurationExceeded = unname(duration[3] >= maxduration)
       , statuscode = statuscode
       , degeneracy = degeneracy
       , lp=lps
       )
}


# scale a numeric vector
scale_fac <- function(x){
   sc <- 1
   m <- max(abs(x),na.rm=TRUE)
   if (is.finite(m) && m >= 1E9 ) sc <- 1/sqrt(m)
   sc   
}

# assumes that E is normalized!
#' Coerces a \code{mip} object into an lpsolve object
#' 
#' \code{as.lp.mip} transforms a mip object into a lpSolveApi object.
#' @param mip object of type \code{mip}.
#' @seealso \code{\link{as.mip}}, \code{\link{make.lp}}
#' @export
as.lp.mip <- function(mip){
#    if (!require(lpSolveAPI)){
#      stop("This function needs lpSolveAPI which can be installed from CRAN")
#    }
   E <- mip$E
   
   A <- getA(E)
   lps <- make.lp(nrow(A), ncol(A))
   
   dimnames(lps) <- dimnames(A)   
   for (v in 1:ncol(A)){
     set.column(lps, v, A[,v])
   }
   
   ops <- getOps(E)
   ops[ops=="=="] <- "="
   ops[ops=="<"] <- "<="
   set.constr.type(lps,types=ops)
   #print(list(lps=lps, objfn=mip$objfn, mip=mip, b=getb(E)))
   set.objfn(lps, mip$objfn)
   set.type(lps, columns=mip$binvars , "binary")
   set.bounds(lps, lower=rep(-Inf, length(mip$numvars)), columns=mip$numvars)
   
   # should improve performance quite a lot: a SOS1 makes bin variables exclusive.
   for (sos in asSOS(colnames(lps))){
     add.SOS( lps, sos$name, 
              type=1, priority=1, 
              columns=sos$columns, 
              weights=sos$weights
            )
   }
   
   b <- getb(E)
   set.constr.value(lps, b)
   lps
}

# splits category names (<variable>:<category>) into variable column groups needed
# for SOS1 constraints
asSOS <- function(vars){
  CAT <- ":\\w+"
  
  idx <- grepl(CAT, vars)
  var <- sub(CAT, "", vars)
  
  sosname <- unique(var[idx])
  sapply(sosname, function(sos){
    columns = which(var == sos)
    list( name=sos
        , columns=columns
        , weights = rep(1, length(columns))
    )
  }, simplify=FALSE)
}

asCat <- function(x, sep=":", useLogicals=TRUE){
  nms <- paste(names(x),x, sep=sep)
  if (useLogicals) {
    idx <- x %in% c("TRUE", "FALSE")
    nms[idx] <- names(x)[idx]
  }
  names(nms) <- names(x)
  nms
}

#' Transform a found solution into a categorical record
#' @keywords internal
asLevels <- function(x){
  vars <- sub(":.+", "", names(x))
  lvls <- sub(".+:", "", names(x))
  logicals <- vars == lvls
  lvls[logicals] <- as.character(x[logicals] > 0)
  names(lvls) <- vars
  lvls[x > 0 | logicals]
}

#testing...

# Et <- editmatrix(expression(
#         p + c == t,
#         c - 0.6*t >= 0,
#         c>=0,
#         p >=0
#         )
#                )
# 
# x <- c(p=755,c=125,t=200)
# 
# errorLocalizer_mip(Et, x)
# 
# Et2 <- editmatrix(expression(
#   p + c == t    
#   ))
# x <- c(p=75,c=125,t=300)
# errorLocalizer_mip(Et2, x)  # random?
# errorLocalizer_mip(Et2, x, weight=c(1,1,1))  # random?
# 
# 
# 
# Es <- c(
#   "age %in% c('under aged','adult')",
#   "maritalStatus %in% c('unmarried','married','widowed','divorced')",
#   "positionInHousehold %in% c('marriage partner', 'child', 'other')",
#   "if( age == 'under aged' ) maritalStatus == 'unmarried'",
#   "if( maritalStatus %in% c('married','widowed','divorced')) !positionInHousehold %in% c('marriage partner','child')"
#   )
# Ec <- cateditmatrix(c(
#   "age %in% c('under aged','adult')",
#   "maritalStatus %in% c('unmarried','married','widowed','divorced')",
#   "positionInHousehold %in% c('marriage partner', 'child', 'other')",
#   "if( age == 'under aged' ) maritalStatus == 'unmarried'",
#   "if( maritalStatus %in% c('married','widowed','divorced')) !positionInHousehold %in% c('marriage partner','child')"
#   ))
# Ec
# r <- c(age = 'under aged', maritalStatus='married', positionInHousehold='child')
# # buildELMatrix(Et,x)
# # buildELMatrix(Ec,r)
#   errorLocalizer_mip(Et, x)
#   errorLocalizer_mip(Ec, r)
# # # asCat(r)
# #  
