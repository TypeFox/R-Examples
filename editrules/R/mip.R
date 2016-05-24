#' Write an editset into a mip representation 
#' 
#' Writes an editset or an object coercable to an editset as a mip problem. 
#' @param E an \code{link{editset}} or an object that is coerciable to an 
#' \code{editset}
#' @param x named \code{list}/\code{vector} with variable values
#' @param weight reliability weights for values of \code{x}
#' @param M Constant that is used for allowing the values to differ from \code{x}
#' @param epsilon Constant that is used for converting '<' into '<='
#' @param prefix prefix for dummy variables that are created
#' @param ... not used 
#' @return a mip object containing al information for transforming it 
#' into an lp/mip problem
#' @export
as.mip <- function( E, x=NULL, weight=NULL, M=1e7, epsilon=1e-3, prefix="delta."
                  , ...){
  
  #TODO add ruleWeights = c(Inf, ..., 1, 1, NA) etc. check with `is.finite`
  #TODO add lambda used in objfn of softEdits
  #TODO add objfn specification (needs reordering in function so it has the same order)
  
  E <- as.editset(E)
  objfn <- NULL
  
  E_mip = c( E$num
           , cateditmatrix(E$mixcat)
           , softEdits(editmatrix(invert(as.character(E$mixnum))), prefix="") 
           )
  
  if (!missing(x)){
    if (is.null(weight)){
      weight <- rep(1.0, length(x))
    } else {
      stopifnot(length(weight) == length(x))
    }
    
    names(weight) <- paste0(prefix, names(x))
    # create expression vector with var == value
    expr <- as.expression( lapply(names(x)
                         , function(v){ 
                              substitute( var == value
                                        , list(var=as.symbol(v), value=x[[v]])
                                        )
                         }))
    
    isna <- sapply(x, is.na)
    
    num_vars <- getVars(E, type="num")
    nums <- !isna & names(x) %in% num_vars
    
    if (any(nums)) { 
      A0_num <- editmatrix(expr[nums])
      rownames(A0_num) <- names(x)[nums]
      A0_num <- softEdits(A0_num, prefix=prefix, postfix="0",M=M)
    } else A0_num <- NULL
    
    cat_vars <- getVars(E, type="cat")
    cats <- !isna & names(x) %in% cat_vars
    
    if (any(cats)){
      A0_cat <- cateditmatrix(expr[cats])
      rownames(A0_cat) <- names(x)[cats]
      A0_cat <- softEdits(A0_cat, prefix=prefix, postfix="0")
    } else A0_cat <- NULL
        
    E_mip <- c(E_mip, A0_num, A0_cat)
    if (any(isna)){
      E_mip <- c( E_mip
                , editmatrix(paste0(prefix, names(x)[isna], "==", 1))
                )
    }
    
  }
  
  # replace strict inequalities...
  A <- getA(E_mip)
  b <- getb(E_mip)
  ops <- getOps(E_mip)
  lt <- ops == "<"
  b[lt] <- b[lt] - epsilon
  ops[lt] <- "<="
  E_mip = as.editmatrix(A=A, b=b, ops=ops)
  
  vars <- colnames(A)
  objfn <- sapply(vars, function(v) 0)
  objfn[names(weight)] <- weight
  
  numvars = sapply(vars, `%in%`, getVars(E, type="num"))
  binvars = !numvars
  
  structure(
    list( E = E_mip
        , objfn = objfn
        , binvars = which(binvars)
        , numvars = which(numvars)
        , M = M
        , epsilon = epsilon
    ),
    class="mip"
  )
}

#' @method print mip
print.mip <- function(x, ...){
  print.editmatrix(x$E, textOnly=T)
  if (!is.null(x$objfn)) {
    idx <- which(x$objfn != 0)
    if (length(idx)) {
      of <- paste0(x$objfn, "*", colnames(x$E))[idx]
      of <- paste(of, collapse=" + ")    
      cat("objective function = min: ", of, "\n")
    }
  }
}

# # quick test 
# E <- editset(c(r1="x > 1","y >= x", r2="if (x>1) y> 2", r3="A %in% c('a', 'b')"))
# as.mip(E)
