#' Rewrite an editset and reported values into the components needed for a mip solver
#' 
#' @param E \code{\link{editset}} or any object that is coercable to an editset.
#' @param x named \code{list} or \code{vector} with data
#' @param weight vector with weights of the variable in the same order as x
#' @param M maximum allowed difference between reported value and corrected value
#' @param epsilon offset needed for writing a strict inequality into a an inequality 
#' @return list with an editmatrix, objfn, binvars, numvars, M and epsilon 
#' @keywords internal
writeELAsMip <- function( E
                      , x
                      , weight = rep(1, length(x))
                      , M = 1e7
                      , epsilon = 1e-3
                      , ...
                      ){
  E <- as.editset(E)
  el.E <- NULL
    
#   soft <- is.finite(editweight)
#   if (any(soft)){
#     soft.E <- E[soft,]
#     soft.weights <- editweight[soft]
#     
#     #TODO process softedits into el.E
#     soft.num <- softEdits(soft.E$num, xlim, prefix=".soft.")
#     
#     #TODO column with diag(1, nrow(soft.E$mixcat))
#     soft.cat <- NULL
#     #soft.cat <- softEdits(cateditmatrix(soft.E$mixcat),xlim,prefix=".softcat.")
#     el.E <- c(soft.num, soft.cat, el.E)
#     E <- E[!soft,]
#   }
    
  # num part
  num.vars <- getVars(E, type="num")
  
  if (!is.null(num.vars)){
    num.idx <- match(num.vars, names(x))
    num.x <- diag(1, nrow=length(num.vars))
    dimnames(num.x) <- list(num.vars,num.vars)
    num.x0 <- unlist(x[num.idx])
    # create an editmatrix x_i == x^0_i
    num.E <- as.editmatrix(num.x, num.x0)
    num.se <- softEdits(num.E, "adapt.")
    el.E <- c(num.se, E$num, el.E)
  }

  # cat part
  cat.vars <- getVars(E, type="cat")
  if ( length(cat.vars) > 0 ){
    cat.idx <- match(cat.vars, names(x))
    cat.x_0 <- unlist(x[cat.idx])
    

    cat.A <- diag(1, nrow=length(cat.x_0))
    cat.A <- cbind(cat.A, cat.A)
    #browser()    
    colnames(cat.A) <- c(asCat(cat.x_0), paste("adapt.", cat.vars, sep=""))
    
    # check for non existing levels (including NA's)
    invalidCats <- !(asCat(cat.x_0, useLogicals=FALSE) %in% getlevels(E$mixcat))    
    if (any(invalidCats)){ # remove invalid categories otherwise they will turn up in the resulting editmatrix...
      cat.A <- cat.A[,-which(invalidCats), drop=FALSE]
    }
    cat.b <- rep(1, nrow(cat.A))
    
    cat.se <- as.editmatrix(cat.A, b=cat.b)
    el.E <- c(cat.se, cateditmatrix(E$mixcat), el.E)
  }
  
  # mix part  
  mix.E <- editmatrix(invert(as.character(E$mixnum)))
  mix.vars <- getVars(mix.E)
  if (!is.null(mix.vars)){
    mix.idx <- match(mix.vars, names(x))
    mix.se <- softEdits(mix.E, prefix="")
    el.E <- c(mix.se, el.E)
  }
  
#  el.E <- c(mix.se, cat.se, num.se, E$num, cateditmatrix(E$mixcat))     
  lt <- getOps(el.E) == "<"
    
  el.vars <- getVars(el.E)
  el.binvars <- sapply(el.vars, is.character)
  el.binvars[el.vars %in% num.vars] <- FALSE
  g <- grepl("delta.", el.vars, fixed=TRUE)
  #print(g)
  el.binvars[g] <- FALSE
  
  objfn <- sapply(el.vars, function(v) 0)
  adapt.idx <- grep("^adapt\\.", el.vars)
  adapt.nms <- names(adapt.idx) <- sub("^adapt\\.", "", el.vars[adapt.idx])
  
  objfn[adapt.idx] <- weight[match(adapt.nms, names(x))]

#   if (any(soft)){
#      soft.idx <- grep("^\\.soft", el.vars)
#      objfn[soft.idx] <- (1-lambda) * soft.weights
#   }
  structure(
    list( E = el.E
        , objfn = objfn
        , binvars = which(el.binvars)
        , numvars = match(num.vars, el.vars)
        , M = M
        , epsilon = epsilon
        ),
    class="mip"
  )
}

buildELMatrix <- writeELAsMip

# E <- editset(expression(
#   x < y,
#   y < z,
#   x < z,
#   a %in% c(TRUE, FALSE),
#   if (a) x > 1
#   ))
# 
# editsetToMip(E)
#testing...

# E <- editset(expression(
#          if (x>0) y > 0
#       ,  maritalstatus %in% c("married", "single")
#       ,  if (maritalstatus == "married") age >= 17
#       ))
# # 
# x <- list(x = 1, y = -1, age=16, maritalstatus="married")
# # #x <- list(x = 1, y = -1, age=16, maritalstatus=NA)
# # # e <- expression( pregnant %in% c(TRUE, FALSE)
# # #                , gender %in% c("male", "female")
# # #                , if (pregnant) gender == "female"
# # #                )
# # # 
# # # cateditmatrix(e)
# # checkXlim(list(age=c(0,200)), x)
# # 
# buildELMatrix(E, x)# -> B
# #errorLocalizer.mip(E, x=x,, xlim=list(age=c(0,200)))
