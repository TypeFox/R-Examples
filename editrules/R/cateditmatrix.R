#' Create an editmatrix with categorical variables
#'
#' \code{cateditmatrix} is an alternative representation of a categorial edit. 
#' The default representation in editrules is \code{\link{editarray}}, but cateditmatrix is useful for
#' transforming and solving categorical edit into a mixed integer programming problem
#'
#' @param x \code{character or expression} with categorical edits
#' @param sep seperator to be used in naming categories within variables.
#' @param env environment in which \code{x} will be evaluated.
#' @return cateditmatrix object, which is a specialized \code{\link{editmatrix}}
#' @keywords internal
cateditmatrix <- function(x, sep=":", env=parent.frame()){
    if (is.editarray(x)) {
      x <- as.character(x)
    }
    
    if ( length(x) == 0 ){
      res <- neweditmatrix(matrix(numeric(0)),ops=character(0),normalized=TRUE)
      class(res) <- c("cateditmatrix", "editmatrix")
      return(res)
    }
    
    edts <- parseEdits(x)
    #names(edts) <- names(x)
    
    catedits <- lapply(edts,parseCat,sep=sep, useLogical=TRUE, env=env)
    catedits <- lapply(catedits, parseCatEdit)
    
    # remove empty rules...
    empty <- sapply(catedits, is.null)
    
    catedits <- catedits[!empty]
    edts <- edts[!empty]
    
    categories <- sort(unique(names(unlist(catedits))))
    categories <- c(categories[categories!="b"],"b")

    A <- matrix( 0
               , ncol=length(categories)
               , nrow=length(catedits)
               , dimnames = list( rules = names(catedits)
                                , cats=categories
                                )
               )
                     
    for (i in seq_along(catedits)){
           A[i,names(catedits[[i]])] <- catedits[[i]]
    }

    b <- A[,ncol(A)]
    A <- A[,-ncol(A), drop=FALSE]
    
    ops <- sapply(edts, function(e){deparse(e[[1]])})
    ops <- ifelse(ops %in% c("if", "||"), "<=", "==")

    E <- as.editmatrix( A
                      , b
                      , ops
                      # create named TRUE vector for dummy variables
                      , binvars=sapply(categories,`!=`, "b")
                      )
    class(E) <- c("cateditmatrix", "editmatrix")
    E
}

#' Coerce an cateditmatrix to a \code{character} vector
#'
#' Derives readable editrules from an editmatrix.
#' 
#' @method as.character cateditmatrix
#'
#' @param x cateditmatrix object to be printed
#' @param asIfStatement if \code{TRUE} the representation will be an \code{if} statement.
#' @param ... further arguments passed to or from other methods.
#' @keywords internal
as.character.cateditmatrix <- function(x, asIfStatement=FALSE, ...){
  # empty cateditmatrix?
  if (length(x) == 0) 
    return(str(x))
  
  if (isTRUE(asIfStatement)){
    A <- getA(x)
    ops <- getOps(x)
    
    ifc <- (A > 0 & ops == "<=")
    thenc <- (A < 0 | ( A > 0 & ops == "=="))

    vars <- getVars(x, "var")
    cats <- getVars(x, "cat")
    
    #generate %in% statement
    inclause <- function(idx, collapse=NULL){
      vc <- split(cats[idx], vars[idx])
      
      vc <- sapply( names(vc), function(var){
        cats <- vc[[var]]
        if (length(cats) == 1){ #this is for a logical variable
          if (cats %in% TRUE)
            return(var)
          paste(var," == '",cats,"'", sep="")  # cosmetic, for one category we generate an "==" statement
        } else {
          cats <- paste("'",cats,"'", sep="", collapse=",")
          paste(var," %in% c(",cats,")", sep="")
        }
      })
      paste(vc, collapse=collapse)
    }
    
    catedits <- rownames(x)
    for (i in 1:length(catedits)){
      if (any(thenc[i,])){
        thenvars <- inclause(thenc[i,], collapse=" || ")
        if (any(ifc[i,])){
          ifvars <- inclause(ifc[i,], collapse=" && ")
          catedits[i] <- paste("if (",ifvars,") ",thenvars, sep="")
        } else {
          catedits[i] <- thenvars
        }
      } else {
        catedits[i] <- paste("!(", inclause(ifc[i,]), ")",sep="",collapse=" || ")
      }
    }
    names(catedits) <- rownames(x)
    log <- sapply(which(cats == "TRUE"), function(v){paste(vars[v]," %in% c(TRUE,FALSE)",sep="")})
    c(catedits, log)
  } else {
    class(x) <- "editmatrix"
    as.character(x)
  }
}

toCat <- function(categories, sep=":"){
  vars <- sub(":.+", "", categories)
  lvls <- sub(".+:", "", categories)
  names(lvls) <- vars
  split(lvls, vars)
}

#' parse categorial edit

#' @param e \code{expression} with a single edit
#' @return named \code{numeric} with coefficients
#' @keywords internal
parseCatEdit <- function(el){
  #el <- parseCat(e, useLogical=TRUE)
  if (any(is.na(el))){
    if (length(el) == 1) return(NULL)
    #browser(expr={length(el) == 1})
    val <- rep(1, length(el)+1)
    names(val) <- c(names(el), "b")    
    #print(el)
  } else {
    vars <- gsub(":.+","",names(el))
    # coefficients in form 
    val <- ifelse(el, 1, -1)
    m <- tapply(val, vars, max)
    b <- sum(m[m>0]) - 1
    val <- c(val, b=b)
  }
  val
}

# ### examples....
# 
# #civilStatusLevels <- c("married","unmarried","widowed","divorced")
# # 
# x <- c( "if (positionInHousehold == 'marriage partner') civilStatus == 'married'"
#       , "if (age == '< 16') civilStatus=='unmarried'"
# #      , "civilStatus %in% civilStatusLevels" #that looks magical, but civilstatusLevels is evaluated
#       , "if (pregnant) gender == 'female'"
#       , "if (nace %in% c('A','B')) valid==TRUE"
#       , "gender %in% c('male','female')"
#       )
# 
# (E <- cateditmatrix(x))
# attributes(E)
