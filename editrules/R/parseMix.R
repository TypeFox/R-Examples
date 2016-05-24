MIXOPS <- c("if", "||", "|", "&&", "&")

#' Parse a mixed edit
#'
#' parseMix replaces all numerical edits with a generated dummy boolean variable and returns the resulting categorical
#' edit plus the list of found of numerical edits. These expressions should be handled further by \code{parseCat} and 
#' \code{parseNum}.
#' @param e expression to be parsed
#' @param numid starting number for dummy name generation
#' @return list with categorical expression (\code{cat}), which is normalized,  a numerical expression (code{nums}) 
#' and a negated version of this  expression (code{negNums})
#' @keywords internal
parseMix <- function(e, editname="", numid=0, negate=TRUE){
  
  # should the expressions be returned or should parseCat and parseNum be called on cat and nums?
  
  if (length(e) < 3) return(NULL)
  op <- as.character(e[[1]])
  if (!op %in% MIXOPS) stop("invalid mix")
    
  cat <- e
  nums <- expression()
  numid <- numid
  
  pm <- function(i, neg=negate){
    # rewrite equality as inequalities
    #e[[i]] <- rewriteEq(e[[i]])
    edit <- e[[i]]
    if (length(edit) ==  1) return()
    cmp <- deparse(edit[[1]])
    
    # remove brackets
    if (cmp == "("){
      edit <- edit[[2]]  
    }
    
    if (isNum(edit)){
      
      if (cmp == "=="){
        stop("edit '", deparse(e),"' contains '==' which is not allowed, in a mixed edit")
      }
      
      numid <<- numid + 1
      numvar <- paste(editname, ".l",numid,sep="")
      #replace numeric edit with generated dummy boolean edit name and handle normal form
      dum <- as.name(numvar)
      if (neg){ 
        dum <- quote(!a)
        dum[[2]] <- as.name(numvar)
        edit <- negateEdit(edit)
      } 
      cat[[i]] <<- dum
      nums[numvar] <<- as.expression(edit)
    } else if (!isCat(edit)){
      l <- parseMix(edit, numid=numid, editname=editname, negate=neg)
      cat[[i]] <<- l$cat
      nums[names(l$nums)] <<- l$nums
      numid <<- l$numid
    }  
  }
  
  # don't negate the premisse (term 2)
  pm(2, neg=ifelse(op == "if", !negate, negate))
  # negate the consequent (term 3)
  pm(3, neg=negate)
  
  negNums <- nums
  negNums[] <- sapply(nums, negateEdit)
  
  if (length(negNums)){
    names(negNums) <- paste("!",names(negNums), sep="")
  }
  
  list( cat  = cat # parseCat(cat)
      , nums = nums # lapply(nums, parseNum)
      , negNums = negNums # lapply(orNums, parseNum)
      , numid=numid 
      )
}

# e has to be an edit
negateEdit <- function(e){
  op <- as.character(e[[1]])
  if (op == "!") 
    return(e[[2]])
  op <- switch( op
              , "<"  = ">="
              , "<=" = ">"
              , ">"  = "<="
              , ">=" = "<"
              , "==" = "!="
              , "!=" = "=="
              )
  if (is.null(op)){
    ne <- quote(!a)
    ne[[2]] <- e
    e <- ne
  } else {
    e[[1]] <- as.symbol(op)
  }
  e
}

# e has to be an numerical inequality!
rewriteInEq <- function(e){
  op <- as.character(e[[1]])
  if (op != "!=") return(e)
  eAnd <- quote(a || b)
  eAnd[[2]] <- e
  eAnd[[2]][[1]] <- as.symbol(">")
  eAnd[[3]] <- e
  eAnd[[3]][[1]] <- as.symbol("<")
  eAnd
}

# quick test
# rewriteInEq(quote(x != y + 1))
# 
# a <- negateEdit(quote(x>2))
# a
# negateEdit(a)
# 
# a <- negateEdit(quote(A %in% "a"))
# a
# negateEdit(a)
# 
# pm <- parseMix( quote( if(x>1 && 
#                        x < 10 && 
#                        A %in% c('a1'))  y > 2                     )
#               , editname="e1")
# pm
# 

#pm <- parseMix(quote(if ((x>0)) A))
# pm
