CATCMP <- c("==", "!=", "%in%")

#' Parse a categorical edit expression 
#'
#' @param x a valid R expression
#' @param val logical (scalar)
#' @param edit logical (vector)
#' @param sep edit separator
#' @param useLogical (logical), should logicals be treated as a factor or as a logical?
#' @keywords internal
parseCat <- function(x, val=NA, edit=logical(0), sep=":", useLogical=FALSE, env=parent.frame()){
    if ( length(x) == 1 ) {
       # corner case: the always FALSE edit (array must be TRUE at every category)
       if ( is.na(val) && !x[[1]] ) return(NULL)
       if (is.logical(x)){
         if (val == x) { return(edit)
         }else { 
           # if this happens the statement is always true, so delete it...
           return (logical())
         }
       }
       var <- if (useLogical) as.character(x)
              else paste(x,"TRUE",sep=sep)
       edit[var] <- val
       return(edit)
    }
    op <- as.character(x[[1]])
    if ( op == "if" ){
        edit <- parseCat(x[[2]],TRUE,  edit, sep, useLogical, env=env)
        edit <- parseCat(x[[3]],FALSE, edit, sep, useLogical, env=env)
    } else if ( op %in% c("(","{") ){
        edit <- parseCat(x[[2]], val,  edit, sep, useLogical, env=env)
    } else if ( op %in% c("%in%","==") ){
        cat <- eval(x[[3]], envir=env)
        if ( is.na(val) && op == "==" ) val <- FALSE
        if (is.logical(cat) && useLogical){
            if (length(cat) > 1){
              val <- NA
            } else {
              if (!cat) val <- !val              
            }
            var <- as.character(x[[2]])
        } else {
            var <- paste(x[[2]],cat,sep=sep)
        }
        edit[var] <- val
    } else if (op == "!=") {
        cat <- eval(x[[3]], envir=env)
        if (is.logical(cat) && useLogical){
          var <- as.character(x[[2]])
          if (!cat) val <- !val
        } else{
          var <- paste(x[[2]],cat,sep=sep)
        }
        edit[var] <- !val
    } else if (op == "!") {
        if (is.na(val)){
          val <- FALSE
        }
        edit <- parseCat(x[[2]],!val,  edit, sep, useLogical, env=env)
    } else if (op %in% c("&", "&&")){
        if (!isTRUE(val)){
            stop( "Operator '",op,"' is not allowed in 'then' clause.\n Edit '"
                , deparse(x)
                ,"' can be split into multiple edits")
        }
        edit <- parseCat(x[[2]],TRUE, edit, sep, useLogical, env=env)
        edit <- parseCat(x[[3]],TRUE, edit, sep, useLogical, env=env)
    } else if (op %in% c("||","|")){
        if (isTRUE(val)){
          stop( "Operator '",op,"' is not allowed in 'if' clause.\n Edit '"
                , deparse(x)
                ,"' can be split into multiple edits")
        }
        edit <- parseCat(x[[2]],FALSE, edit, sep, useLogical, env=env)
        edit <- parseCat(x[[3]],FALSE, edit, sep, useLogical, env=env)
    } else {
        stop("Operator '",op,"' not implemented")
    }
    edit
}

isCat <- function(e){
  if (length(e)==1) {
    return(is.symbol(e))
  }
  
  cmp <- deparse(e[[1]])
  return( cmp %in% c(CATCMP,"!"))  
}

