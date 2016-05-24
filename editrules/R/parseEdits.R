#' Parse a character vector of edits
#'
#' This function wraps the native \code{parse} function in a \code{tryCatch}.
#' The function is \code{editrules} internal. It tries to give a meaningfull error message when
#' parsing fails for some reason.
#'
#' @param E \code{character}
#' @param type optional filter for type of edit, may be \code{"all"}, \code{"num"}
#' \code{"cat"} or \code{"mix"}
#' @return The edits in \code{E} parsed to R expressions.
#'
#' @keywords internal
parseEdits <- function(E, type=c("all", "num", "cat", "mix")){
     if (is.expression(E)){
       edits <- E
     } else {
     edits <- 
        tryCatch(parse(text=E), 
            error=function(e){
                stop(paste("Not all edits can be parsed, parser returned", e$message,sep="\n"))
            })
     }
     
     #TODO fix this: should work, but currently generates errors in the tests.
#      nms <- names(E)
#      if (!is.null(nms)){
#        names(edits) <- make.unique(nms, sep="")
#      }
     
     type <- match.arg(type)
     if (type=="all"){
       return(edits)
     } 
     return(edits[editTypes(edits) == type])
}

parseTree <- function(expr,prefix=NULL){
   if (length(expr) == 1){
      indent <- paste("[", prefix,"]", sep="", collapse="")
      cat(indent, expr,"\n")
   }
   else {
       for (i in 1:length(expr)){
          parseTree(expr[[i]], c(prefix,i)) 
       }
   }
}

# basic test for type of edit
# edts \code{expression}
# 
editTypes <- function(edts){
  ops <- sapply(edts, function(e){deparse(e[[1]])})
  
  type  <- ifelse(ops %in% NUMCMP, "num", "cat")
  # todo add check for "=="

  iff <- ops == "if"
  mix <- sapply(edts[iff], hasNum)
  type[iff] <- ifelse(mix, "mix", "cat")
  as.factor(type)
}



