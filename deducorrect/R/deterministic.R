
##-------------------------------------------------------------------------
# default symbols allowed to define correction rules.
.onLoad <- function(libname,pkgname){
   options(allowedSymbols = c(
      'if', 'else', 'is.na','is.finite',
      '==','<','<=','=','>=','>','!=','!', '%in%',
      'identical','sign','abs',
      '||', '|', '&&', '&', 
      '(','{','<-','=',
      '+', '-', '*', '^', '/', '%%', '%/%'
      )
   )
}

.onUnload <- function(libpath){
   options(allowedSymbols=NULL)
}


##-------------------------------------------------------------------------
# define correction rules.



#' Rules for deterministic correction
#'
#' @section Details:
#' Data editing processes are rarely completely governed by in-record consistency rules.
#' Many \emph{ad-hoc} rules are commonly used to impute empty or erroneous values. 
#' Such rules are often applied manually or hidden in source code. This
#' function, together with \code{\link{correctWithRules}} allows for easy definition and execution 
#' of simle deterministic replacement rules.
#'
#' These functions are ment to support very simple rules, such as \emph{if variable x is missing, then
#' set it to zero}. Such actions usually basically model-free corrections stemming from subject-matter knowledge.
#' Given the nature of such rules, the type of rules are by default limited to R-statements containing
#' conditionals (\code{if}-\code{else}), arithmetic and logical operators, and brackets and assignment operators.
#' see \code{getOption('allowedSymbols')} for a complete list.
#'
#' If you cannot execute your 'simple' corrections with just these functions, we strongly recommend to 
#' write a separate imputation or correction routine. However, it's a free world, so you may alter the list of allowed symbols
#' as you wish. 
#' 
#' @section Note:
#' \code{getVars} is overloaded from the \code{editrules} package.
#'
#' @param x \code{character} or \code{expression} vector. 
#' @param strict If \code{TRUE} an error is thrown if any forbidden symbol is used (see details).
#' @param allowed A \code{character} vector of allowed symbols
#' @param ... Currently unused.
#' @return \code{correctionRules} returns an object of class \code{correctionRules}
#' 
#'
#' @export
#' @seealso \code{\link{correctWithRules}}
correctionRules <- function(x, strict=TRUE, allowed=getOption('allowedSymbols'), ...){
   UseMethod('correctionRules')
}


#' @method correctionRules character
#' @param file If \code{file=TRUE}, \code{x} is treated as a filename from which the rules are read.
#' @rdname correctionRules
#' @export
correctionRules.character <- function(x, strict=TRUE, allowed=getOption('allowedSymbols'), file=TRUE, ...){
   if ( file ){ 
      x <- parse(file=x)
      return(correctionRules.expression(x,strict,allowed))
   } 
    
   
   i <- 0
   L <- sapply(x, function(y){
      i <<- i+1
      tryCatch(parse(text=y), 
            error = function(e){
               msg <- sprintf("\nCould not process rule [%d]:\n '%s'\nparser returned:\n%s\n",i,y,e$message)
               stop(msg)
            })
   }
   )
   if (strict){ 
      M <- checkRules(L,allowed=allowed)
      if ( any(M$error) ){ 
         printErrors(x,M)
         stop('Forbidden symbols found')
      }
   }
   structure(L,class='correctionRules')
}


#' @method correctionRules expression
#' @rdname correctionRules
#' @export
correctionRules.expression <- function(x,strict=TRUE, allowed=getOption('allowedSymbols'), ...){

   if (strict){ 
      M <- checkRules(x,allowed=allowed)
      if ( any(M$error) ){ 
         printErrors(x,M)
         stop("Forbidden symbols found")
      }
   }
   structure(x,class='correctionRules')
}

#' @method print correctionRules 
#' @export
#' @rdname correctionRules
print.correctionRules <- function(x,...){
   cat("Object of class 'correctionRules'")
   v <- as.character(x)
   v <- gsub("^","  ",v)
   v <- gsub("\n","\n  ",v)
   cat(sprintf("\n## %2d-------\n%s",1:length(v),v),'\n')
}

#' @method as.character correctionRules
#' @param oneliner Coerce to oneliner
#' @export
#' @rdname correctionRules
as.character.correctionRules <- function(x, oneliner=FALSE,...){
   # this seems to be the easiest way to retain formatting information (mvdl)
   v <- sapply(x,function(r) as.character(noquote(list(r))))
   if ( oneliner ){
      v <- gsub("\n"," ",v)      # remove line end
      v <- gsub("[ ]+"," ",v)    # collapse multiple spaces to one
   }
   names(v) <- NULL
   v
}


##-------------------------------------------------------------------------
# check if expression contains forbidden symbols (boolean)
checkSymbols <- function(x, b=TRUE, allowed=getOption('allowedSymbols'), ...){
   if (is.expression(x)){ 
      x <- x[[1]]
   }

   if (length(x) == 1 ) return(TRUE)
   if (  is.symbol(x[[1]]) && 
         !(as.character(x[[1]]) %in% allowed) ){
            return(FALSE)
   }
   for ( i in 1:length(x) ){b <- b & checkSymbols(x[[i]],b)}
   b
}

##-------------------------------------------------------------------------
# extract forbidden symbols from an expression
extractSymbols <- function(x, allowed, L=character(0), ...){
   if ( is.expression(x)){
      x <- x[[1]]
   }
   if ( length(x) == 1 ) return(NULL)
   if ( is.symbol(x[[1]]) && 
      !( as.character(x[[1]]) %in% allowed ) ){
      return(as.character(x[[1]]))
   }
   for ( i in 1:length(x) )  L <- c(L, extractSymbols(x[[i]], allowed, L))
   L
}

checkRules <- function(x, allowed, ...){
   M <- lapply(x,extractSymbols, allowed=allowed,...)
   list(error = sapply(M,function(m) length(m) > 0), symbols=M)
}


##-------------------------------------------------------------------------
# print rules and their errors. 
# x : list of rules, M result of checkRules
printErrors <- function(x, M){
   ix <- which(M$error)
   v <- as.character(x[ix])
   v <- gsub("^","  ",v)
   v <- gsub("\n","\n  ",v)
   S <- lapply(M$symbols[ix], paste, collapse=", ")
   cat('\nForbidden symbols found:')
   cat(sprintf('\n## ERR %2d ------\nForbidden symbols: %s\n%s',ix,S,v),'\n')
}





#' @method getVars correctionRules
#' @rdname correctionRules
#' @param E object of class \code{\link{correctionRules}}
#'
#' @return \code{getVars} returns a character vector of variable names.
#' @export
getVars.correctionRules <- function(E, ...){
   unique(do.call(c,lapply(E,getvrs)))
}

getvrs <- function(x, L=character(0), ...){
   if ( is.expression(x) ){
      x <- x[[1]]
   }  
   if ( length(x) == 1){
      if ( is.symbol(x) ) return(as.character(x))
      return(NULL)
   }

   for ( i in 2:length(x) ) L <- c(L, getvrs(x[[i]],L))
   unique(L)
}



#' Deterministic correction
#'
#' Apply simple replacement rules to a \code{data.frame}.
#'
#' @section Details:
#' This function applies the the \code{rules} one by one to \code{dat} and logs
#' their actions. Rules are excuted in order of occurrence in the \code{\link{correctionRules}}
#' so order may matter for the final result. Rules are applied to one record at the time, so 
#' the use of statistical funtions such as \code{mean} is useless, and forbidden by default.
#' See \code{\link{correctionRules}} for details on the type of rules that are possible.
#'
#' @param rules object of class \code{\link{correctionRules}} 
#' @param dat \code{data.frame}
#' @param strict If \code{TRUE}, an error is produced when the rules use variables other than in the \code{data.frame}.
#' @seealso \code{\link{correctionRules}}
#'
#' @return list with altered data (\code{$corrected}) and a list of alterations (\code{$corrections}).
#'
#' @example ../examples/correctWithRules.R
#'
#' @export
correctWithRules <- function(rules, dat, strict=TRUE){
   if (strict){
      vars <- getVars(rules)
      I <- vars %in% names(dat)
      if (!all(I)) stop(
         sprintf("Variables '%s' in rules do not occur in data",paste(vars[!I],sep=", "))
      )
   }

   out <- dat
   m <- nrow(dat)
   n <- length(rules)
   row <- numeric(0)
   variable <- character(0)
   old <- character(0)
   new <- character(0)
   how <- character(0)
   vars <- colnames(dat)
   tr <- as.character(rules,oneliner=TRUE)
   
   for ( i in 1:m ){
      for ( j in 1:n ){
         #d <- out[i,,drop=FALSE]
         d <- within(out[i,,drop=FALSE],eval(rules[[j]]))

         if ( !all(equal(d,out[i,])) ){
            rule <- tr[j]
            w <- which(!equal(d,out[i,]))
            row <- c(row, rep(i,length(w)))
            variable <- c(variable,vars[w])
            old <- c(old, format(unlist(out[i,w])) )
            new <- c(new, format(unlist(d[1,w])) )
            how <- c(how,rep(rule,length(w)))
            out[i,] <- d
         }
      }
   }
   list(
      corrected=out, 
      corrections=data.frame(
            row=row,
            variable=variable,
            old=old,
            new=new,
            how=how,
            stringsAsFactors=FALSE
         )
   )
}

##-------------------------------------------------------------------------'
# NA-robust pairwise comparison of data.frames
equal <- function(d,e){
   dNA <- is.na(d)
   eNA <- is.na(e)
   d == e & !(dNA != eNA) | (dNA & eNA) 
}



