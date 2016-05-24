# Author: Babak Naimi, naimi.b@gmail.com
# Date :  April 2015
# Version 1.0
# Licence GPL v3

if (!isGeneric(".testMethod")) {
  setGeneric(".testMethod", function(x, ...)
    standardGeneric(".testMethod"))
}  


setMethod('.testMethod', signature(x='.methodTemplate'), 
          function(x,template,arguments,outputs,test.args,...) {
            fo <- formals(template)
            xfo <- names(formals(x@Function))
            t1 <- all(unlist(lapply(names(fo),function(x) x %in% xfo)))
            t2 <- all(unlist(lapply(names(x@arguments),function(n) n %in% names(arguments))))
            if (t2) t3 <- all(unlist(lapply(names(x@arguments),function(n) x@arguments[n] == arguments[n])))
            else t3 <- FALSE
            
            if (!is.null(x@user.arguments)) {
              if (!is.null(x@user.argument.values)) {
                test.args <- c(test.args,x@user.argument.values)
                o <- try(do.call(x@Function,test.args),TRUE)
                t4 <- !inherits(o, "try-error")
                if (t4) t5 <- class(o) == class(outputs)
                else t5 <- FALSE
              } else t4 <- t5 <- TRUE
            } else {
              o <- try(do.call(x@Function,test.args),TRUE)
              t4 <- !inherits(o, "try-error")
              if (t4) t5 <- class(o) == class(outputs)
              else t5 <- TRUE
            }
            
            if (t4) {
              t5 <- class(o) == class(outputs)
            } else t5 <- TRUE
            
            if (!t1) print('Error: the function should have the same frame as the template!\n')
            if (!t2) print('Error: the reserved argument names should be used for the function; or new arguments by the used should be explicitely defined\n')
            if (!t3) print('Error: the type of function arguments does not match with the template\n')
            if (!t4) print('Error: function raised an error in a test call, check the template for the expected input and output!\n')
            if (!t5) print('Error: output class does not match with the template\n')
            if (t1 & t2 & t3 & t4 & t5) return(TRUE)
            else return(FALSE)
          }
          
          )
#--------

.newMethod <- function(name,f,Help=NULL) {
  n <- new('.methodTemplate',name=name[1])
  if (length(name) > 1) n@aliases <- name[2:length(name)]
  n@Help <- Help
  a <- formals(f)
  v <- lapply(a,function(x) class(eval(x)))
  n@arguments <- unlist(v)
  n@user.argument.values <- as.list(a)
  for (i in 1:length(a)) formals(f)[[i]] <- quote(expr =)
  n@Function <- f
  n
}
#---------
.templateMatch <- function(f,template) {
  a <- names(formals(f))
  b <- names(formals(template))
  
  if ('...' %in% b) {
    b <- b[b != '...']
    if (length(b) > 0) {
      if (!all(b %in% a)) stop('arguments in template do not match with the function...')
      a <- a[!a%in%b]
    }
    bo <- as.character(body(f))
    bo <- bo[2:length(bo)]
    ad <- c('dot <- list(...)')
    ad <- c(ad,paste(a,paste('<- dot',paste(paste('[["',a,sep=''),'"]]',sep=''),sep='')))
    bo <- c(ad,bo)
    for(i in 1:length(bo)) body(f)[i+1] <- parse(text=eval(as.expression(bo[i])))
    formals(f) <- formals(template)
  } else {
    if (!a %in% b) stop('arguments in template do not match with the function...')
  }
  f
}

#------

