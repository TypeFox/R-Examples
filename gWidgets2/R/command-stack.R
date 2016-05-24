##' @include misc.R
NULL


##'  helper function to bypass lack of cached value in method call
##'
##' @param meth method name
##' @param obj method of object's class
##' @return the method
##' @note use as do.call(call_meth, args)
call_meth <- function(meth, obj) {
  out <- eval(substitute(obj$x, list(x=meth)))
  out
}


## Command and CommandStack classes

##' Command class
##'
##' Class for commands. Has methods do, redo, undo
##' @exportClass Command
##' @aliases Command
##' @rdname S4-classes
##' @name Command-class
Command <- setRefClass("Command",
                       fields=list(
                         receiver="ANY",
                         meth="character",
                         params="list",
                         old_params="list"
                         ),
                       methods=list(
                         initialize=function(receiver="", meth="", ...) {
                           l <- list(...)
                           initFields(receiver=receiver, meth=meth, params=l, old_params=l)
                           callSuper()
                         },
                         execute=function(args, meth_name=meth) {
                           do.call(call_meth(meth_name, receiver), args)
                         },
                         do=function() {
                           out <- execute(params)
                           old_params$value <<- out
                         },
                         redo=function() execute(params),
                         undo=function() execute(old_params)
                         ))
                       
                       
## Sample subclass
## cmd <- setRefClass("OtherCommand", contains="Command",
##                    methods=list(
##                      undo=function() message("huh")
##                    ))$new("Fred", "meth_name", "value")
## cmd$undo()

##' Class for multple commands
##'
##' @exportClass CommandList
##' @aliases CommandList
##' @rdname S4-classes
##' @name CommandList-class
CommandList <- setRefClass("CommandList",
                           fields=list(
                             l="list"
                             ),
                           methods=list(
                             initialize=function(..., lst) {
                               if(missing(lst))
                                 lst <- list(...)
                               initFields(l=lst)
                               callSuper()
                             },
                             add=function(cmd) {
                               l[[length(l) + 1]] <<- cmd
                             },
                             do=function() {
                               lapply(l, function(i) i$do())
                             },
                             redo=function() {
                               lapply(l, function(i) i$redo())
                             },
                             undo=function() {
                               lapply(l, function(i) i$undo())
                             }
                             ))

##' Stack to hold commands
##' 
##' A list with ptr. delegates call of do or undo to appropriate command
##' @exportClass CommandStack
##' @aliases CommandStack
##' @rdname S4-classes
##' @name CommandStack-class
CommandStack <- setRefClass("CommandStack",
                            fields=list(
                              l="list",
                              ptr="integer"
                              ),
                            methods=list(
                              initialize=function() {
                                initFields(l=list(), ptr=0L)
                                callSuper()
                              },
                              do=function() {
                                if(!can_do()) return()
                                cmd <- l[[ptr]]
                                ptr <<- ptr + 1L
                                cmd$do()
                              },
                              undo=function() {
                                if(!can_undo()) return()
                                cmd <- l[[ptr-1]]
                                ptr <<- ptr - 1L
                                cmd$undo()
                              },
                              redo=function() {
                                if(!can_do()) return()
                                cmd <- l[[ptr]]
                                ptr <<- ptr + 1L
                                cmd$redo()
                              },
                              can_do=function() {
                                ptr > 0 && ptr <= length(l)
                              },
                              can_undo=function() {
                                ptr > 1
                              },
                              add=function(cmd, call=TRUE) {
                                if(ptr <= 1) {
                                  l <<- list(cmd)
                                  ptr <<- 1L
                                } else {
                                  l <<- l[1:(ptr-1)]
                                  l[[length(l) + 1]] <<- cmd
                                }
                                if(call)
                                  do()
                              },
                              clear=function(cmd) {
                                l <<- list(); ptr <<- 0L
                              }
                              ))
                            

