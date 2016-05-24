##' @include BasicInterface.R
NULL


##' A class for monitoring workspace changes
##'
##' A reference class to create a model that monitors the global
##' workspace. The class has method \code{update_state} and the
##' "getting" methods \code{get_by_class}, \code{get_by_function}
##' (filter), \code{get_changes}. Use with a \code{gtimer} instance to
##' keep up to date with changes to the workspace.
##' @aliases WSWatcherModel
##' @exportClass WSWatcherModel
##' @rdname S4-classes
##' @name WSWatcherModel-class
WSWatcherModel <-  setRefClass("WSWatcherModel",
                               contains="Observable",
                               fields=list(
                                 nms="character",
                                 digests="ANY",
                                 old_nms = "character",
                                 old_digests="ANY",
                                 changes = "list" # should signal when updated
                                 ),
                               methods=list(
                                 initialize=function( ...) {
                                   "Initialze state of cached objects"

                                   
                                   update_state() # initial
                                   old_nms <<- nms
                                   old_digests <<- digests
                                   callSuper(...)
                                 },
                                 update_state=function(...) {
                                   "update cache of names/digests, then notify observers if there are changes"
                                   nms <<- ls(envir=.GlobalEnv)
                                   ## avoid certain classes
                                   skip_these <- function(x) !(is(x, "DigestClass") || is(x, "envRefClass"))
                                   digests <<- sapply(Filter(skip_these, mget(nms, .GlobalEnv)), digest)
                                   if(any_changes()) {
                                     ## Update the "changes"
                                     is_changed <- function(i) digests[i] != old_digests[i]
                                     changes <<- list(removed=setdiff(old_nms, nms),
                                                      added=setdiff(nms, old_nms),
                                                      changed=Filter(is_changed, intersect(old_nms, nms))
                                                      )
                                     old_nms <<- nms
                                     old_digests <<- digests
                                     ## notify any observers if there are any changes
                                     notify_observers()
                                   }
                                 },
                                 any_changes=function(...) {
                                   "Report  if any changes"
                                   if(length(old_nms) == 0) {
                                     out <- TRUE
                                   } else  {
                                     out <- (length(old_digests) != length(digests)) || any(old_digests != digests)
                                   }
                                   out
                                 },
                                 ## get
                                 get_by_class = function(classes=character(0)) {
                                   "Return objects matching any of classes"
                                   if(length(classes) == 0)
                                     return(nms)
                                   f <- function(x) Reduce("||", sapply(classes, is, object=x))
                                   get_by_function(f)
                                 },
                                 get_by_function= function(f) {
                                   "Filter objects by function f"
                                   
                                   objs <- mget(nms, .GlobalEnv, ifnotfound=list(function(x) {}))
                                   Filter(f, objs)
                                 },
                                 filter_names=function(f) {
                                   "Filter the names by f"
                                   objs <- mget(nms, .GlobalEnv, ifnotfound=list(function(x) {}))
                                   objs[Filter(f, nms)]
                                 },
                                 get_changes=function() {
                                   "Return list of changes"
                                   changes
                                 },
                                 pop_changes= function() {
                                   "pop changes, reset"
                                   x <- get_changes()
                                   changes <<- list(removed=character(0),
                                                    added=character(0),
                                                    changed=character(0))
                                   x
                                 }
                                 ))
