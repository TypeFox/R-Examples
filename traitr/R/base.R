##  Copyright (C) 2010 John Verzani
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  A copy of the GNU General Public License is available at
##  http://www.r-project.org/Licenses/


##' @include helpers.R
roxygen()
## to quiet down some errors and notes in R CMD check:
require(proto)
assign(".super", NULL)

##' Base Trait to place common properties and methods
##'
##' @export
BaseTrait <- proto(
                    ## we add a class to our objects to mimic class behavior (not dispatch though)
                    .doc_class=paste(
                      desc("There is no class information from <code>class()</code> for <code>proto</code> objects",
                           "this propertys allows us, in a lightweight manner, to keep track of a class for them.")
                      ),
                    class="TraitR",
                   ## identify as a traitr proto object
                   traitr=TRUE,
                    ## Convenience method to add to the class.
                    ## Not documented
                    add_class = function(., newclass) .$class <- c(newclass, .$class),
                    ## some methods to simplify tasks
                    ## make a new child object and call its initialization method if present
                    .doc_new=paste(
                      desc("Creates clone of proto object and calls the <code>init</code> method, if defined."),
                      param("...", "passed to <code>proto</code> call.")
                      ),
                    new = function(., ...) {
                      obj <- .$proto(...)
                      obj$do_call("init")
                      obj
                    },
                   ## assign if not null
                   .doc_assign_if_null=paste(
                     desc("Assign value or method to proto object if not present as local slot"),
                     param("key","Key to assign to"),
                     param("value", "Value to assign")
                     ),
                   assign_if_null=function(., key, value) {
                     if(!.$has_local_slot(key)) {
                       assign(key, value, envir=.)
                     } 
                   },
                    ## append to a property list
                    ## if a list, then optional key will set key for lookup
                    .doc_append=paste(
                      desc("Append a value to a property list."),
                      param("name","Property name"),
                      param("value","value"),
                      param("key","If a list, this adds a key for lookup by name. Otherwise, appended to end of list.")
                      ),
                    append=function(., name, value, key) {
                      val <- get(name, envir=.)
                      if(is.list(val)) {
                        if(!missing(key))
                          val[[key]] <- value
                        else
                          val[[length(val)+1]] <- value
                      } else {
                        val <- c(val, value)
                      }
                      assign(name, val, envir=.)
                    },

                    ## is method to check if the proto object has a class
                    ## the class is something we set, not an R class
                    .doc_is=paste(
                      desc("A function to test if object has any of classes specified."),
                      param("class","A character vector containing class names. (Not real classes, but those",
                            "defined within this package.)",
                            "Function returns <code>TRUE</code> if any mathches exists.")
                      ),
                    is = function(., class=NULL) {
                      if(!is.null(class))
                        any(class %in% .$class)
                      else
                        TRUE
                    },
                   ## some mutatr like methods
                   ## has_slot <--> exists
                   .doc_has_slot=paste(
                     desc("Does slot (property) exist?"),
                     param("key","Name of property"),
                     param("inherits","Look into ancestry?")
                   ),
                   has_slot=function(., key) {
                     key <- as.character(key)
                     exists(key, envir=.) && !is.null(get(key, envir=.))
                   },
                   .doc_has_local_slot=paste(
                     desc("Same as has_slot, only looks only in present object, not ancestors")
                     ),
                   has_local_slot=function(.,key) {
                     exists(as.character(key), envir=., inherits=FALSE)
                   },
                   ## get slot
                   .doc_get_slot=paste(
                     desc("Get slot (property) by name. Can be used to get a method to call with different context."),
                     param("key", "slot name"),
                     returns("Returns the method or NULL")
                     ),
                   get_slot=function(., key) {
                     if(.$has_slot(key))
                       get(as.character(key), envir=.)
                     else
                       NULL
                   },
                   .doc_get_local_slot=paste(
                     desc("Get slot if in local object, not ancestors"),
                     param("key", "slot name")
                     ),                     
                   get_local_slot = function(., key) {
                     if(.$has_local_slot(key))
                       get(as.character(key), envir=.)
                     else
                       NULL
                   },
                   .doc_set_slot=paste(
                     desc("Set slot (property) value"),
                     param("key","Property name"),
                     param("value","Value to set in slot"),
                     param("initialize","If <code>TRUE</code> assign even if not already present.")
                     ),
                   set_slot=function(., key, value, initialize=TRUE) {
                     key <- as.character(key)
                     if(!initialize && !.$has_slot(key)) {
                       warning(sprintf("property %s not intialized. Call with initialize=TRUE to set", key))
                       return()
                     } 
                     assign(key, value, envir=.)
                   },
                   ## return all objects in the widget
                   .doc_list_objects=paste(
                     desc("List all objects in the widget."),
                     returns("Returns a list with components methods and properties")
                     ),
                   list_objects = function(., all.names=FALSE) {
                     s <- .
                     is_method <- function(lst, env) {
                       if(length(lst) == 0)
                         return(logical())
                       ind <- sapply(lst, function(i) {
                         is.function(get(i,envir=env))
                       })
                       ind
                     }
                     out <- ls(s, all.names=all.names)
                     ind <- is_method(out, s)
                     methods <- out[ind]
                     properties <- out[!ind]
                     while(is.proto(s <- s$parent.env())) {
                       if(s$has_slot("class")) {
                         out <- ls(s, all.names=all.names)
                         ind <- is_method(out, s)
                         methods <- unique(c(methods, out[ind]))
                         properties <- unique(c(properties, out[!ind]))
                       }
                     }
                     return(list(methods=methods,
                                 properties=setdiff(properties, c("class","traitr","slot"))))
                   },
                   
                   ## return non functions (properties)
                   ## will return objects or just names if return_names=TRUE
                   ## can pass class= value if desired
                   .doc_list_properties=paste(
                     desc("Return all properties (non-methods) for this proto object."),
                     param("class","If non-NULL, returns only from objects of this class.")
                     ),
                   list_properties=function(., all.names=FALSE) .$list_objects(all.names)$properties,
                   ## list methods
                   .doc_list_methods=paste(
                     desc("Method to list all possible methods for object")
                     ),
                   list_methods=function(., all.names=FALSE) .$list_objects(all.names)$methods,
                   ## find next method that is different or return NULL
                   .doc_next_method=paste(
                     desc("Find next method for this object with same name by recursing up through parent",
                          "environments. Checks equivalence of methods after stripping out environment info."),
                     param("meth_name", "Name of method to find next one of"),
                     returns("The method or <code>NULL</code> if unable to find one.")
                     ),
                   next_method = function(.,meth_name) {
                     e <- new.env()
                     get_meth <- function(.) {
                       f <- get(meth_name,envir=.)
                       environment(f) <- e
                       f
                     }
                     cur_meth <- get_meth(.)

                     find_meth <- function(s,meth_name) {
                       p <- s$.super
                       if(exists(meth_name, envir=p)) {
                         f1 <- get_meth(p)
                         if(digest(cur_meth) != digest(f1)) {
                           return(get(meth_name, p))
                         } else {
                           find_meth(p, meth_name)
                         }
                       } else {
                         return(NULL)
                       }
                     }
                     find_meth(., meth_name)
                   },
                   
                   ## call a method if it is present. Basically do.call with
                    ## a check that the function exists
                    .doc_do_call=paste(
                      desc("Function to call method if the method exists."),
                      param("fun","method name as character"),
                      param("lst","List of arguments. Default is empty list")
                      ),
                    do_call = function(., fun, lst=list()) {
                      fun <- as.character(fun)
                      if(.$has_slot(fun) && is.function(.$get_slot(fun))) {
                        do.call(.$get_slot(fun), c(., lst))
                      }
                    },
                    ## doc stuff
                    ## produce a list with doc strings for each object
                    create_doc_list = function(.) {
                      ## get properties/methods with .doc_ entries
                      f <- function(s) {
                        tmp <- ls(s, all.names=TRUE)
                        ind <- grepl("^\\.doc_",tmp)
                        tmp[ind]
                      }
                      l <- list(); l$properties <- list(); l$methods <- list()
                      s <- .
                      while(is.proto(s)) {
                        for(i in f(s)) {
                          nm <- gsub("^\\.doc_","",i)
                          obj <- get(nm,envir=s)
                          if(is.function(obj) || is.null(obj))
                            type <- "methods"
                          else
                            type <- "properties"
                          if(is.null(l[[i]]))
                            l[[type]][[i]] <- list(class = s$class,
                                                   doc = get(i, envir=s))
                        }
                        s <- s$parent.env()
                      }
                      l
                    },
                    .doc_show_help = paste(
                      desc("Create HTML web page listing documented methods and properties for",
                           "the  given <code>proto</code> object"),
                      returns("Opens browser to help page.")
                      ),
                    show_help = function(.) {
                      tmpfile <- tempfile()
                      l <- .$create_doc_list()
                      cat("<html><head><style type=text/css>", file=tmpfile)
                      cat(readLines(system.file("css","traitr.css",package="traitr")),
                          sep="\n", file=tmpfile, append=TRUE)
                      ## cat("add css here", file=tmpfile, append=TRUE)
                      cat("</style></head><body>",
                          sep="\n", file=tmpfile, append=TRUE)
                      for(type in c("methods","properties")) {
                        cat("<h2>List of",type,"for object</h2>",
                            file=tmpfile, append=TRUE)
                        for(i in names(l[[type]])) {
                          if(!grepl("^\\.doc\\_Class",i)) {
                            nm <- gsub("^\\.doc\\_", "", i)
                            cat("<div class='traitrdoc'>",
                                "<h1 class='traitr'>",nm,"(",l[[type]][[i]]$class[1],")","</h1>",
                                "<div class='traitrdoc'>",l[[type]][[i]]$doc,"</div>",
                                "</div>",
                                "<hr></hr>",
                                sep="\n", file=tmpfile, append=TRUE)
                          }
                        }
                      }
                      cat("</body></html>",
                          sep="\n", file=tmpfile, append=TRUE)
                      ## browse
                      browseURL(tmpfile)
                    }
                    )

