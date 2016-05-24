##' @include methods.R
##' @include BasicInterface.R
NULL

## Put here until we can figure out how to get past R CMD check

##' A widget for filtering a data frame
##'
##' This widget provides a simple means to subset, or filter, a data
##' frame. 
##' @param DF a data frame or \code{GDf} instance to look variables up within.
##' @param allow.edit logical. If \code{TRUE} a user may add new
##' variables to filter by. If FALSE, then one should specify the
##' variables a user can filter by to \code{initial.vars}.
##' @param initial.vars When given, this is a data frame whose first
##' column specifies the variables within \code{DF} to filter by and
##' whose second column indicates the type of filter desired. The
##' available types are \code{single} to select one from many, 
##' \code{multiple}, for multiple selection; and \code{range}, to
##' specify a from and to value.
##' @inheritParams gwidget
##' @return returns \code{GFilter} object
##' @export
##' @examples
##' \dontrun{
##' DF <- mtcars[, c("mpg", "cyl", "hp", "am", "wt")]
##' w <- gwindow("Example of gfilter", visible=FALSE)
##' pg <- ggroup(container=w)
##' df <- gtable(DF, container=pg)
##' a <- gfilter(df, initial.vars=data.frame(names(DF), names(DF),
##'                    c("single", "multiple", "range", "single", "range"),
##'                    stringsAsFactors=FALSE),
##'              allow.edit=TRUE,
##'              container=pg,
##'              handler=function(h,...) {
##'                visible(df) <- h$obj$get_value()
##'              }
##'              )
##' size(w) <- c(600, 600)
##' visible(w) <- TRUE
##' }
gfilter <- function(DF,
                    allow.edit=TRUE, initial.vars=NULL,
                    handler=NULL, action=NULL,
                    container=NULL,
                    ...,
                    toolkit=guiToolkit()){
  
  
  obj <- .gfilter (toolkit,
                   DF=DF,
                   allow.edit=allow.edit, initial.vars=initial.vars,
                   handler=handler, action=action,
                   container=container,
                   ...
                   )
  check_return_class(obj, "GFilter")
  obj
}




##' generic for toolkit dispatch
##'
##' @export
##' @rdname gfilter
.gfilter <-  function(toolkit,
                      DF,
                      allow.edit=TRUE, initial.vars=NULL,
                      handler=NULL, action=NULL,
                      container=NULL,
                    ... )
           UseMethod( '.gfilter' )


##' svalue method
##'
##' The \code{svalue} method for a filter object returns a logical
##' containing which rows are selected. There is no assignment method.
##' @inheritParams svalue
##' @export
##' @rdname gfilter
##' @method svalue GFilter
##' @S3method svalue GFilter
svalue.GFilter <- function(obj, index=NULL, drop=NULL, ...)   NextMethod()


##' The extraction method returns the child items of the filter, like a container object.
##'
##' @param x the \code{GFilter} object
##' @param i passed to \code{get_items}
##' @param j passed to \code{get_items}
##' @param ... dots argument
##' @export
##' @rdname gfilter
##' @method [ GFilter
##' @S3method [ GFilter
"[.GFilter" <- function(x, i, j, ..., drop=TRUE) {
  if(isExtant(x)) {
    val <- x$get_items(i, j, ..., drop=drop)
    if(!is.null(drop) && drop && length(val) == 1)
      val <- val[[1]]
    return(val)
  } else {
    return(NULL)
  }
}



## Default implementation

##' Default
##'
##' @export
##' @rdname gfilter
##' @method .gfilter default
##' @S3method .gfilter default
.gfilter.default <-  function(
                              toolkit=guiToolkit(),
                              DF,
                              allow.edit=TRUE, initial.vars=NULL,
                              handler=NULL, action=NULL,
                              container=NULL,
                              ... ) {

  obj <- GFilter$new(toolkit,  DF=DF,
                     allow_edit=allow.edit, initial_vars=initial.vars,
                     handler=handler, action=action,
                     container=container,...)

  return(obj)
}



GFilter <- setRefClass("GFilter",
                       contains="GDefaultWidget",
                       fields=list(
                         DF="ANY", 
                         initial_vars="ANY",
                         allow_edit="logical",
                         na_filter = "ANY",
                         container="ANY",
                         l="list",
                         types="ANY",
                         presets="ANY",
                         handler="ANY",
                         action="ANY",
                           new_item_handler="ANY"
                         ),
                       methods=list(
                         initialize=function(
                           toolkit, 
                           DF=NULL,
                           initial_vars=NULL,
                           allow_edit=is.null(initial_vars),
                           handler=NULL, action=NULL,
                             new_item_handler=NULL,
                           container=NULL,
                           ...) {
                           initFields(DF=DF,
                                      initial_vars=initial_vars,
                                      allow_edit=allow_edit,
                                      types  = c("single"="Select one level", "multiple" = "Select multiple levels", "range"="Select range"),
                                      presets = c("preset"="Use a head()/tail()/some() filter", "generic"="Use a generic filter"),
                                      l=list(),
                                      handler=handler,
                                      action=action,
                                      new_item_handler=new_item_handler,
                                      toolkit=toolkit
                                      )
                           connect_df()                              
                           init_ui(container, ...)
                           callSuper()
                         },
                         add_handler_changed=function(handler, action=NULL, ...) {
                           ## we use a different handler mechanism
                           ## here. Each items calls the invoke
                           ## change handler withc in turn uses the
                           ## handler and action properties
                           handler <<- handler
                           action <<- action
                         },
                         invoke_change_handler=function() {
                           "Some value was changed"
                           if(is.null(handler))
                             return()
                           
                           h <- list(obj=.self, action=action)
                           handler(h)
                         },
                         init_ui=function(container, ..., use.scrollwindow) {
                           block <<- ggroup(container=container, horizontal=FALSE, ..., use.scrollwindow=FALSE)
                           container <<- ggroup(container=block, expand=FALSE, horizontal=FALSE)
                           
                           if(allow_edit && !is.null(DF)) {
                             bg <- ggroup(container=block)
                             addSpring(bg)                                
                             btn_add.preset <- gbutton(gettext("Preset"), container=bg, handler=function(h,...) {
                               w <- gbasicdialog(gettext("Use a preset or generic filter"),
                                                 handler=function(h,...) {
                                                   nms <- names(DF)
                                                   var <- nms[1]

                                                   preset <- svalue(preset, index=TRUE)
                                                   add_item(var, names(presets)[preset], type=names(presets)[preset])
                                                 }, parent=h$obj)
                               #presets <- c("preset"="Use a head()/tail()/some() filter", 
                               #             "generic"="Use a generic filter")
                               lyt <- glayout(container=w, expand=TRUE, fill=TRUE)
                               lyt[1,1] <- gettext("Select:")
                               
                               ##disable 'generic' temporarily
                               lyt[1,2, expand=TRUE, fill="x"] <- (preset <- gradio(presets[1], selected=1L, 
                                                      container=lyt))
                               visible(w) <- TRUE
                             })
                             btn_add.preset$set_icon("add")
                             btn_add.item <- gbutton(gettext("Variable"), container=bg, handler=function(h,...) {
                               w <- gbasicdialog(gettext("Select a variable and selector type"),
                                                 handler=function(h,...) {
                                                   var <- svalue(varname)
                                                   if(!(var %in% names(DF))) {
                                                     message(sprintf("There is no variable %s", var))
                                                     return()
                                                   }
                                                   
                                                   type <- svalue(type, index=TRUE)
                                                   names(types)[type]
                                                   add_item(var, var, type=type)
                                                 }, parent=h$obj)
                               nms <- names(DF)
                               lyt <- glayout(container=w, expand=TRUE, fill=TRUE)
                               lyt[1,1] <- gettext("Variable:")

                               lyt[1,2, expand=TRUE, fill="x"] <- (varname <- gcombobox(nms, selected=1L, # have a selection as otherwise can have issue
                                                                 editable=length(nms) > 20, use_completion=length(nms) > 20,
                                                                 container=lyt, handler=function(h,...) {
                                 nm <- svalue(h$obj)
                                 if(! (nm %in% names(DF))) {
                                   message("Name is not a match")
                                   visible(type) <- FALSE
                                 } else {
                                   visible(type) <- TRUE
                                 }

                                 var = DF[[nm]]
                                 lt <- get_avail_types(nm)
                                 type[] <- lt$types
                                 svalue(type, index=TRUE) <- lt$ind
                                 ## print(lt)
                                 ## if(is.numeric(var)) {
                                 ##   type[] <- lt$types
                                 ##   print(list("update type to ", lt$types))
                                 ##   svalue(type, index=TRUE) <- lt$ind
                                 ## } else if(is.factor(var) || is.character(var)) {
                                 ##   type[] <- lt$types                           
                                 ##   svalue(type, index=TRUE) <- lt$ind
                                 ## } else if(is.logical(var)) {
                                 ##   type[] <- lt$types
                                 ##   svalue(type, index=TRUE) <- lt$ind
                                 ## } else {
                                 ##   svalue(type, index=TRUE) <- lt$ind
                                 ## }
                                 enabled(type) <- TRUE
                               }))

                               
##                               lyt[2,1, anchor=c(-1,1)] <- gettext("Edit by:")

                               
                               lyt[2,2, expand=TRUE, fill="x"] <- (type <- gradio(types, selected=2, container=lyt))
                               lt <- get_avail_types(nms[1])
                               type[] <- lt$types

                               visible(w) <- TRUE
                             })
                             btn_add.item$set_icon("add")
                           }
                           addSpring(block) ## push to top
                           
                           ## add initial
                           if(!is.null(initial_vars))
                             sapply(seq_len(nrow(initial_vars)), function(i) {
                               add_item(initial_vars[i,1], name=initial_vars[i,2], type=initial_vars[i,3])
                             })
                           invoke_change_handler()
                         },
                         get_avail_types=function(nm) {
                           var <- DF[[nm]]
                           if(is.numeric(var)) {
                             list(types=types, ind=1)
                           } else if(is.factor(var) || is.character(var)) {
                             list(types=types[-3], ind=1)      
                           } else if(is.logical(var)) {
                             list(types=types[-3], ind=1)
                           } else {
                             list(types=types, ind=1)
                           }
                         },
                         connect_df=function() {
                           "connect DF to filter"
                           if(is(DF, "GDf"))
                             addHandlerChanged(DF, function(h,...) .self$invoke_change_handler())
                         }, 
                         get_x=function(x) {
                           if(is.character(x) && length(x) == 1)
                             x <- DF[,x]
                           x
                         },
                         add_item=function(x, name=deparse(substitute(x)), 
                                           type=c("single", "multiple", "range", "preset", "generic")) {
                           if(missing(type)) 
                             if(is.numeric(get_x(x)))
                               type <- "range"
                             else
                               type <- "multiple"
                           tmp <- type
                           if(is.numeric(type))
                             type <- c("single", "multiple", "range", "preset", "generic")[type]
                           else
                             type <- match.arg(type)

                           ## dispatch on type
                           if(type == "single")
                             item <- RadioItem$new(x, name=name, parent=.self)
                           else if(type == "multiple")
                             item <- ChoiceItem$new(x, name=name, parent=.self)
                           else if(type == "range")
                             item <- RangeItem$new(x, name=name, parent=.self)
                           else if(type == "preset")
                             item <- PresetItem$new(x, name=name, parent=.self)
                           else if(type == "generic")
                             item <- GenericItem$new(x, name=name, parent=.self)

                           l <<- c(l, item)
                           item$make_ui(visible=TRUE)
                           if (is(new_item_handler, "function"))
                               new_item_handler()
                           invoke_change_handler()
                         },
                         remove_item=function(child) {
                           ## remove from GUI
                           delete(container, child$frame)
                           ## remove from list
                           ind <- get_index_in_list(l, child)
                           l[[ind]] <<- NULL
                           invoke_change_handler()
                         },
                         show=function(...) cat("Object for showing filter list\n"),
                         ##
                         ## gWidgets interface: svalue and [
                         ##
                         get_value = function(drop=NULL, ...) {
                           "Return logical of length nrow(df)"
                           if(length(l)) {
                             m <- sapply(l, function(i) i$get_value())
                             apply(m, 1, all)
                           } else {
                             rep(TRUE, nrow(DF))
                           }
                         },
                         get_items=function(i,j, ..., drop=NULL) {
                           "Allow [ to access items"
                           l[i]
                         },
                         get_visible=function() visible(block),
                         set_visible=function(value) visible(block) <<- value,
                         get_enabled=function() enabled(block),
                         set_enabled=function(value) enabled(block) <<- value,
                         set_size=function(value) block$set_size(value)
                         ))



## Filter items
BasicFilterItem <- setRefClass("BasicFilterItem",
                               contains="GDefaultWidget",
                               fields=list(
                                 x="ANY",
                                 name="character",
                                 includeNA="ANY",
                                 parent="ANY",
                                 frame="ANY",
                                 widgetc="ANY",
                                 widget="ANY",
                                 prevnext="ANY"
                                 ),
                               method=list(
                                 initialize=function(x="", name=x, parent=NULL, includeNA=TRUE, prevnext=list(),...) {
                                   initFields(x=x,
                                              name=name,
                                              parent=parent,
                                              includeNA=includeNA
                                              )
                                   callSuper(...)
                                 },
                                 show=function(...) cat("A filter item\n"),
                                 get_x=function() {
                                   "Get value of x from looking for a data frame"
                                   if(is.character(x) && length(x) == 1) {
                                     DF <- parent$DF
                                     DF[,x]
                                   } else {
                                     x
                                   }
                                 },
                                 make_ui=function(visible=TRUE) {
                                   parent_container <- parent$container # or something else here
##                                   frame <<- gexpandgroup(name, container=parent_container, horizontal=FALSE,
##                                                          expand=TRUE, anchor=c(-1,1))
                                   ##frame$set_visible(visible)
                                   frame <<- gframe(name,  horizontal=FALSE,
                                                    container=parent_container
                                                    )
                                                    ##,
                                                    ##expand=TRUE, fill="x",
                                                    ##anchor=c(-1,1))

                                   make_item_type(container=frame)
                                   f <- function(i) addHandlerChanged(i, handler=function(h,...) {
                                     .self$invoke_change_handler()
                                   })
                                   if(is.list(widget))
                                     sapply(widget, f)
                                   else
                                     f(widget)

                                   make_buttons(frame)
                                   enabled(includeNA) <<- any(is.na(get_x()))
                                 },
                                 ## need to subclass these
                                 make_item_type=function(container) {
                                   "Make editor for item"
                                 },
                                 make_buttons=function(frame) {
                                   g <- ggroup(container=frame, horizontal=TRUE)

                                   includeNA <<- gcheckbox("NA", checked=FALSE, cont=g)
                                   tooltip(includeNA) <<- "Check to include NA"
                                   enabled(includeNA) <<- FALSE
                                   addHandlerChanged(includeNA, function(...) {
                                     parent$invoke_change_handler()
                                   })
                                   
                                   addSpring(g) # right justify
                                   b_reset <- gbutton("Reset", container=g, handler=function(h,...) {
                                     initialize_item()
                                     .self$invoke_change_handler()
                                   })

                                    addSpace(g, 5) # right justify
                                   if(parent$allow_edit) {
                                     b_rm <- gbutton("Remove", container=g, handler=function(h,...) {
                                       parent$remove_item(.self)
                                     })
                                     b_rm$set_icon("remove")
                                     tooltip(b_rm) <- "Remove filter"
                                   }
                                   disableFilter <- gcheckbox("", checked=FALSE, container=g, handler=function(h,...) {
                                     if(svalue(disableFilter)){ 
                                       enabled(widgetc) <<- FALSE 
                                       enabled(includeNA) <<- FALSE 
                                       enabled(b_reset) <- FALSE 
                                     } else if(!svalue(disableFilter)){ 
                                       enabled(widgetc) <<- TRUE
                                       enabled(includeNA) <<- any(is.na(get_x()))
                                       enabled(b_reset) <- TRUE
                                     }
                                     .self$invoke_change_handler()
                                   })
                                   tooltip(disableFilter) <- "Check to temporarily disable filter"
                                 },
                                 ## need to subclass these
                                 make_item_type=function(container) {
                                   "Make selector for item"
                                 },
                                 initialize_item=function() {
                                   "Method to initialize the item values"
                                 },
                                 invoke_change_handler=function() {
                                   ## pass along
                                   parent$invoke_change_handler()
                                 },
                                 ## gWidgets methods for later use.
                                 get_value=function(...) {
                                   "Return logical of length x"
                                 },
                                 do_na = function() svalue(includeNA),
                                 ## pass off to frame
                                 get_visible=function(...) visible(frame),
                                 set_visible=function(value, ...) visible(frame) <<- value,
                                 get_enabled=function(...) enabled(frame),
                                 set_enabled=function(value, ...) enabled(frame) <<- value
                                 ))

RadioItem <- setRefClass("RadioItem",
                         contains="BasicFilterItem",
                         methods=list(
                           make_item_type=function(container) {
                             "Select one from many"
                             u_x <- sort(unique(get_x()))
                             lots_o_them <- length(u_x) > 10
                             widgetc <<- ggroup(container=container, horizontal = FALSE, expand=TRUE, fill="y")
                             widget <<- gcombobox(u_x, container=widgetc, anchor=c(-1,0), editable=lots_o_them, use_completion=lots_o_them)
                             if(is.numeric(u_x))
                               widget$coerce_with <<- as.numeric
                             
                             initialize_item()
                           },
                           initialize_item = function() {
                             svalue(widget, index=TRUE) <<- 1L
                           },
                           make_buttons=function(frame) {
                             g <- ggroup(container=frame, horizontal=TRUE)

                             includeNA <<- gcheckbox("NA", checked=FALSE, cont=g)
                             tooltip(includeNA) <<- "Check to include NA"
                             enabled(includeNA) <<- FALSE
                             addHandlerChanged(includeNA, function(...) {
                               parent$invoke_change_handler()
                             })
                             
                             addSpring(g) # right justify
                             prevnext <<- list()
                             prevnext[["b_prev"]] <<- gbutton("", container=g, handler=function(h,...) {
                               svalue(widget, index=TRUE) <<- (svalue(widget, index=TRUE) - 1L)
                               .self$invoke_change_handler()
                             })
                             prevnext[["b_prev"]]$set_icon("go-up")
                             tooltip(prevnext[["b_prev"]]) <<- "Choose the previous level"
                             enabled(prevnext[["b_prev"]]) <<- FALSE
                             prevnext[["b_next"]] <<- gbutton("", container=g, handler=function(h,...) {
                               svalue(widget, index=TRUE) <<- (svalue(widget, index=TRUE) + 1L)
                               .self$invoke_change_handler()
                             })
                             prevnext[["b_next"]]$set_icon("go-down")
                             tooltip(prevnext[["b_next"]]) <<- "Choose the next level"
                             b_reset <- gbutton("Reset", container=g, handler=function(h,...) {
                               initialize_item()
                               .self$invoke_change_handler()
                             })

                             addSpace(g, 5) # right justify
                             if(parent$allow_edit) {
                               b_rm <- gbutton("Remove", container=g, handler=function(h,...) {
                                 parent$remove_item(.self)
                               })
                               b_rm$set_icon("remove")
                               tooltip(b_rm) <- "Remove filter"
                             }
                             disableFilter <- gcheckbox("", checked=FALSE, container=g, handler=function(h,...) {
                               if(svalue(disableFilter)){ 
                                 enabled(widgetc) <<- FALSE 
                                 enabled(includeNA) <<- FALSE 
                                 enabled(b_reset) <- FALSE 
                                 enabled(prevnext[["b_next"]]) <<- FALSE
                                 enabled(prevnext[["b_prev"]]) <<- FALSE
                               } else if(!svalue(disableFilter)){ 
                                 enabled(widgetc) <<- TRUE
                                 enabled(includeNA) <<- any(is.na(get_x()))
                                 enabled(b_reset) <- TRUE
                                 validx <- svalue(widget, index=TRUE)
                                 if(is.na(validx)){ 
                                    enabled(prevnext[["b_prev"]]) <<- FALSE
                                    enabled(prevnext[["b_next"]]) <<- FALSE
                                 } else {
                                    if(validx>1) enabled(prevnext[["b_prev"]]) <<- TRUE
                                    if(validx<length(na.omit(unique(get_x())))) enabled(prevnext[["b_next"]]) <<- TRUE
                                 }
                               }
                               .self$invoke_change_handler()
                             })
                             tooltip(disableFilter) <- "Check to temporarily disable filter"
                           },
                           get_value=function(...) {
                             if(enabled(widgetc)){
                               val <- svalue(widget)
                               if(length(val) == 0) # might have no choices in widget (all NA), This helps..
                                 val <- NA
                               
                               validx <- svalue(widget, index=TRUE)
                               if(is.na(validx)){ 
                                    enabled(prevnext[["b_prev"]]) <<- FALSE
                                    enabled(prevnext[["b_next"]]) <<- FALSE
                               } else {
                                   if(validx==1) enabled(prevnext[["b_prev"]]) <<- FALSE else 
                                      enabled(prevnext[["b_prev"]]) <<- TRUE
                                   if(validx==length(na.omit(unique(get_x())))) enabled(prevnext[["b_next"]]) <<- FALSE else 
                                      enabled(prevnext[["b_next"]]) <<- TRUE
                               }

                               out <- get_x() == val
                               out[is.na(out)] <- do_na()
                               out
                             } else if(!enabled(widgetc)){
                               out <- rep(TRUE, length(get_x()))
                               out
                             }
                           }
                           ))

ChoiceItem <- setRefClass("ChoiceItem",
                          contains="BasicFilterItem",
                          fields=list(
                            "old_selection"="ANY",
                            "search_type"="ANY"
                            ),
                           methods=list(
                           make_item_type=function(container) {
                             "Select one from many"
                             u_x <- as.character(sort(unique(get_x(), na.rm=TRUE)))
                             use.table <- length(u_x) > 4 # XXX make 4 part of parent so it can be configure
                             widgetc <<- ggroup(container=container, horizontal = FALSE, expand=TRUE, fill="y")
                             vb <- gvbox(container=widgetc)
                             search_type <<-  list(ignore.case=TRUE, perl=FALSE, fixed=FALSE)
                             if(use.table) {
                               gp <- ggroup(cont=vb)
                               ed <- gedit("", initial.msg="Filter values by...", expand=TRUE, container=gp)
                               ed$set_icon("ed-search", "start")
                               ed$set_icon("ed-remove", "end")
                               ed$set_icon_handler(function(h,...) {
                                 svalue(ed) <- ""
                                 focus(ed) <- TRUE
                               }, where="end")
                               ed$widget$setIconActivatable("primary", FALSE)
                               
                               search_handler <- function(h,..., do_old=TRUE) {
                                 ## we keep track of old selection here
                                 ## that updates only when user changes selection, not when filter does
                                 
                                 cur_sel <- old_selection
                                 blockHandlers(widget)
                                 on.exit(unblockHandlers(widget))
                                 val <- svalue(ed)

                                 if(val == "") {
                                   widget[] <<- u_x
                                 } else {
                                   l <- c(list(pattern=val, x=u_x), search_type)
                                   new_vals = u_x[do.call(grepl, l)]
                                   if (length(new_vals)) {
                                     widget[] <<- new_vals
                                   } else {
                                     widget[] <<- character(0) 
                                     return()
                                   }
                                 }

                                 svalue(widget) <<- cur_sel
                                 ## XXX
#                                 if(do_old)
#                                   old_selection <<- cur_sel
                               }

                               b <- gbutton("opts", cont=gp)
                               cbs <- list(gcheckbox("Ignore case", checked=TRUE, handler=function(h,...) {
                                                     search_type[["ignore.case"]] <<- svalue(h$obj)
                                                     search_handler(do_old=FALSE)
                                                     }),
                                           gcheckbox("Regex", checked=TRUE, handler=function(h,...) {
                                             search_type[["fixed"]] <<- !svalue(h$obj)
                                             search_handler(do_old=FALSE)                                                     
                                           }),
                                           gcheckbox("Perl compatible", checked=FALSE, handler=function(h,...) {
                                             search_type[["perl"]] <<- svalue(h$obj)
                                             search_handler(do_old=FALSE)                                                     
                                           })
                                           )
                               
                               addPopupMenu(b, gmenu(cbs, popup=TRUE))

                               addHandlerKeystroke(ed, search_handler)
                               addHandlerChanged(ed, search_handler)
                             }
                             
                             widget <<- gcheckboxgroup(u_x, container=widgetc,
                                                       use.table=use.table,
                                                       expand=TRUE, fill=TRUE
                                                       )
                             old_selection <<- svalue(widget)
                             
                             addHandlerChanged(widget, function(h,...) {
                               ### XXX selection
                               ## have to be careful, as items may be narrowed
                               visible_items = widget[]
                               new <- svalue(h$obj)
                               old <- intersect(visible_items, old_selection)

                               
                               added <- setdiff(new, old)
                               removed <- setdiff(old, new)

                               ## This is sort of tricky, not sure it is correct
                               if(length(added) > 0) {
                                 old_selection <<- unique(c(old_selection, added))
                               }
                               if(length(removed) > 0) {
                                 old_selection <<- setdiff(old_selection, removed)
                               }
                             })
                             if(length(u_x) >= 4){
                                size(widget) <<- list(height= 4 * 25)
                             } else {
                                size(widget) <<- list(height= length(u_x) * 25)
                             }
                             if(is.numeric(u_x))
                               widget$coerce_with <<- as.numeric
                             
                             initialize_item()
                           },
                           initialize_item = function() {
                             svalue(widget, index=TRUE) <<- FALSE # all NOT selected (unchecked)
                           },
                           make_buttons=function(frame) {
                             g <- ggroup(container=frame, horizontal=TRUE)

                             includeNA <<- gcheckbox("NA", checked=FALSE, cont=g)
                             tooltip(includeNA) <<- "Check to include NA"
                             addHandlerChanged(includeNA, function(...) {
                               parent$invoke_change_handler()
                             })
                             enabled(includeNA) <<- any(is.na(get_x()))

                             addSpring(g) # right justify
                             b_invert <- gbutton("", cont=g, handler = function(h,...) {
                                svalue(widget, index=TRUE) <<- setdiff(1:length(widget[]), 
                                                                       which(widget[] %in% old_selection))
                               .self$invoke_change_handler()
                             })
                             tooltip(b_invert) <- 'Invert selection'
                             b_invert$set_icon("jump-to")
                             b_selall <- gbutton("Select all", container=g, handler=function(h,...) {
                               #initialize_item()
                               svalue(widget, index=TRUE) <<- TRUE
                               .self$invoke_change_handler()
                             })
                             tooltip(b_selall) <- 'Select all'
                             b_selall$set_icon("select-all")
                             b_clear <- gbutton("Clear", container=g, handler=function(h,...) {
                               ## XXX
                               svalue(widget) <<- FALSE
                               .self$invoke_change_handler()
                             })
                             tooltip(b_clear) <- 'Select none'
                             
                             addSpace(g, 5) # right justify
                             if(parent$allow_edit) {
                               b_rm <- gbutton("", container=g, handler=function(h,...) {
                                 parent$remove_item(.self)
                               })
                               b_rm$set_icon("remove")
                               tooltip(b_rm) <- "Remove filter"
                             }
                             disableFilter <- gcheckbox("", checked=FALSE, container=g, handler=function(h,...) {
                               if(svalue(disableFilter)){ 
                                 enabled(widgetc) <<- FALSE 
                                 enabled(includeNA) <<- FALSE 
                                 enabled(b_selall) <- FALSE 
                                 enabled(b_clear) <- FALSE 
                                 enabled(b_invert) <- FALSE 
                               } else if(!svalue(disableFilter)){ 
                                 enabled(widgetc) <<- TRUE
                                 enabled(includeNA) <<- any(is.na(get_x()))
                                 enabled(b_selall) <- TRUE
                                 enabled(b_clear) <- TRUE
                                 enabled(b_invert) <- TRUE 
                               }
                               .self$invoke_change_handler()
                             })
                             tooltip(disableFilter) <- "Check to temporarily disable filter"
                           },
                           get_value=function(...) {
                             ## out <- get_x() %in% svalue(widget)
                             if(enabled(widgetc)){
                               out <- get_x() %in% old_selection
                               na_vals <- is.na(get_x())
                               out[na_vals] <- do_na()
                               out
                             } else if(!enabled(widgetc)){
                               out <- rep(TRUE, length(get_x()))
                               out
                             }
                           }
                           ))

RangeItem <- setRefClass("RangeItem", 
                         contains="BasicFilterItem",
                         methods=list(
                           make_item_type=function(container) {
                             "a <= widget <= b"
                             widget <<- list()

                             g <- ggroup(container=container, expand=TRUE, fill="y")
                             widgetc <<- g
                             g1 <- ggroup(container=widgetc, horizontal=FALSE, expand=TRUE)
                             widget[[1]] <<- gedit("", container=g1, width=10)

                             glabel(gettext("to"), container=widgetc)
                             
                             g2 <- ggroup(container=widgetc, horizontal=FALSE, expand=TRUE)
                             widget[[2]] <<- gedit("", container=g2, width=10)
                             initialize_item()
                             
                           },
                           initialize_item=function() {
                             sapply(widget, function(i) {
                               svalue(i) <- ""
                               i[] <- sort(unique(get_x(), na.rm=TRUE))
                             })
                             widget[[1]]$set_value(min(get_x(), na.rm=TRUE))
                             widget[[2]]$set_value(max(get_x(), na.rm=TRUE))

                           },
                           get_value=function(...) {
                             if(enabled(widgetc)){
                               vals <- get_x()
                               a <- svalue(widget[[1]])
                               b <- svalue(widget[[2]])
                               asx <- function(x, a) UseMethod("asx")
                               asx.default <- function(x, a) as.numeric(a)
                               asx.integer <- function(x, a) as.integer(a)
                               asx.character <- function(x, a) as.character(a)
                               asx.POSIXct <- function(x, a) as.POSIXct(a)
                               asx.Date <- function(x, a) as.Date(a)
                               a <- asx(vals,a); b <- asx(vals, b)
                               no_a <- is.null(a) || is.na(a)
                               no_b <- is.null(b) || is.na(b)


                               if(no_a & no_b) {
                                 out <- rep(TRUE, length(vals))
                               } else if(no_a) {
                                 out <- (vals <= b)
                               } else if (no_b) {
                                 out <- (a <= vals)
                               } else {
                                 out <- (a <= vals & vals <= b)
                               }
                               out[is.na(out)] <- do_na()
                               out
                             } else if(!enabled(widgetc)){
                               out <- rep(TRUE, length(get_x()))
                               out
                             }
                             }
                           ))

seq_sane <- (function(){
  master_seq <- unlist(lapply(seq_len(.Machine$double.digits), function(x){
    by <- 10^x
    seq.int(by, 10*by-1, by)
  }))
  #master_seq <- c(seq.int(1, 4, 1), seq.int(5, 50, 5), master_seq[-seq_len(5)])
  #master_seq <- c(seq.int(5, 20, 5), master_seq[-seq_len(2)])
  function(y){
    #if(any(is.na(y), is.nan(y))) y <- 0
    if(y<=10) seq.int(1, y, 1) else 
      #if(all(y>=200, y<=200)) else
         c(master_seq[master_seq < y], y)
  }
})()
# lapply(c(9,10,11,49,50,51,99,100,101,999,1000,1001,9999,10000,10001,385744), seq_sane)
# lapply(c(99,523,4548,27304), seq_sane)
# seq_sane(93)

PresetItem <- setRefClass("PresetItem",
                         contains="BasicFilterItem",
                         methods=list(
                           make_item_type=function(container) {
                             "Filter using head()/tail()/some()"
                             hts <- c("head", "tail", "some")
                            u_x <- get_x()
                            n_x <- length(u_x) 
                            seqn_x <- seq_sane(n_x)
                            seqn_x.def <- if(any(seqn_x %in% 1000)) which(seqn_x %in% 1000) else 
                              which(seqn_x %in% n_x)
                             #if(is.numeric(u_x))
                              # widget$coerce_with <<- as.numeric
                            widget <<- list()
                            g <- ggroup(container=container, expand=TRUE, fill="y")
                            g1 <- ggroup(container=g, horizontal=TRUE, expand=TRUE, fill='x')
                            widgetc <<- g1
                            ##default to 1000 as in RStudio
                             #widget[[1]] <<- gslider(1, n_x, by=10, value=1001, container=g1, 
                             #                   anchor=c(-1,0))
                            widget[[1]] <<- gradio(hts, selected=1L, horizontal=FALSE, container=widgetc)
                            addSpring(widgetc)
                            g2 <- ggroup(container=widgetc, horizontal=FALSE, expand=TRUE, fill='y')
                            widget[[2]] <<- gcombobox(seqn_x, selected=seqn_x.def, editable=TRUE, 
                                                      use_completion=TRUE, coerce.with=function(x1){
                                                          ##FIXME allow for negative values, too
                                                          x1 <- as.numeric(x1[1])  ##take only 1st element
                                                          if(any(is.na(x1), is.nan(x1))) x1 <- 0
                                                          x1 <- abs(x1)
                                                          return(x1)
                                                        }, container=g2)
                                
                             initialize_item()
                           },
                           initialize_item = function() {
                             ##FIXME how is it possible to avoid redoing this computation
                             u_x <- get_x()
                             n_x <- length(u_x) 
                             seqn_x <- seq_sane(n_x)
                             seqn_x.def <- if(n_x < 1000) which(seqn_x %in% n_x) else which(seqn_x %in% 1000)
                             svalue(widget[[1]], index=TRUE) <<- 1L
                             svalue(widget[[2]], index=TRUE) <<- seqn_x.def
                           },
                           make_buttons=function(frame) {
                             g <- ggroup(container=frame, horizontal=TRUE)

                             includeNA <<- gcheckbox("NA", checked=FALSE, cont=g)
                             tooltip(includeNA) <<- "Check to include NA"
                             enabled(includeNA) <<- FALSE
                             addHandlerChanged(includeNA, function(...) {
                               parent$invoke_change_handler()
                             })
                             visible(includeNA) <<- FALSE
                             
                             addSpring(g) # right justify
                             b_reset <- gbutton("Reset", container=g, handler=function(h,...) {
                               initialize_item()
                               .self$invoke_change_handler()
                             })

                             addSpace(g, 5) # right justify
                             if(parent$allow_edit) {
                               b_rm <- gbutton("Remove", container=g, handler=function(h,...) {
                                 parent$remove_item(.self)
                               })
                               b_rm$set_icon("remove")
                               tooltip(b_rm) <- "Remove filter"
                             }
                             disableFilter <- gcheckbox("", checked=FALSE, container=g, handler=function(h,...) {
                               if(svalue(disableFilter)){ 
                                 #if(is.list(widget)){
                                 #  for(i in 1:length(widget)) widget[[i]] <<- FALSE
                                 #}
                                 #enabled(widget[[1]]) <<- FALSE 
                                 #enabled(widget[[2]]) <<- FALSE 
                                 enabled(widgetc) <<- FALSE 
                                 enabled(b_reset) <- FALSE 
                               } else if(!svalue(disableFilter)){ 
                                 enabled(widgetc) <<- TRUE 
                                 enabled(b_reset) <- TRUE 
                               }
                               .self$invoke_change_handler()
                             })
                             tooltip(disableFilter) <- "Check to temporarily disable filter"
                           },
                           get_value=function(...) {
                               meth <- svalue(widget[[1]])
                               val <- svalue(widget[[2]])
                               #if(any(is.na(val), is.nan(val))) val <- 0
                               n_x <- length(get_x())
                               #if(length(val) == 0) # might have no choices in widget (all NA), This helps..
                                # val <- NA
                               
                               #out <- get_x() == val
                             if(enabled(widgetc)){
                               if(meth=="head") out <- 1:n_x %in% seq_len(val) else
                                 if(meth=="tail") out <- 1:n_x %in% seq.int(to=n_x, length.out=val) else
                                 if(meth=="some") out <- 1:n_x %in% sort(sample(n_x, min(val, n_x)))
                               out[is.na(out)] <- do_na()
                               out
                             } else if(!enabled(widgetc)){
                               out <- rep(TRUE, n_x)
                               out
                             }
                           }
                         ))

##currently not functional
GenericItem <- setRefClass("GenericItem",
                          contains="BasicFilterItem",
                          methods=list(
                            make_item_type=function(container) {
                              "Filter using a user-supplied logical vector"
                              u_x <- get_x()
                              n_x <- length(u_x) 
                              widget <<- gedit("", container=container, width=10)
                              initialize_item()
                            },
                            initialize_item = function() {
                              widget$set_value(TRUE)
                            },
                            get_value=function(...) {
                              val <- svalue(widget)
                              n_x <- length(get_x())
                              
#                                 if(!is.logical(val)) {
#                                   message(sprintf("The expression doesn't evaluate to a logical vector"))
#                                   return()
#                                 } else if(length(val)!=n_x){
#                                   message(sprintf("The logical vector has different length than the number of rows of the data frame"))
#                                   return()
#                                 }
                             out <- val
                              out[is.na(out)] <- do_na()
                              out
                            }
                          ))



  
