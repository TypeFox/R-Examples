##' Checkboxgroup with table

GCheckboxGroupTable <- setRefClass("GCheckboxGroupTable",
                                   contains=c("TcltkWidget"),
                                   fields=list(
                                     df="ANY",
                                     checked="logical",
                                     tr="ANY",
                                     frame="ANY"
                                     ),
                                    
                                   methods= list(
                                     init_widget=function(parent, items,
                                       checked=rep(TRUE, nrow(df)), ...) {
                                       if(is.data.frame(items)) {
                                         df <<- items
                                       } else {
                                         df <<- as.data.frame(items, stringsAsFactors=FALSE)
                                       }
                                       checked <<- rep(checked, length.out=nrow(df))
                                       m <- df
                                       for(i in seq_len(ncol(m))) m[,i] <- as.character(m[,i])

                                       a <- populate_rectangular_treeview(parent, m)
                                       frame <<- a$frame
                                       widget <<- tr <<- a$tr

                                       tkconfigure(tr, selectmode="none")
                                       tkconfigure(tr, show="tree headings")
                                       tcl(tr, "column", "#0", minwidth=40, width=40, anchor="center") 
                                       
                                       ## add in checked off
                                       show_checked()

                                       ## put in handler to toggle
                                       tkbind(tr, "<Button-1>", function(W, x, y) {
                                         if(!.self$is_enabled()) return()

                                         row <- as.character(tkidentify(W, "row", x, y))
                                         children <- as.character(tcl(W, "children", ""))
                                         i <- match(row, children)
                                         .self$toggle_checked(i)
                                       })
                                       .self
                                     },
                                     do_icon = function(i) {
                                       icon <- ifelse(.self$checked[i], "::image::on", "::image::off")
                                       ind <- as.character(tcl(tr, "children", ""))[i]
                                       tcl(tr, "item", ind, image=icon)
                                     },
                                     show_checked=function() {
                                       children <- as.character(tcl(tr, "children", ""))
                                       for(i in seq_len(no_items())) {
                                         icon <- ifelse(checked[i], "::image::on", "::image::off")
                                         tcl(tr, "item", children[i], image=icon)
                                       }

                                     },
                                     toggle_checked = function(i) {
                                       checked[i] <<- !checked[i]
                                       do_icon(i)
                                     },
                                     get_value = function(index=FALSE) {
                                       if(index)
                                         return(which(checked))
                                       else
                                         return(get_items()[checked])
                                     },
                                     set_value = function(val, index=FALSE) {
                                       ## index, logical or by name
                                       if(is.logical(val)) {
                                         checked <<- rep(val, length=no_items())
                                         show_checked()
                                         return()
                                       }

                                       if(is.logical(index) && !index) {
                                         val <- match(val, get_items())
                                         if(length(val) == 1 && is.na(val))
                                           val <- integer(0)
                                       }
                                       tmp <- rep(FALSE, length=no_items())
                                       tmp[val] <- TRUE

                                       checked <<- tmp
                                       show_checked()
                                     },
                                     get_items = function(drop=TRUE) {
                                       if(drop)
                                         df[,1,drop=TRUE]
                                       else
                                         df
                                     },
                                     set_items = function(new_items) {
                                       ## clear tree, add in items
                                       ##  nms in case not a data frame
                                       if(!is.data.frame(new_items)) {
                                         new_items <- data.frame(new_items, stringsAsFactors=FALSE)
                                       }
                                       if(ncol(new_items) != ncol(df))
                                         stop("Wrong number of columns")
                                       
                                       df <<- new_items
                                       checked <<- rep(checked, length=nrow(df))
                                       m <- df
                                       for(i in 1:ncol(m))
                                         m[,i] <- as.character(m[,i])
                                       
                                       ## clear old
                                       all_ind <- as.character(tcl(tr, "children", ""))
                                       sapply(all_ind, function(i) tcl(tr, "detach", i))
                                       ## add values
                                       apply(m, 1, function(vals) {
                                         if(length(vals) == 1) vals <- paste("{", vals, "}", sep="")
                                         tcl(tr, "insert", "", "end", values=vals)
                                         })
                                       show_checked()
                                     },
                                     no_items = function() length(get_items())
                                   ))


##' Helpers

## Now I need the gWIdgets interface

## build widget based on table
setClass("gCheckboxgroupTabletcltk",
         contains="gComponentR5tcltk",
         prototype=prototype(new("gComponentR5tcltk"))
         )

##' Constructor.
##'
##' Need not be a method, as only called internally
##' Standard parameters, but we don't need horizontal argument. (We do need toolkit.)
.gcheckboxgrouptable <- function(toolkit, items, checked = FALSE,
                   horizontal=FALSE, use.table=TRUE,
                   handler = NULL, action = NULL, container = NULL, ...) {

            tt = getWidget(container)

            cbg_widget <- getRefClass("GCheckboxGroupTable")$new(parent=tt,
                                                                 items=items,
                                                                 checked=checked)

            obj <- new("gCheckboxgroupTabletcltk",
                       block=cbg_widget$frame,
                       widget=cbg_widget$get_widget(),
                       R5widget=cbg_widget,
                       toolkit=toolkit,
                       e = new.env())

            svalue(obj) <- checked
            
            
            
            ## add to container
            add(container,  obj,...)
  
            ## add handler
            if(!is.null(handler))
              tag(obj, "handler.id") <- addhandlerchanged(obj, handler, action)

            invisible(obj)
          }


## ### methods
## setMethod(".svalue",
##           signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxgroupTabletcltk"),
##           function(obj, toolkit, index=NULL, drop=NULL, ...) {
            
##             cbg_widget <- obj@R5widget
##             index <- getWithDefault(index, FALSE)
##             if(index) {
##               return(cbg_widget$get_value(index=TRUE))
##             } else {
##               val <- cbg_widget$get_value()
##               return(val)
##             }
          
            

##           })

## ## toggles state to be T or F
## setReplaceMethod(".svalue",
##                  signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxgroupTabletcltk"),
##                  function(obj, toolkit, index=NULL, ..., value) {

##                    cbg_widget <- obj@R5widget
##                    index <- getWithDefault(index, FALSE)

##                    cbg_widget$set_value(value, index=index)
##                    return(obj)
##                  })

## ## [ and [<- refer to the names -- not the TF values

## setMethod("[",
##           signature(x="gCheckboxgroupTabletcltk"),
##           function(x, i, j, ..., drop=TRUE) {
##             .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
##           })
## setMethod(".leftBracket",
##           signature(toolkit="guiWidgetsToolkittcltk",x="gCheckboxgroupTabletcltk"),
##           function(x, toolkit, i, j, ..., drop=TRUE) {
##             cbg_widget <- x@R5widget
##             items <- cbg_widget$get_items()
##             if(missing(i))
##               items
##             else
##               items[i]
        
##           })

## ## assigns names
## setReplaceMethod("[",
##                  signature(x="gCheckboxgroupTabletcltk"),
##                  function(x, i, j,..., value) {
##                    .leftBracket(x, x@toolkit, i, j, ...) <- value
##                    return(x)
##                  })

## setReplaceMethod(".leftBracket",
##           signature(toolkit="guiWidgetsToolkittcltk",x="gCheckboxgroupTabletcltk"),
##           function(x, toolkit, i, j, ..., value) {

##             cbg_widget <- x@R5widget
##             if(!missing(i)) {
##               items <- cbg_widget$get_items()
##               items[i] <- value
##               value <- items
##             }
##             cbg_widget$set_items(value)

##              return(x)
##           })


## setMethod(".length",
##           signature(toolkit="guiWidgetsToolkittcltk",x="gCheckboxgroupTabletcltk"),
##           function(x,toolkit) {

##             cbg_widget <- x@R5widget
##             cbg_widget$no_items()
##           })


## ## inherited enabled isn't workgin                
## setReplaceMethod(".enabled",
##                  signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxgroupTabletcltk"),
##                  function(obj, toolkit, ..., value) {

##                    cbg_widget <- obj@R5widget                   
##                    cbg_widget$set_enabled(value)
##                    return(obj)
                  
##                  })


## This handler code is common to gradio and gcheckboxgroup. Should abstract out into a superclass.
## IF we do that, we should also use CheckButton bit
setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxgroupTabletcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {

            cbg_widget <- obj@R5widget            
            user.data=list(obj=obj, handler=handler, action=action)
            id <- cbg_widget$add_handler("<ButtonRelease-1>",
                                        handler=function(user.data) {
                                          h <- user.data[c("obj", "action")]
                                          user.data$handler(h)
                                  },
                                        user.data=user.data)
            invisible(id)
            
          })
## clicked is changed
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gCheckboxgroupTabletcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            .addhandlerclicked(obj, toolkit, handler, action, ...)
          })

