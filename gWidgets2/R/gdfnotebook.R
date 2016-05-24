##' @include methods.R
##' @include BasicInterface.R
NULL


##' A notebook container for many \code{gdf} instances
##'
##' @export
##' @param items data frame for initial page, when given
##' @param container parent container
##' @param ... passed to \code{add} method of parent container
##' @param toolkit toolkit
gdfnotebook <- function(
                        items = NULL, container = NULL, ... ,
                        toolkit=guiToolkit()){

  if(is.character(toolkit))
    toolkit <- guiToolkit(toolkit)


  

  obj <- .gdfnotebook (toolkit,
                       items=items, container=container ,...
                       )
  
  check_return_class(obj, "GDfNotebook")
  return(obj)
}



##' S3 generic whose methods are implemented in the toolkit packages
##'
##' @rdname gdfnotebook
##' @export
.gdfnotebook <- function(toolkit, items,  container, ...) UseMethod(".gdfnotebook")

## Default


##' Toolkit constructor
##'
##' @inheritParams gdfnotebook
##' @export
##' @rdname gWidgets2-undocumented
##' @method .gdfnotebook default
##' @S3method .gdfnotebook default
.gdfnotebook.default <-  function(toolkit,
                                                 items,
                                                 container = NULL, ... ) {
  GDfNotebook$new(toolkit,
                  items, container = container, ...)
}

## ## class inherited here, and by ggrpahics..
 GNotebookOfPages <- setRefClass("GNotebookOfPages",
                                 contains="GDefaultWidget",
                                 fields=list(
                                   pages="list",
                                   nms="character"
                                   ),
                                 methods=list(
                                   initialize=function(toolkit=NULL,
                                     container=NULL, ...) {
                                  
                                  
                                  
                                     initFields(pages=list(),
                                                nms=character(0),
                                                default_expand=TRUE,
                                                default_fill=TRUE,
                                                toolkit=toolkit
                                                )

                                     callSuper(...)
                                   },
                                   make_ui=function(container) {
                                     g <- ggroup(expand=TRUE, horizontal=FALSE)
                                     tb_container <- ggroup(cont=g)
                                     add_toolbar(tb_container)
                                     widget <<- gnotebook(container=g, expand=TRUE, fill=TRUE)
                                     block <<- g$block
                                   },
                                   page_change_handler=function(page.no) {
                                     "Called when page is changed"
                                   },
                                   add_toolbar=function(tb_container) {
 ### XXX("Subclass")
                                   },
                                   get_index_from_page=function(page) {
                                     "get page index in the pages list"
                                     which(sapply(pages, function(i) identical(i, page)))
                                   },
                                   add_page=function(...) {
 ### XXX("sublass")
                                   },
                                   remove_page=function(i) {
                                     if(!is.numeric(i))
                                       i <- get_index_from_page(i)

                                     pages[[i]] <<- NULL
                                     widget$remove_page_by_index(i) ## remove from notebook
                                     if(length(nms) >= i)
                                       nms <<- nms[-i]
                                   },
                                   set_cur_page=function(i) {
                                     "Set current page to page i, a number or a page reference"
                                     if(!is.numeric(i))
                                       i <- get_index_from_page(i)
                                     widget$set_index(i)
                                   },
                                   get_cur_page=function() {
                                     cur_page <- svalue(widget)
                                     get_page(cur_page)
                                   },
                                   get_page=function(i) {
                                     pages[[i]]
                                   },
                                   ## These are passthroughs
                                   get_value=function( ...) {
                                     widget$get_value(...)
                                   },
                                   set_value=function(value, ...) {
                                     widget$set_value(value, ...)
                                   },
                                   get_index = function(...) {
                                     widget$get_index(...)
                                   },
                                   set_index = function(value,...) {
                                     widget$set_index(value, ...)
                                   },
                                   get_items = function(...) {
                                     widget$get_items(...)
                                   },
                                   set_items = function(value, ...) {
                                     widget$set_items(value, ...)
                                   },
                                   get_length=function(...) {
                                     widget$get_length()
                                   },
                                   get_names=function(...) {
                                     widget$get_names()
                                   },
                                   set_names = function(value, ...) {
                                     widget$set_names(value, ...)
                                   },
                                   ## delegate
                                   add_handler_changed=function(handler, action=NULL, ...) {
                                     widget$add_handler_changed(handler, action, ...)
                                   }
                                   ))
        
GDfNotebook <- setRefClass("GDfNotebook",
                           contains="GNotebookOfPages",
                            methods=list(
                              initialize=function(toolkit=NULL,
                                container=NULL, ...) {
                                ## put into subclass, otherwise we get an error
                                make_ui(container)
                           
                                ## set current page when page is changed
                                ## addHandlerChanged(widget, handler=function(h,...) {
                                ##   ##widget$add_handler_changed(handler=function(h,...) {
                                ##   set_cur_page(h$page.no)
                                ##   page_change_handler(h$page.no)
                                ## })

                                callSuper(toolkit)
                           
                              },
                              make_ui=function(container) {
                                g <- ggroup(expand=TRUE, horizontal=FALSE, container=container)
                                tb_container <- ggroup(cont=g, spacing=0)
                                add_toolbar(tb_container)
                                widget <<- gnotebook(container=g, expand=TRUE, fill=TRUE)
                                block <<- g$block

                           
                              },
                              add_toolbar=function(tb_container) {
                                gbutton("new", container=tb_container, handler=function(h, ...) {
                                  blank_df <- data.frame(lapply(1:10, function(i) rep("", 100)), stringsAsFactors=FALSE)
                                  names(blank_df) <- sprintf("X%s", 1:10)
                                  nm <- ginput("Variable name:", text="x", title="Name of new page?")
                                  if(nchar(nm))
                                    add_page(blank_df, make.names(nm))
                                })
                                gbutton("open", container=tb_container, handler=function(h,...) {
                                  ## present data frames in a list
                                  cur_dfs <- Filter(function(x) is.data.frame(get(x, .GlobalEnv)), ls(.GlobalEnv))
                                  if(length(cur_dfs) == 0) {
                                    galert(gettext("No data frames to choose from"), parent=block)
                                    return()
                                  } else if (length(cur_dfs) == 1) {
                                    add_page(get(cur_dfs, .GlobalEnv), cur_dfs)
                                  }
                                  if(length(cur_dfs) >= 2) {
                                    w <- gbasicdialog(gettext("Select a data frame to edit"), parent=block,
                                                      handler=function(h,...) {
                                                        if(length(val <- svalue(tbl))) {
                                                          add_page(get(val, .GlobalEnv), val)
                                                        }
                                                      })
                                    tbl <- gtable(cur_dfs, cont=w)
                                    size(tbl) <- c(300, 300)
                                    visible(w, set=TRUE)
                                  }
                                })
                                gbutton("close", container=tb_container, handler=function(h,...) {
                                  df <- get_cur_page()
                                  ## undo only in RGtk2
                                  if(gtoolkit() == "RGtk2" && df$can_undo()) {
                                    if(!gconfirm(gettext("Really close? There are unsaved changes"), parent=block))
                                      return()
                                  }
                                  remove_page(df)
                                })
                                gbutton("save", container=tb_container, handler=function(h,...) {
                                  save_DF(.GlobalEnv)
                                })
                              },
                              add_page=function(new_df, name=deparse(substitute(new_df))) {
                                page <- gdf(new_df, container=widget, label=name, expand=TRUE)
                                pages <<- c(pages, page)
                                nms <<- c(nms, name)
                                set_cur_page(length(pages))
                              },
                              undo = function(...) {
                                get_cur_page()$undo(...)
                              },
                              redo = function(...) {
                                get_cur_page()$undo(...)
                              },
                              save_DF=function(where) {
                                df <- get_cur_page()
                                nm <- nms[get_index_from_page(df)]
                                if(exists(nm, .GlobalEnv))
                                  if(!gconfirm(c("Variable exists", "really overwrite?")))
                                    return()
                           
                                df$save_data(nm, where)
                              }
                              ))

