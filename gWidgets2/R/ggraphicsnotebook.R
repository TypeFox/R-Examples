##' @include ggraphics.R
##' @include gdfnotebook.R
NULL

##' A notebook widget holding plot devices
##'
##' @param width width in pixels
##' @param height height in pixels
##' @param dpi screen resolution
##' @inheritParams gwidget
##' @export
ggraphicsnotebook <- function(
                              width = dpi * 6, height = dpi * 6, dpi = 75, container = NULL,      ... ,
                              toolkit=guiToolkit()){

  
  if(is.character(toolkit))
    toolkit <- guiToolkit(toolkit)

  obj <- .ggraphicsnotebook (toolkit,
                            width=width, height=height, dpi=dpi, container=container ,...
                            )
  
  check_return_class(obj, "GGraphicsNotebook")
  return(obj)
}


##' S3 generic whose methods are implemented in the toolkit packages
##'
##' @rdname ggraphicsnotebook
##' @export
.ggraphicsnotebook <- function(toolkit, width, height, dpi, container, ...) UseMethod(".ggraphicsnotebook")




## Default notebook
##' toolkit implementation
##'
##' @rdname ggraphicsnotebook
##' @method .ggraphicsnotebook default
##' @S3method .ggraphicsnotebook default
.ggraphicsnotebook.default <- function(toolkit, width, height, dpi,  container, ...) {
  GGraphicsNotebook$new(toolkit, width=width, height=height, dpi=dpi, container=container, ...)
}

## basic subclass
GGraphicsNotebook <- setRefClass("GGraphicsNotebook",
                             contains="GNotebookOfPages",
                                 fields=list(
                                   width="numeric",
                                   height="numeric",
                                   dpi="numeric",
                                   ctr="numeric",
                                   close_btn="ANY"
                                   ),
                             methods=list(
                               initialize=function(toolkit=NULL, width, height, dpi,
                                 container=NULL, ...) {

                                 initFields(width=width,
                                            height=height,
                                            dpi=dpi,
                                            ctr=0)
                                 make_ui(container)
                                 add_handler_changed(function(h,...) {
                                   ## need page.no here -- not svalue(h$obj) which may be behind the count
                                   i <- h$page.no
                                   if(1 <= i && i <= length(pages)) {
                                     page <- pages[[i]]
                                     visible(page) <- TRUE
                                   }
                                   close_btn$set_enabled(get_length() > 0)
                                 })
                                 callSuper(toolkit)
                               },
                               make_ui=function(container) {
                                 g <- ggroup(expand=TRUE, horizontal=FALSE, container=container)
                                 tb_container <- ggroup(cont=g, spacing=0)
                                 add_toolbar(tb_container)
                                 widget <<- gnotebook(container=g, expand=TRUE, fill=TRUE)
                                 block <<- g$block
                                 ## close button
                               },
                               add_toolbar=function(tb_container) {
                                 gbutton("new", container=tb_container, handler=function(h, ...) {
                                   add_page()
                                 })
                                 close_btn <<- gbutton("close", container=tb_container, handler=function(h,...) {
                                   remove_page(get_cur_page())
                                 })
                                 close_btn$set_enabled(FALSE)

                                 
                                 addSpring(tb_container)
                               },
                               get_index_from_page=function(page) {
                                 "get page index in the pages list"
                                 which(sapply(pages, function(i) identical(i, page)))
                               },
                               add_page=function(name="Plot") {
                                 ctr <<- ctr + 1
                                 page <- ggraphics(container=widget,
                                                   label=sprintf("Plot %s", ctr),
                                                   expand=TRUE,
                                                   width=width, height=height, dpi=dpi
                                                   ) ## XXX modify name
                                 pages <<- c(pages, page)
                                 nms <<- c(nms, name)
                                 set_cur_page(length(pages))
                               },
                               page_change_handler=function(page.no) {
                                 "Called when page is changed"
                               },
                               save_plot=function() {
                                 
                               }
                               ))
