##' @include methods.R
NULL

##' Constructor for grid layout container
##'
##' The grid layout container uses matrix notation to position its
##' child components. This allows one to align widgets both
##' horizontally and vertically, as desired. There is some support for
##' matrix methods, such as \code{dim} and \code{[} to reference the
##' child objects.
##' @param homogeneous are cells all the same size
##' @param spacing between cell spacing
##' @inheritParams gwidget
##' @seealso \code{\link{gformlayout}} for a more convenient means to layout forms.
##' @export
##' @examples
##' \dontrun{
##' 
##' w <- gwindow("glayout example", visible=FALSE)
##' g <- gvbox(container=w)
##' lyt <- glayout(container=g)
##' lyt[1,1] <- "a label"
##' lyt[1,2] <- gedit("A widget", container=lyt)
##' lyt[2, 1:2] <- gcombobox(state.name, cont=lyt)
##' g1 <- ggroup(container=g)
##' addSpring(g1)
##' gbutton("values", container=g1, handler=function(h, ...) {
##'   print(sapply(lyt[,2], svalue))
##' })
##' visible(w) <- TRUE
##' 
##' }
glayout <- function(
                    homogeneous = FALSE, spacing = 10, container = NULL,      ... ,
                    toolkit=guiToolkit()){
  obj <- .glayout (toolkit,
                   homogeneous=homogeneous, spacing=spacing, container=container ,...
                   )


  check_return_class(obj, "GLayout")
  obj   
    
}

##' generic for toolkit dispatch
##'
##' @export
##' @rdname glayout
.glayout <- function(toolkit,
                     homogeneous = FALSE, spacing = 10, container = NULL,
                     ... )
           UseMethod( '.glayout' )

##' refer to children
##'
##' The \code{[} method for the grid layout allows one to reference
##' the child objects by index. The return value is non standard. It
##' may be the item, a list (if one dimensonaL) or an array. The list
##' format is convenient to refer to all the child objects in a
##' column.
##' @param x object
##' @param i row index
##' @param j column index
##' @param drop drop return type?
##' @export
##' @rdname glayout
##' @method [ GLayout
##' @S3method [ GLayout
"[.GLayout" <- function(x, i, j, ..., drop=TRUE) {
  getWithDefault(drop, TRUE)

  if(missing(i)) i <- seq_len(dim(x)[1])
  if(missing(j)) j <- seq_len(dim(x)[2])
  
  if(length(i) == 1 && length(j) == 1)
    return(x$get_items(i, j, ..., drop=drop))

  ## a matrix or list
  out <- sapply(j, function(col) lapply(i, function(row) x[row, col]))
  if(is.matrix(out))
    return(out[,,drop=drop])
  else if(is.list(out) && length(out) == 1 && drop)
    return(out[[1]])
  else
    return(out)
}



##' add child components to layout using matrix notation
##'
##' The matrix notation allows for spanning of multiple rows and or columns, but no holes.
##' The \code{...} argument is used to pass in values for expand, fill, anchor (see
##' the \code{add} method of \code{ggroup}) for their meaning).
##' @param value constructor for a widget using this object as the
##' parent container
##' @export
##' @usage \method{[}{GLayout} (x, i ,j, ...) <- value
##' @rdname glayout
##' @method [<- GLayout
##' @S3method [<- GLayout
"[<-.GLayout" <- function(x, i, j, ..., value) {

  theArgs <- list(...)

  ## get expand, anchor, fill
  expand <- getWithDefault(theArgs$expand, FALSE)
  if(!is.null(theArgs$align))
    theArgs$anchor <- theArgs$align
  fill <- getWithDefault(theArgs$fil, "x") # "", x, y or both
  anchor <- getWithDefault(theArgs$anchor, NULL)

  x$set_items(value, i, j, expand, fill, anchor)
  x
}
