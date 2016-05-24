##' @include methods.R
NULL

##' constructor for widget to display heirarchical data
##'
##' The \code{gtree} widget is used to present structured heirarchical
##' data. This data may be specified through a data frame with some
##' accompanying columns by which to split the data, or dynamically
##' through a function (\code{offspring}).
##'
##' In the former case, the data frame is split up by the columns
##' specified by INDICES. The first index is used to give the intial
##' branches, the second index the second, etc. The end leaves are the
##' data associated with a given path, with key given by that column
##' specified by \code{chosen.col}
##'
##' For the latter case, the "path" of the current node (the node and
##' its ancestors) is passed to the \code{offspring} function which
##' computes the next level in the heirarchy. This level is specified
##' through a data frame. This data frame has special columns. The
##' \code{chosen.col} specifies which column is used as the key in the
##' path, the \code{icon.col} (when given) points to a stock icon name
##' to decorate the column. Similarly, the \code{tooltip.columns}. The
##' fact that a row in the data frame has offspring is specified
##' through the \code{offspring.col} column, again specifed by index
##' or column name.
##' @param x Data frame. Optional, if given specify INDICES
##' value to split data into heirarchical data structure
##' @param INDICES Integers or column names, referring to columns of \code{x}. Used to form heirarchical structure. Order is important.
##' @param offspring function. A function passed values \code{path} and \code{data}, the latter from \code{offspring.data}. The path is the current position of the parent item using the named keys from the chosen column.
##' @param offspring.data Passed to second argument of \code{offspring} function. Used to parameterize a function call.
##' @param chosen.col integer or one of column names of data frame
##' returned by \code{offspring}. The chosen column gives the key and
##' value of the path.
##' @param offspring.col integer or column name. Points to column containing logical values indicating if a row has offspring.
##' @param icon.col integer of one of the column names of the data
##' frame. If provided (non-NULL), then this column should provide a
##' stock icon name to be placed in the row for the given data.
##' @param tooltip.col integer or one of the column names of the data frame. If provided (non-NULL), then the row for this item will have a tooltip given by the pointed to value.
##' @param multiple logical. Is multiple selection allowed?
##' @inheritParams gwidget
##' @export
##' @examples
##' ##################################################
##' ## This tree reads a list
##' offspring <- function(path=character(0), lst, ...) {
##'   if(length(path))
##'     obj <- lst[[path]]
##'     else
##'       obj <- lst
##'   nms <- names(obj)
##'   hasOffspring <- sapply(nms, function(i) {
##'     newobj <- obj[[i]]
##'     is.recursive(newobj) && !is.null(names(newobj))
##'     })
##'   data.frame(comps=nms, hasOffspring=hasOffspring, ## fred=nms,
##'              stringsAsFactors=FALSE)
##' }
##' l <- list(a="1", b= list(a="21", b="22", c=list(a="231")))
##' 
##' \dontrun{
##' w <- gwindow("Tree test")
##' t <- gtree(offspring=offspring, offspring.data=l, cont=w)
##' }
##' 
##' ##################################################
##' ## This tree looks at recursive objects
##' describe <- function(x) UseMethod("describe")
##' describe.default <- function(x) sprintf("An object with class %s", class(x)[1])
##' describe.integer <- function(x) sprintf("An integer with %s value%s", length(x),
##'    ifelse(length(x) > 1, "s", ""))
##' describe.numeric <- function(x) sprintf("A numeric with %s value%s", length(x),
##'    ifelse(length(x) > 1, "s", ""))
##' describe.factor <- function(x) sprint("A factor with %s level%s", length(levels(x)),
##'    ifelse(length(levels(x)) > 1, "s", ""))
##' 
##' offspring <- function(path, obj) {
##'   if(length(path) > 0)
##'     x <- obj[[path]]
##'   else
##'     x <- obj
##' 
##'   nms <- names(x)
##'   recursive <- sapply(x, function(i) {
##'     is.recursive(i) &&
##'     !is.null(attr(i, "names")) &&
##'     length(i) > 0
##'     })
##'   descr <- sapply(x, describe)
##'   
##'   data.frame(Variable=nms, offspring=recursive, Description=descr, stringsAsFactors=FALSE)
##' }
##' 
##' l <- lm(mpg ~ wt, mtcars)
##' \dontrun{
##' w <- gwindow("test")
##' gtree(offspring=offspring, offspring.data=l, cont=w)
##' }
##' 
gtree <- function(x=NULL, INDICES=NULL,
                  offspring = x, offspring.data = NULL,
                  chosen.col = 1, offspring.col=2, icon.col=NULL, tooltip.col=NULL,
                  multiple = FALSE,
                  handler = NULL, action = NULL, container = NULL, ... ,
                  toolkit=guiToolkit()){


  deprecated_args=list(
    "hasOffspring"="Offspring designation now done by specifying a column",
    "chosencol"="Now chosen.col",
    "icon.FUN"="Use icon.col and add to the data frame"
  )

  if(!is.null(x) && is.data.frame(x)) {
    
    ## kludgy means to place data.frame + INDICES into offspring framework
    index_to_name <- function(i, x) {
      if(is.null(i)) return(NULL)
      ifelse(is.numeric(i), names(x)[i], i)
    }
    name_to_index <- function(nm, x) {
      if(is.null(nm))  return(NULL)
      ifelse(is.character(nm), match(nm, names(x)), nm)
    }
    move_column <- function(x, from ,to) {
      stopifnot(is(x, "data.frame"))
      from <- name_to_index(from, x)
      to <- name_to_index(to, x)
      
      if(from > to) {
        x[to:from] <- x[c(from, to:(from-1))]
        names(x)[to:from] <- names(x)[c(from, to:(from-1))]
      } else if(from < to) {
        x[from:to] <- x[c((from + 1):to, from)]
        names(x)[from:to] <- names(x)[c((from + 1):to, from)]
      }
      x
    }
    subset_vals <- function(items, path, INDICES) {
      if(length(path) == 0)
        return(rep(TRUE, nrow(items)))
      
      ind <- INDICES[seq_along(path)]
      ind <- sapply(ind, name_to_index, x=items)
      nms <- names(items)[ind]
      f <- function(varname, value) items[[varname]] == value
      apply(mapply(f, nms, path),1, function(x) Reduce("&&", x))
    }
    
    ## make things names
    INDICES <- sapply(INDICES, index_to_name, x)
    chosen.col <- index_to_name(chosen.col, x)
    icon.col <- index_to_name(icon.col, x)
    tooltip.col <- index_to_name(tooltip.col, x)
    
    ## we rearrange items: indices, chosen.col, ..offspring, [icons],[tooltip], rest
    if(!is.null(tooltip.col)) {
      x <- move_column(x, tooltip.col, 1)
      tooltip.col <- 3 + !is.null(icon.col)
    }
    if(!is.null(icon.col)) {
      x <- move_column(x, icon.col, 1)
      icon.col <- 3
    }
    
    x$..offspring <- TRUE
    x <- move_column(x, "..offspring", 1)
    offspring.col <- 2

    x <- move_column(x, chosen.col, 1)
    for(i in rev(INDICES)) {
      x <- move_column(x, i, 1)
    }
    ## make numeric
    INDICES <- seq_along(INDICES)
    chosen.col <-  1L
    
    ## define offspring.data; offspring
    offspring.data <- list(x=x, no_indices=length(INDICES), extra_columns=length(unlist(list(icon.col, tooltip.col))))
    
    offspring <- function(path, data) {
      ## we have chosen_col, ..offspring, [icon],[tooltip], others
      items <- data$x; no_indices <- data$no_indices; extra_columns <- data$extra_columns
      ## restrict rows to match path
      items <- items[subset_vals(items, path, INDICES), ]

      ## restrict columns 
      n <- no_indices + 1
      this_one <- length(path) + 1
      not_these <- setdiff(seq_len(n), this_one)
      items <- items[-not_these]
      
      if(length(path) == no_indices) {
        items$..offspring <- FALSE
        items
      } else {
        ## reduce by unique values in first row
        items <- items[order(items[,1]), ]
        a <- items[,1]; n <- length(a)
        if(n > 1)
          items <- items[c(TRUE, a[1:(n-1)] != a[2:n]), ]

        ## clear out rows
        clear_out <- function(x) UseMethod("clear_out")
        clear_out.default <- function(x) ""
        clear_out.numeric <- function(x) NaN
        
        if(ncol(items) > 2 + extra_columns) {
          for(j in (3 + extra_columns):ncol(items)) {
            items[[j]] <- clear_out(items[[j]])
          }
        }
        names(items)[1] <- gettext("Key")
      }
      return(items)
    }
  }
  obj <- .gtree (toolkit,
                 offspring=offspring, offspring.data=offspring.data,
                 chosen.col=chosen.col, offspring.col=offspring.col, icon.col=icon.col, tooltip.col=tooltip.col,
                 multiple=multiple,
                 handler=handler, action=action, container=container ,...
                 )
  check_return_class(obj, "GTree")
  return(obj)
  
}

##' generic for toolkit dispatch
##'
##' @export
##' @rdname gtree
.gtree <-  function(toolkit,
                    offspring = NULL, offspring.data = NULL,
                    chosen.col = 1, offspring.col=2,
                    icon.col=NULL, tooltip.col=NULL, 
                    multiple = FALSE,
                    handler = NULL, action = NULL, container = NULL, ... )
  UseMethod( '.gtree' )


## methods


##' svalue method
##'
##' For a \code{GTree} object, svalue refers to the path specified as
##' a  vector of keys or (if \code{INDEX=TRUE}) by an integer vector
##' of offspring positions. The \code{drop} argument is used to
##' indicate if the terminus of the path is returned or the entire
##' path, defaults=TRUE. To get the data associated with a row, use the \code{[} method.
##' @param obj object
##' @param index index
##' @param drop do we return tip or path
##' @export
##' @rdname gtree
##' @method svalue GTree
##' @S3method svalue GTree
svalue.GTree <-  function(obj, index=FALSE, drop=TRUE, ...) NextMethod()

##' \code{svalue<-} method
##'
##' For a \code{GTree} object, svalue refers to the path specified as
##' a  vector of keys . For the assignment method, one assigns by
##' index. That is \code{svalue(tr, index=TRUE) <- svalue(tr,
##' index=TRUE)} should not change the state of the widget. (The
##' \code{index=TRUE} argument is the case for setting, but not
##' getting.)
##' @param value vector of indices
##' @export
##' @usage \method{svalue}{GTree} (obj, index=TRUE, ...) <- value
##' @rdname gtree
##' @method svalue<- GTree
##' @S3method svalue<- GTree
"svalue<-.GTree" <-  function(obj, index=TRUE,  ..., value) NextMethod()


##' extract method for gtree
##'
##' The \code{[} method is used to return the data associated with a
##' selected row. The \code{svalue} method returns the path or its
##' endpoint, the \code{[} method returns the row data associated with
##' the path.
##' @param i ignored
##' @param j ignored
##' @export
##' @rdname gtree
##' @method [ GTree
##' @S3method [ GTree
"[.GTree" <- function(x, i, j, ..., drop=FALSE) {
  if(isExtant(x))
    x$get_items(i, j, ..., drop=drop)
}

##' update method
##'
##' The update method for \code{GTree} recomputes the base nodes, then reopens the given node if still available
##' @param object object to update
##' @param ... passed to update method
##' @export
##' @rdname gtree
##' @method update GTree
##' @S3method update GTree
update.GTree <- function(object, ...) NextMethod()
