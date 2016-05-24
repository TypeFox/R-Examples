## Functions to grab the corner of data similar to head or tail

#' @title WhichCorner
#' @description Function to build the right row selection depending on the desired corner.
#' @details Function to build the right row selection depending on the desired corner.  Helper function for getting the indexing for data.frame's, matrices
#' @author Jared P. Lander
#' @aliases WhichCorner
#' @param corner (character) which corner to display c("topleft", "bottomleft", "topright", "bottomright")
#' @param r (numeric) the number of rows to show
#' @param c (numeric) the number of columns to show
#' @param object The name of the object that is being subsetted
#' @return An expression that is evaluated to return the proper portion of the data
#' @examples 
#' \dontrun{
#' WhichCorner('topleft')
#' WhichCorner('bottomleft')
#' WhichCorner('topright')
#' WhichCorner('bottomright')
#' 
#' WhichCorner('topleft', r=6)
#' WhichCorner('bottomleft', r=6)
#' WhichCorner('topright', r=6)
#' WhichCorner('bottomright', r=6)
#' 
#' WhichCorner('topleft', c=7)
#' WhichCorner('bottomleft', c=7)
#' WhichCorner('topright', c=7)
#' WhichCorner('bottomright', c=7)
#' 
#' WhichCorner('topleft', r=8, c=3)
#' WhichCorner('bottomleft', r=8, c=3)
#' WhichCorner('topright', r=8, c=3)
#' WhichCorner('bottomright', r=8, c=3)
#' }
#' 
WhichCorner <- function(corner=c("topleft", 'bottomleft', 'topright', 'bottomright'), r=5L, c=5L, object="x")
{
    corner <- match.arg(corner)
    
    theCorners <- list(
                    topleft=sprintf("list(rows=1:%s, cols=1:%s)", r, c),
                    bottomleft=sprintf("list(rows=(nrow(%s)-%s+1):nrow(%s), cols=1:%s)", object, r, object, c),
                    topright=sprintf("list(rows=1:%s, cols=(ncol(%s)-%s+1):ncol(%s))", r, object, c, object),
                    bottomright=sprintf("list(rows=(nrow(%s)-%s+1):nrow(%s), cols=(ncol(%s)-%s+1):ncol(%s))", object, r, object, object, c, object)
                )
    
    return(parse(text=theCorners[[corner]]))
}
# eval(WhichCorner(corner="topleft", 5, 5, "testFrame"))
# eval(WhichCorner(corner="topright", 5, 5, "testFrame"))
# eval(WhichCorner(corner="bottomleft", 5, 5, "testFrame"))
# eval(WhichCorner(corner="bottomright", 5, 5, "testFrame"))


## S3 generic function for getting the corner of data



#' Grabs a corner of a data set
#' 
#' Display a corner section of a rectangular data set
#' 
#' Displays a
#' corner of a retangular data set such as a data.frame, martrix or table.  If showing
#' the right side or bottom, the order of the data is preserved.
#' 
#' The default method reverts to simply calling \code{\link{head}}
#' 
## @aliases corner corner.default corner.data.frame corner.matrix
#' @title corner
#' @description Display a corner section of a rectangular data set
#' @details corner of a retangular data set such as a data.frame, martrix or table.  If showing the right side or bottom, the order of the data is preserved.
#' @author Jared P. Lander
#' @rdname corner
#' @param x The data
#' @param r Number of rows to display
#' @param c Number of columns to show
#' @param corner Which corner to grab.  Posibble values are c("topleft", "bottomleft", "topright", "bottomright")
#' @param \dots Arguments passed on to other functions
#' @return ... The part of the data set that was requested.  The size depends on r and c and the position depends on corner.
#' @seealso \code{\link{head}} \code{\link{tail}} \code{\link{topleft}} \code{\link{topright}} \code{\link{bottomleft}} \code{\link{bottomright}} \code{\link{left}} \code{\link{right}}
#' @export corner
#' @keywords corner head tail display subsection view
#' @examples
#' 
#' data(diamonds)
#' head(diamonds)      # displays all columns
#' corner(diamonds)    # displays first 5 rows and only the first 5 columns
#' corner(diamonds, corner="bottomleft")       # displays the last 5 rows and the first 5 columns
#' corner(diamonds, corner="topright")       # displays the first 5 rows and the last 5 columns
#'
corner <- function(x, ...)
{
    UseMethod("corner")
}


## gets the corner for a data.frame
## @x (data.frame) the data
## @r (numeric) the number of rows to show
## @c (numeric) the number of columns to show
## @corner (character) which corner to return, c("topleft", "bottomleft", "topright", "bottomright")
#' @rdname corner
#' @export
#'
corner.data.frame <- function(x, r=5L, c=5L, corner="topleft", ...)
{
    r <- if(nrow(x) < r) nrow(x) else r
    c <- if(ncol(x) < c) ncol(x) else c
    
    seqs <- eval(WhichCorner(corner=corner, r=r, c=c, object="x"))

    return(x[seqs$rows, seqs$cols, drop=FALSE])
}


## gets the corner for a matrix
## @x (matrix) the data
## @r (numeric) the number of rows to show
## @c (numeric) the number of columns to show
## @corner (character) which corner to return, c("topleft", "bottomleft", "topright", "bottomright")
#' @rdname corner
#' @export
corner.matrix <- function(x, r=5L, c=5L, corner="topleft", ...)
{
    r <- if(nrow(x) < r) nrow(x) else r
    c <- if(ncol(x) < c) ncol(x) else c
    
    seqs <- eval(WhichCorner(corner=corner, r=r, c=c, object="x"))
                 
    return(x[seqs$rows, seqs$cols, drop=FALSE])
}

#' @rdname corner
#' @export
corner.table <- function(x, r=5L, c=5L, corner="topleft", ...)
{
    r <- if(nrow(x) < r) nrow(x) else r
    c <- if(ncol(x) < c) ncol(x) else c
    
    seqs <- eval(WhichCorner(corner=corner, r=r, c=c, object="x"))
                 
    return(x[seqs$rows, seqs$cols, drop=FALSE])
}

## gets the corner for default
## @x (data) the data
## @r (numeric) the number of rows to show
#' @rdname corner
#' @importFrom utils head
#' @export
corner.default <- function(x, r=5L, ...)
{
    head(x, n=r, ...)
}


#' Grabs the top left corner of a data set
#' 
#' Display the top left corner of a rectangular data set
#' 
#' Displays the top left corner of a retangular data set.
#'
#' This is a wrapper function for \code{\link{corner}}
#' 
#' @aliases topleft
#' @rdname topleft
#' @param x The data
#' @param r Number of rows to display
#' @param c Number of columns to show
#' @param \dots Arguments passed on to other functions
#' @return ... The top left corner of the data set that was requested.  The size depends on r and c.
#' @author Jared P. Lander www.jaredlander.com
#' @seealso \code{\link{head}} \code{\link{tail}} \code{\link{corner}} \code{\link{topright}} \code{\link{bottomleft}} \code{\link{bottomright}} \code{\link{left}} \code{\link{right}}
#' @export topleft
#' @keywords corner head tail display subsection view
#' @examples
#' 
#' data(diamonds)
#' head(diamonds)      # displays all columns
#' topleft(diamonds)    # displays first 5 rows and only the first 5 columns
#' 
topleft <- function(x, r=5L, c=5L, ...)
{
    corner(x, r=r, c=c, corner="topleft", ...)
}



#' Grabs the top right corner of a data set
#' 
#' Display the top right corner of a rectangular data set
#' 
#' Displays the top right corner of a retangular data set.
#'
#' This is a wrapper function for \code{\link{corner}}
#' 
#' @aliases topright
#' @rdname topright
#' @param x The data
#' @param r Number of rows to display
#' @param c Number of columns to show
#' @param \dots Arguments passed on to other functions
#' @return ... The top right corner of the data set that was requested.  The size depends on r and c.
#' @author Jared P. Lander www.jaredlander.com
#' @seealso \code{\link{head}} \code{\link{tail}} \code{\link{corner}} \code{\link{topleft}} \code{\link{bottomleft}} \code{\link{bottomright}} \code{\link{left}} \code{\link{right}}
#' @export topright
#' @keywords corner head tail display subsection view
#' @examples
#' 
#' data(diamonds)
#' head(diamonds)      # displays all columns
#' topright(diamonds)    # displays first 5 rows and only the last 5 columns
#' 
topright <- function(x, r=5L, c=5L, ...)
{
    corner(x, r=r, c=c, corner="topright", ...)
}



#' Grabs the bottom left corner of a data set
#' 
#' Display the bottom left corner of a rectangular data set
#' 
#' Displays the bottom left corner of a retangular data set.
#'
#' This is a wrapper function for \code{\link{corner}}
#' 
#' @aliases bottomleft
#' @rdname bottomleft
#' @param x The data
#' @param r Number of rows to display
#' @param c Number of columns to show
#' @param \dots Arguments passed on to other functions
#' @return ... The bottom left corner of the data set that was requested.  The size depends on r and c.
#' @author Jared P. Lander www.jaredlander.com
#' @seealso \code{\link{head}} \code{\link{tail}} \code{\link{corner}} \code{\link{topright}} \code{\link{topleft}} \code{\link{bottomright}} \code{\link{left}} \code{\link{right}}
#' @export bottomleft
#' @keywords corner head tail display subsection view
#' @examples
#' 
#' data(diamonds)
#' head(diamonds)      # displays all columns
#' bottomleft(diamonds)    # displays last 5 rows and only the first 5 columns
#' 
bottomleft <- function(x, r=5L, c=5L, ...)
{
    corner(x, r=r, c=c, corner="bottomleft", ...)
}



#' Grabs the bottom right corner of a data set
#' 
#' Display the bottom right corner of a rectangular data set
#' 
#' Displays the bottom right corner of a retangular data set.
#'
#' This is a wrapper function for \code{\link{corner}}
#' 
#' @aliases bottomright
#' @rdname bottomright
#' @param x The data
#' @param r Number of rows to display
#' @param c Number of columns to show
#' @param \dots Arguments passed on to other functions
#' @return ... The bottom right corner of the data set that was requested.  The size depends on r and c.
#' @author Jared P. Lander www.jaredlander.com
#' @seealso \code{\link{head}} \code{\link{tail}} \code{\link{corner}} \code{\link{topright}} \code{\link{bottomleft}} \code{\link{topleft}} \code{\link{left}} \code{\link{right}}
#' @export bottomright
#' @keywords corner head tail display subsection view
#' @examples
#' 
#' data(diamonds)
#' head(diamonds)      # displays all columns
#' bottomright(diamonds)    # displays last 5 rows and only the last 5 columns
#' 
bottomright <- function(x, r=5L, c=5L, ...)
{
    corner(x, r=r, c=c, corner="bottomright", ...)
}



#' Grabs the left side of a data set
#' 
#' Display the left side of a rectangular data set
#' 
#' Displays the left side of a retangular data set.
#'
#' This is a wrapper function for \code{\link{corner}}
#' 
#' @aliases left
#' @rdname left
#' @param x The data
#' @param c Number of columns to show
#' @param \dots Arguments passed on to other functions
#' @return ... The left side of the data set that was requested.  The size depends on c.
#' @author Jared P. Lander www.jaredlander.com
#' @seealso \code{\link{head}} \code{\link{tail}} \code{\link{corner}} \code{\link{topright}} \code{\link{bottomleft}} \code{\link{bottomright}} \code{\link{topleft}} \code{\link{right}}
#' @export left
#' @keywords corner head tail display subsection view
#' @examples
#' 
#' data(diamonds)
#' head(diamonds)      # displays all columns
#' left(diamonds)    # displays all rows and only the first 5 columns
#' 
left <- function(x, c=5L, ...)
{
    corner(x, r=nrow(x), c=c, corner="topleft", ...)
}


#' Grabs the right side of a data set
#' 
#' Display the right side of a rectangular data set
#' 
#' Displays the right side of a retangular data set.
#'
#' This is a wrapper function for \code{\link{corner}}
#' 
#' @aliases right
#' @rdname right
#' @param x The data
#' @param c Number of columns to show
#' @param \dots Arguments passed on to other functions
#' @return ... The left side of the data set that was requested.  The size depends on c.
#' @author Jared P. Lander www.jaredlander.com
#' @seealso \code{\link{head}} \code{\link{tail}} \code{\link{corner}} \code{\link{topright}} \code{\link{bottomleft}} \code{\link{bottomright}} \code{\link{topleft}} \code{\link{topleft}}
#' @export right
#' @keywords corner head tail display subsection view
#' @examples
#' 
#' data(diamonds)
#' head(diamonds)      # displays all columns
#' right(diamonds)    # displays all rows and only the last 5 columns
#' 
right <- function(x, c=5L, ...)
{
    corner(x, r=nrow(x), c=c, corner="topright", ...)
}
