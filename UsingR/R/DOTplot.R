## make DOTplot generic

##' Make a DOTplot of a data set
##'
##' @param x a list, vector, formula, or matrix
##' @param ... passed along
##' @return produces graphic as side effect
##' 
##' @rdname DOTplot
##' @export DOTplot
DOTplot <- function(x, ...) UseMethod("DOTplot")

## Basic usage, x is a vector or list of vectors (eg. dataframe)

##' Make a  dotplot of a data set
##'
##' @param x list, vector, formula, matrix
##' @param pch plot character
##' @param names labels
##' @param main main title
##' @param delta offset
##' @param ... passed to \code{points}
##'
##' @rdname DOTplot
##' @export
DOTplot.default <-
  function(x,                           # list, vector, formula, matrix
           pch = 16,                    # what plot character
           names = NULL,                # labels for x (n=1), y axis
           main = NULL,                 # title
           delta = 0,                   # offset
           ...                          # passed onto points()
           ) 
  {

    if(is.matrix(x)) {                  # make x a list
      x = as.data.frame(x)
    } else if (is.numeric(x)) {
      tmp = deparse(substitute(x))
      x = list(x)
      names(x) = tmp
    }

    if(is.list(x)) {
      ## get n, names
      n = length(x)
      if(is.null(names)) {
        if(is.null(names(x))) {
          names(x) = 1:n
        } 
        names = names(x)
      }
    } else {
      stop("x must be a vector, matrix, data frame, list or formula")
    }

    ## max.dots are max for each group, heights are starting points of each group
    if(n > 1) {
      tmp = sapply(x,table)
      if (is.matrix(tmp)) tmp = as.data.frame(tmp)
      max.dots = sapply(tmp,max)  # a vector for each variable
      heights = cumsum(max.dots) + 1:length(max.dots)
      heights = c(0,heights[-length(heights)])
    } else {
      max.dots = max(table(x[[1]]))
      heights = 0
    }
    
    ## create window
    plot.new()
    plot.window(xlim=range(x,na.rm=TRUE),
                ylim=c(0-delta,heights[length(heights)] + max.dots[length(max.dots)]),
                bty="L"
                )

    if(!is.null(main)) title(main)      # title
    axis(1)                             # x axis
    if(n > 1) {                         # y axis
      axis(2,at=heights,labels=names)
    } else {
      title(xlab=names[1])
    }


    for(i in 1:n) {
      DOTplt(x[[i]],heights[i],pch=pch,...)
    }


}                                     # end of default

##' Make a  dotplot of a data set
##'
##' @param formula formula
##' @param data data frame for variable lookup
##' @param ... passed to default method
##' @param subset logical condition to subset data with
##'
##' @rdname DOTplot
##' @export
DOTplot.formula <- function(formula, data = NULL, ..., subset)
{
    if(missing(formula) || (length(formula) != 3))
        stop("formula missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    DOTplot(split(mf[[response]], mf[-response]), ...)
}

## Make a simple dotplot
"DOTplt" <-
  function(x,y,pch,...
           ) {

    tbl = table(x)
    x.vals = as.numeric(names(tbl))

    points(rep(x.vals,tbl),
           y+unlist(sapply(tbl,function(x) 0:(x-1))),
           pch = pch,...)

  }
