
## update elements of a list recursively. 

updateList <-
    function (x, val)
{
    if (is.null(x))
        x <- list()
    modifyList(x, val)
}


## common operations that only make sense in certain contexts


useOuterStrips <-
    function(x,
             strip = strip.default,
             strip.left = strip.custom(horizontal = FALSE),
             strip.lines = 1,
             strip.left.lines = strip.lines)
{
    dimx <- dim(x)
    stopifnot(inherits(x, "trellis"))
    stopifnot(length(dimx) == 2)
    as.table <- x$as.table
    
    opar <- if (is.null(x$par.settings)) list() else x$par.settings
    par.settings <-
        modifyList(opar,
                   list(layout.heights =
                        if (as.table) list(strip = c(strip.lines, rep(0, dimx[2]-1)))
                        else list(strip = c(rep(0, dimx[2]-1), strip.lines)),
                        layout.widths =
                        list(strip.left = c(strip.left.lines, rep(0, dimx[1]-1)))))
    if (is.character(strip))
        strip <- get(strip)
    if (is.logical(strip) && strip)
        strip <- strip.default
    new.strip <-
        if (is.function(strip))
        {
            function(which.given, which.panel, var.name, ...) {
                row.to.keep <- if (as.table) 1 else nrow(trellis.currentLayout())
                if (which.given == 1 && current.row() == row.to.keep)
                    strip(which.given = 1,
                          which.panel = which.panel[1],
                          var.name = var.name[1],
                          ...)
            }
        }
        else strip # This could reasonable happen only if strip == FALSE
    if (is.character(strip.left))
        strip.left <- get(strip.left)
    if (is.logical(strip.left) && strip.left)
        strip.left <- strip.custom(horizontal = FALSE)
    new.strip.left <-
        if (is.function(strip.left))
        {
            function(which.given, which.panel, var.name, ...) {
                if (which.given == 2 && current.column() == 1)
                    strip.left(which.given = 1,
                               which.panel = which.panel[2],
                               var.name = var.name[2],
                               ...)
            }
        }
        else strip.left
    update(x,
           par.settings = par.settings,
           strip = new.strip,
           strip.left = new.strip.left,
           par.strip.text = list(lines = 0.5),
           layout = dimx)
}




resizePanels <-
    function(x, h = 1, w = 1)
{
    if (!missing(x))
        return(update(x,
                      par.settings =
                      list(layout.heights = list(panel = h),
                           layout.widths = list(panel = w))))

    cl <- trellis.currentLayout()
    if (all(dim(cl) > 1))
        stop("layout must have single column or single row.")
    if (all(dim(cl) == 1))
    {
        message("Nothing to be done.")
        return()
    }
    if (any(cl == 0)) stop("missing panels not allowed")
    if (dim(cl)[2] == 1) ## single column
    {
        pos <- seq(length = dim(cl)[1])
        heights <- 
            sapply(pos,
                   function(i) {
                       trellis.focus("panel", 1, i, highlight = FALSE)
                       ylim <- current.panel.limits()$ylim
                       trellis.unfocus()
                       diff(range(ylim))
                   })
        return(trellis.last.object(par.settings =
                                   list(layout.heights = list(panel = heights))))
    }
    else if (dim(cl)[1] == 1) ## single row
    {
        pos <- seq(length = dim(cl)[2])
        widths <- 
            sapply(pos,
                   function(i) {
                       trellis.focus("panel", i, 1, highlight = FALSE)
                       xlim <- current.panel.limits()$xlim
                       trellis.unfocus()
                       diff(range(xlim))
                   })
        return(trellis.last.object(par.settings =
                                   list(layout.widths = list(panel = widths))))
    }
    print(dim(cl))
    stop("shouldn't come here")
}




## utility functions to extract components of a formula.  Don't work
## reliably with unusual symbols

.responseName <- function(formula)
{
    if (length(formula) == 3) as.character(formula[2])
    else stop("invalid formula")
}

.covariateName <- function(formula)
{
    RHS <- 
        if (length(formula) == 3) as.character(formula[3])
        else if (length(formula) == 2) as.character(formula[2])
        else stop("invalid formula")
    RHS <- strsplit(RHS, " | ", fixed = TRUE)[[1]]
    RHS[1]
}

.groupsName <- function(formula)
{
    RHS <- 
        if (length(formula) == 3) as.character(formula[3])
        else if (length(formula) == 2) as.character(formula[2])
        else stop("invalid formula")
    RHS <- strsplit(RHS, " | ", fixed = TRUE)[[1]]
    RHS[2]
}


