condtour <-
function(data, model, path, response = NULL, S = NULL, C = NULL, sigma = NULL, 
    distance = "euclidean", cex.axis = NULL, cex.lab = NULL, tck = NULL, view3d 
    = FALSE, conf = FALSE, select.colour = "blue")
{
    xold <- NULL
    yold <- NULL 
    mousemove <- function ()
    {
        function (buttons, x, y)
        {
            if (all(findInterval(x, xscoords[1:2]) == 1, identical(
                xsplot$plot.type, "ccc"), xsplot$view3d, 0 %in% buttons)){
                if (!is.null(xold))
                    xsplot <<- update(xsplot, theta3d = xsplot$theta3d + 1 * 
                        (xold > x) - 1 * (xold < x), phi3d = xsplot$phi3d + 1 * 
                        (yold > y) - 1 * (yold < y), xs.grid = xsplot$xs.grid, 
                        prednew = xsplot$prednew)
                xold <<- x
                yold <<- y                    
            }
            points(NULL)
        }
    }
    mouseclick <- function ()
    {
        function (buttons, x, y)
        {
            if (0 %in% buttons){
                pathindex <<- max(min(pathindex + 1, max(pathindexrange)), 
                    min(pathindexrange))
                applot <<- update(applot, pathindex = pathindex)    
                xc.cond <- path[pathindex, , drop = FALSE]
                k.order <- order(k[pathindex, ])
                k.order.trimmed <- k.order[k[pathindex, ][k.order] > 0]
                xsplot <<- update(xsplot, xc.cond = xc.cond, data.colour = 
                    rgb(1 - k[pathindex, ], 1 - k[pathindex, ], 1 - k[pathindex, 
                    ]), data.order = k.order.trimmed) 
                for (i in seq_along(C)){
                    xcplots[[i]] <<- update(xcplots[[i]], xc.cond = path[
                        pathindex, colnames(data)[C[i]]])
                }                  
            }
            points(NULL)
        }
    }
    keystroke <- function ()
    {
        function (key)
        {
            if (identical(key, "q")){
                cat("\nInteractive session ended.\n")
                return(invisible(1))            
            } 
            if (identical(xsplot$plot.type, "ccc") & xsplot$view3d & 
                key %in% c("Up", "Down", "Left", "Right")){
                xsplot <<- update(xsplot, theta3d = xsplot$theta3d - 2 * 
                    (key == "Right") + 2 * (key == "Left"), phi3d = xsplot$phi3d 
                    - 2 * (key == "Up") + 2 * (key == "Down"), xs.grid = 
                    xsplot$xs.grid, prednew = xsplot$prednew)                
            }
            if (key %in% c("[", "]")){
                pathindex <<- max(min(pathindex + 1 * (key == "]") - 1 * (key == 
                    "["), max(pathindexrange)), min(pathindexrange))
                applot <<- update(applot, pathindex = pathindex)    
                xc.cond <- path[pathindex, , drop = FALSE]
                #vw <<- visualweight(xc = data[, uniqC, drop = FALSE], 
                #    xc.cond = xc.cond, sigma = sigma, distance = distance)
                k.order <- order(k[pathindex, ])
                k.order.trimmed <- k.order[k[pathindex, ][k.order] > 0]
                xsplot <<- update(xsplot, xc.cond = xc.cond, data.colour = 
                    rgb(1 - k[pathindex, ], 1 - k[pathindex, ], 1 - k[pathindex, 
                    ]), data.order = k.order.trimmed) 
                for (i in seq_along(C)){
                    xcplots[[i]] <<- update(xcplots[[i]], xc.cond = path[
                        pathindex, colnames(data)[C[i]]])
                }               
            }            
            points(NULL)
        }
    } 
    data <- na.omit(data)
    model <- if (!identical(class(model), "list"))
        list(model)
    else model
    model.name <- if (!is.null(names(model)))
        names(model)
    else NULL
    varnamestry <- try(getvarnames(model[[1]]), silent = TRUE)
    response <- if (is.null(response))
        if (class(varnamestry) != "try-error")
           which(colnames(data) == varnamestry$response[1])
        else stop("could not extract response from 'model'.")
    else if (is.character(response))
            which(colnames(data) == response)
        else response
    S <- if(is.null(S)){
         (1:ncol(data))[-response][1L]
        } else if (is.character(S))
            vapply(S, function(x) which(colnames(data) == x), numeric(1))
            else S       
    C <- if (is.null(C))
        arrangeC(data[, -c(response, S)])
    else C
    try(
        if (class(varnamestry) != "try-error"){
            possibleC <- unique(unlist(lapply(lapply(model, getvarnames), `[[`, 
                2)))
            possibleC <- possibleC[possibleC %in% colnames(data)]
            C <- arrangeC(data[, possibleC[!(possibleC %in% colnames(data)[S])], 
                drop = FALSE])
        }     
    , silent = TRUE) 
    C <- if (all(vapply(C, is.numeric, logical(1))))
        as.list(C)
    else if (all(vapply(C, is.character, logical(1))))
            lapply(C, match, table = colnames(data))
        else
            stop("'C' should be a vector or list (containing vectors of length",
                 " 1 or 2) with integer column indices or character variable",
                 " names from 'data'.")
    uniqC <- unique(unlist(C))
    C <- uniqC
    pathindex <- 1
    pathindexrange <- c(1, nrow(path))
    xc.cond <- path[pathindex, , drop = FALSE]
    if (any(response %in% uniqC))
        stop("cannot have 'response' variable in 'C'")
    if (any(response %in% S))
        stop("cannot have 'response' variable in 'S'")
    if (!identical(length(intersect(S, uniqC)), 0L))
        stop("cannot have variables common to both 'S' and 'C'")
    xcplots <- list()
    coords <- matrix(ncol = 4L, nrow = length(C))    
    plotlegend <- length(S) == 2
    n.selector.cols <- ceiling(length(C) / 4L)
    selector.colwidth <- 2
    height <- 8    
    width <- height + 0.5 * plotlegend
    k <- matrix(ncol = nrow(data), nrow = nrow(path))
    for (i in 1: nrow(path)){
        k[i, ] <- visualweight(xc.cond = path[i, , drop = F], xc = data[, 
            colnames(path), drop = FALSE], sigma = sigma, basicoutput = T)
    }
    opendev(width = width, height = height)
    #if (identical(version$os, "linux-gnu"))
    #    x11(type = "Xlib", height = height, width = height + 0.5 * plotlegend)
    #else x11(height = height, width = height + 0.5 * plotlegend)
    devexp <- dev.cur()    
    close.screen(all.screens = TRUE)    
    
    legendwidth <- 1.15 / height
    xsscreens <- if (plotlegend){
        split.screen(figs = matrix(c(0, 1 - legendwidth, 1 - legendwidth, 1, 
            0, 0, 1, 1), ncol = 4))
    } else split.screen()
    if (plotlegend){
        screen(xsscreens[2L])
        xslegend(data[, response], colnames(data)[response])
    }
    screen(xsscreens[1L])
    k.order <- order(k[pathindex, ])
    k.order.trimmed <- k.order[k[pathindex, ][k.order] > 0]
    par(mar = c(3, 3, 3, 3))
    xsplot <- plotxs1(xs = data[, S, drop = FALSE], data[, response, 
        drop = FALSE], xc.cond = xc.cond, model = model, data.colour = rgb(1 - 
        k[pathindex, ], 1 - k[pathindex, ], 1 - k[pathindex, ]), data.order = 
        k.order.trimmed, view3d = view3d, conf = conf)
    xscoords <- par("fig") 
    xcwidth <- selector.colwidth * n.selector.cols
    n.selector.rows <- ceiling(length(C) / n.selector.cols)
    xcheight <- selector.colwidth * n.selector.rows
    setGraphicsEventHandlers(
        onMouseMove = mousemove(),
        onKeybd = keystroke()) 
    opendev(width = xcwidth, height = xcheight)    
    #if (identical(version$os, "linux-gnu"))
    #    x11(type = "Xlib", height = xcheight, width = xcwidth)
    #else x11(height = 5, width = 3)
    devdiag <- dev.cur()    
    close.screen(all.screens = TRUE)
    diagscreens <- split.screen(c(2, 1))
    screen(diagscreens[1L])
    par(mar = c(4, 4, 2, 2))
    plotmaxk(apply(k, 2, max))
    screen(diagscreens[2L])
    par(mar = c(4, 4, 2, 2))
    applot <- plotap(k)
    setGraphicsEventHandlers(
        onMouseDown = mouseclick(),
        onKeybd = keystroke())
    opendev(width = xcwidth, height = xcheight)    
    #if (identical(version$os, "linux-gnu"))
    #    x11(type = "Xlib", height = xcheight, width = xcwidth)
    #else x11(height = xcheight, width = xcwidth)
    devcond <- dev.cur()    
    close.screen(all.screens = TRUE)    
    xcscreens <- split.screen(c(n.selector.rows, n.selector.cols))
    for (i in seq_along(uniqC)){
        screen(xcscreens[i])
        xcplots[[i]] <- plotxc(xc = data[, C[[i]]], xc.cond = path[
                        pathindex, colnames(data)[C[i]]], 
            name = colnames(data[, C[[i]], drop = FALSE]), select.colour = 
            select.colour)
        coords[i, ] <- par("fig")
    }  
    setGraphicsEventHandlers(
        onMouseDown = mouseclick(),
        onKeybd = keystroke())
    getGraphicsEventEnv()
    getGraphicsEvent()    
}