ceplot.static <-
function (data, model, response = NULL, S = NULL, C = NULL, sigma = NULL, 
    distance = "euclidean", cex.axis = NULL, cex.lab = NULL, tck = NULL, 
    view3d = FALSE, theta3d = 45, phi3d = 20, Corder = "default", 
    Xc.cond = NULL, select.colour = "blue", select.cex = 1)
{
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
        #cat(paste("\n'S' was not specified, picked", colnames(data)[S]))
        } else if (is.character(S))
            vapply(S, function(x) which(colnames(data) == x), numeric(1))
            else S
    C <- if (is.null(C))
        arrangeC(data[, -c(response, S)])
    else C
    C <- if (all(vapply(C, is.numeric, logical(1))))
        as.list(C)
    else if (all(vapply(C, is.character, logical(1))))
            lapply(C, match, table = colnames(data))
        else
            stop("'C' should be a vector or list (containing vectors of length",
                 " 1 or 2) with integer column indices or character variable",
                 " names from 'data'.")
    uniqC <- unique(unlist(C))
    if (any(response %in% uniqC))
        stop("cannot have 'response' variable in 'C'")
    if (any(response %in% S))
        stop("cannot have 'response' variable in 'S'")
    if (!identical(length(unique(vapply(lapply(model, getvarnames), 
        `[[`, character(1), 1))), 1L))
        stop("cannot compare models with different response variables") 
    if (!identical(length(intersect(S, uniqC)), 0L))
        stop("cannot have variables common to both 'S' and 'C'")    
    Xc.cond <- if (is.null(Xc.cond))
        data[1, uniqC, drop = FALSE]
    else Xc.cond  
    xcplots <- list()
    close.screen(all.screens = T)
    n.selector.cols <- ceiling(length(C) / 4L)
    selector.colwidth <- 2
    height <- 8
    width <- height + 0.5 + selector.colwidth * n.selector.cols
    xcwidth <- selector.colwidth * n.selector.cols / width
    main <- split.screen(figs = matrix(c(0, 1 - xcwidth, 1 - xcwidth, 1, 
        0, 0, 1, 1), ncol = 4))
    selectors <- split.screen(figs = c(4, n.selector.cols), screen = main[2])
    dev.hold()
    if (length(uniqC) > 0){
        for(C.index in seq_along(C)){
            screen(selectors[C.index])
            xcplots[[C.index]] <- plotxc(xc = data[, C[[C.index]]], xc.cond = 
                Xc.cond[1L, colnames(data)[C[[C.index]]]], name = colnames(data)
                [C[[C.index]]], select.colour = select.colour, select.lwd = 2, 
                cex.axis = cex.axis, cex.lab = cex.lab, tck = tck, select.cex = 
                select.cex)
        }
    }
    screen(main[1])
    Xc <- data[, uniqC, drop = FALSE]
    vw <- visualweight(xc = Xc, xc.cond = Xc.cond, sigma = sigma, distance = 
        distance)
    k <- vw$k
    data.colour <- rgb(1 - k, 1 - k, 1 - k)
    data.order <- vw$order
    xsplot <- plotxs(xs = data[, S, drop = FALSE],
        y = data[, response, drop = FALSE], xc.cond = Xc.cond, model = model,
        model.colour = NULL, model.lwd = NULL, model.lty = NULL,
        model.name = model.name, yhat = NULL, mar = NULL, 
        data.colour = data.colour, data.order = data.order, view3d = view3d, 
        theta3d = theta3d, phi3d = phi3d)
    dev.flush()
    structure(list(Xc = Xc, sigma = sigma, distance = distance, xcplots = 
        xcplots, xsplot = xsplot, screens = list(main = main, selectors = 
        selectors)), class = "ceplot")
}