.mousemove.separate.full <-
function (buttons, x, y)
{
    if (version$os != "linux-gnu"){
    if (exists("buttons")){
	if (0 %in% buttons){
    expectationwindow <- get("expectationwindow", envir = parent.frame())
    selectorwindow <- get("selectorwindow", envir = parent.frame())
    plotxcobject <- get("plotxcobject", envir = parent.frame())
    plotxsobject <- get("plotxsobject", envir = parent.frame())
    Xc <- get("Xc", envir = parent.frame())
    Xc.cond <- plotxsobject$xc.cond
    factorindex <- vapply(Xc, is.factor, logical(1))
    Xc.num <- vapply(Xc, as.numeric, numeric(nrow(Xc)))
    Xc.cond.num <- vapply(Xc.cond, as.numeric, numeric(1L))

    scr2 <- plotxcobject$scr2
    rows <- plotxcobject$rows
    cols <- plotxcobject$cols
    coords <- plotxcobject$coords

    tmp <- xy2index(x, y, coords)
    i <- which(tmp == scr2)
    if (!identical(rows[i], cols[i])){
    par(bg = "white")
    screen(tmp)
    plot(Xc.num[,cols[i]], Xc.num[,rows[i]], cex = 0.6, xlab = "", ylab = "", 
        xaxt = "n", yaxt = "n", col = if (identical(rows[i], cols[i])) NULL 
        else "black")
    click <- c(grconvertX(x, "ndc", "user"), grconvertY(y, "ndc", "user"))  

    if (cols[i] %in% which(!factorindex)){
        Xc.cond[cols[i]] <- Xc.cond.num[cols[i]] <- click[1]
    } else {
        Xc.cond[cols[i]] <- factor(levels(Xc[, cols[i]])[which.min(abs(click[1]
            - 1:length(levels(Xc[, cols[i]]))))], 
            levels = levels(Xc[, cols[i]]))
        Xc.cond.num[cols[i]] <- as.numeric(Xc.cond[cols[i]])  
    }
    if (rows[i] %in% which(!factorindex)){
        Xc.cond[rows[i]] <- Xc.cond.num[rows[i]] <- click[2]
    } else {
        Xc.cond[rows[i]] <- factor(levels(Xc[, rows[i]])[which.min(abs(click[2]
            - 1:length(levels(Xc[, rows[i]]))))], 
            levels = levels(Xc[, rows[i]]))
        Xc.cond.num[rows[i]] <- as.numeric(Xc.cond[rows[i]])  
    }  
    for (sc in which(rows == rows[i] | cols == rows[i] | rows == cols[i] | cols 
        == cols[i])){
        if (!identical(rows[sc], cols[sc])){
            screen(scr2[sc])
            plot(Xc.num[,cols[sc]], Xc.num[,rows[sc]], cex = 0.6, xlab = "", 
                ylab = "", xaxt = "n", yaxt = "n", col = if (identical(rows[i], 
                cols[i])) NULL else "black")
            abline(v = Xc.cond.num[cols[sc]], h = Xc.cond.num[rows[sc]], col = 
                plotxcobject$select.colour)
        }
    }
    dev.flush()
    
    dev.set(expectationwindow)
    par(bg = "white")
    close.screen(all.screens = TRUE)
    plotxsobject$xc.cond <- Xc.cond
    vwargs <- get("vwargs", envir = parent.frame())
    vw <- visualweight(xc = Xc, xc.cond = Xc.cond, 
        sigma = vwargs$sigma, distance = vwargs$distance)
    k <- vw$k
    plotxsobject$data.colour <- rgb(1-k,1-k,1-k)
    plotxsobject$data.order <- vw$order
    par(bg = "white")
    close.screen()
    do.call(plotxs, plotxsobject)
    assign(x = "plotxsobject", value = plotxsobject, envir = parent.frame())  
    assign(x = "plotxcobject", value = plotxcobject, envir = parent.frame())
    dev.set(selectorwindow)
    }
    points(NULL)
    }}}
}
