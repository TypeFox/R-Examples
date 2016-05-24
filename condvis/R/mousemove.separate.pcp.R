.mousemove.separate.pcp <-
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
    xcoord <- plotxcobject$xcoord
    ycoord <- plotxcobject$ycoord

    click <- c(grconvertX(x, "ndc", "user"), grconvertY(y, "ndc", "user"))
    d <- abs((xcoord - click[1L]))
    index <- which.min(d)
    if (min(d) < 0.3){
        if (index %in% which(!factorindex)){
            ycoord[index] <- max(0, min(click[2], 1))
            Xc.cond[index] <- ycoord[index] * (max(Xc.num[, index]) - 
                min(Xc.num[, index])) + min(Xc.num[, index])
            xc.cond.new <- Xc.cond[index]
            varnames <- names(Xc.cond)[index]            
        } else {
            tmp <- click[2] * (max(Xc.num[, index]) - min(Xc.num[, index])) + 
                min(Xc.num[, index])
            Xc.cond[index] <- factor(levels(Xc[, index])[which.min(abs(tmp - 
                1:length(levels(Xc[, index]))))], levels = levels(Xc[, index]))
            tmp2 <- as.numeric(Xc.cond[index])
            tmp3 <- (tmp2 - min(Xc.num[, index])) / diff(range(Xc.num[, index]))  
            ycoord[index] <- tmp3          
            xc.cond.new <- Xc.cond[index]
            varnames <- names(Xc.cond)[index] 
            } 
    }
    parcoord(Xc.num, main = "Condition selector")
    points(xcoord, ycoord, col = plotxcobject$select.colour, type = "l", lwd = 
        2)
    points(xcoord, ycoord, col = plotxcobject$select.colour, pch = 16)

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
        plotxcobject$ycoord <- ycoord
        assign(x = "plotxcobject", value = plotxcobject, envir = parent.frame())
    dev.set(selectorwindow)
    points(NULL)
    }}}
}
