BOTHlabels2 <- function (text, sel, XX, YY, adj, frame, pch, thermo, pie, piecol, 
    col, bg, horiz, width, height, ...) 
{
    if (missing(text)) 
        text <- NULL
    if (length(adj) == 1) 
        adj <- c(adj, 0.5)
    if (is.null(text) && is.null(pch) && is.null(thermo) && is.null(pie)) 
        text <- as.character(sel)
    frame <- match.arg(frame, c("rect", "circle", "none"))
    args <- list(...)
    CEX <- if ("cex" %in% names(args)) 
        args$cex
    else par("cex")
    if (frame != "none" && !is.null(text)) {
        if (frame == "rect") {
            width <- strwidth(text, units = "inches", cex = CEX)
            height <- strheight(text, units = "inches", cex = CEX)
            if ("srt" %in% names(args)) {
                args$srt <- args$srt%%360
                if (args$srt == 90 || args$srt == 270) {
                  tmp <- width
                  width <- height
                  height <- tmp
                }
                else if (args$srt != 0) 
                  warning("only right angle rotation of frame is supported;\n         try  `frame = \"n\"' instead.\n")
            }
            width <- xinch(width)
            height <- yinch(height)
            xl <- XX - width * adj[1] - xinch(0.03)
            xr <- xl + width + xinch(0.03)
            yb <- YY - height * adj[2] - yinch(0.02)
            yt <- yb + height + yinch(0.05)
            rect(xl, yb, xr, yt, col = bg)
        }
        if (frame == "circle") {
            radii <- 0.8 * apply(cbind(strheight(text, units = "inches", 
                cex = CEX), strwidth(text, units = "inches", 
                cex = CEX)), 1, max)
            symbols(XX, YY, circles = radii, inches = max(radii), 
                add = TRUE, bg = bg)
        }
    }
    if (!is.null(thermo)) {
        parusr <- par("usr")
        if (is.null(width)) {
            width <- CEX * (parusr[2] - parusr[1])
            width <- if (horiz) 
                width/15
            else width/40
        }
        if (is.null(height)) {
            height <- CEX * (parusr[4] - parusr[3])
            height <- if (horiz) 
                height/40
            else height/15
        }
        if (is.vector(thermo)) 
            thermo <- cbind(thermo, 1 - thermo)
        thermo <- if (horiz) 
            width * thermo
        else height * thermo
        if (is.null(piecol)) 
            piecol <- rainbow(ncol(thermo))
        xl <- XX - width/2 + adj[1] - 0.5
        xr <- xl + width
        yb <- YY - height/2 + adj[2] - 0.5
        yt <- yb + height
        if (horiz) {
            rect(xl, yb, xl + thermo[, 1], yt, border = NA, col = piecol[1])
            for (i in 2:ncol(thermo)) rect(xl + rowSums(thermo[, 
                1:(i - 1), drop = FALSE]), yb, xl + rowSums(thermo[, 
                1:i]), yt, border = NA, col = piecol[i])
        }
        else {
            rect(xl, yb, xr, yb + thermo[, 1], border = NA, col = piecol[1])
            for (i in 2:ncol(thermo)) rect(xl, yb + rowSums(thermo[, 
                1:(i - 1), drop = FALSE]), xr, yb + rowSums(thermo[, 
                1:i]), border = NA, col = piecol[i])
        }
        s <- apply(thermo, 1, function(xx) any(is.na(xx)))
        xl[s] <- xr[s] <- NA
        rect(xl, yb, xr, yt, border = "black")
        if (!horiz) {
            segments(xl, YY, xl - width/5, YY)
            segments(xr, YY, xr + width/5, YY)
        }
    }
    if (!is.null(pie)) {
        if (is.vector(pie)) 
            pie <- cbind(pie, 1 - pie)
        xrad <- CEX * diff(par("usr")[1:2])/50
        xrad <- rep(xrad, length(sel))
        XX <- XX + adj[1] - 0.5
        YY <- YY + adj[2] - 0.5
        for (i in 1:length(sel)) {
            if (any(is.na(pie[i, ]))) 
                next
            floating.pie.asp(XX[i], YY[i], pie[i, ], radius = xrad[i], 
                col = piecol)
        }
    }
    if (!is.null(text)) 
        text(XX, YY, text, adj = adj, col = col, ...)
    if (!is.null(pch)) 
        points(XX + adj[1] - 0.5, YY + adj[2] - 0.5, pch = pch, 
            col = col, bg = bg, ...)
}
