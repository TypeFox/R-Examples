plot.GPvam <-
function (x, ..., alpha=.1) 
{
#plotCI is taken from R package gplots
plotCI<-function (x, y = NULL, uiw, liw = uiw, ui, li, err = "y", ylim = NULL, 
    xlim = NULL, type = "p", col = par("col"), barcol = col, 
    pt.bg = par("bg"), sfrac = 0.01, gap = 1, lwd = par("lwd"), 
    lty = par("lty"), labels, add = FALSE, xlab, ylab, 
    minbar, maxbar, ...) 
{
    if (is.list(x)) {
        y <- x$y
        x <- x$x
    }

    if (is.null(y)) {
        if (is.null(x)) 
            stop("both x and y NULL")
        y <- as.numeric(x)
        x <- seq(along = x)
    }
    li<-y - uiw
    ui<-y + uiw
    if (err == "y") 
        z <- y
    else z <- x
    if (FALSE) 
        uiw <- NA
    if (FALSE) 
        liw <- NA
    if (FALSE) 
        ui <- z + uiw
    if (FALSE) 
        li <- z - liw
    if (FALSE) 
        li <- ifelse(li < minbar, minbar, li)
    if (FALSE) 
        ui <- ifelse(ui > maxbar, maxbar, ui)
    if (err == "y") {
        if (is.null(ylim)) 
            ylim <- range(c(y, ui, li), na.rm = TRUE)
        if (is.null(xlim) && !is.R()) 
            xlim <- range(x, na.rm = TRUE)
    }
    else if (err == "x") {
        if (is.null(xlim)) 
            xlim <- range(c(x, ui, li), na.rm = TRUE)
        if (is.null(ylim) && !is.R()) 
            ylim <- range(x, na.rm = TRUE)
    }
    if (!add) {
        if (FALSE) 
            plot(x, y, ylim = ylim, xlim = xlim, col = col, xlab = xlab, 
                ylab = ylab, ...)
        else {
            plot(x, y, ylim = ylim, xlim = xlim, col = col, type = "n", 
                xlab = xlab, ylab = ylab, ...)
            text(x, y, label =labels, col = col, ...)
        }
    }
    if (err == "y") {
        if (gap != FALSE) 
            gap <- strheight("O") * gap
        smidge <- par("fin")[1] * sfrac
        if (!is.null(li)) 
            suppressWarnings(arrows(x, li, x, pmax(y - gap, li), col = barcol, 
                lwd = lwd, lty = lty, angle = 90, length = smidge, 
                code = 1))
        if (!is.null(ui)) 
            suppressWarnings(arrows(x, ui, x, pmin(y + gap, ui), col = barcol, 
                lwd = lwd, lty = lty, angle = 90, length = smidge, 
                code = 1))
    }
    else {
        if (gap != FALSE) 
            gap <- strwidth("O") * gap
        smidge <- par("fin")[2] * sfrac
        if (!is.null(li)) 
            arrows(li, y, pmax(x - gap, li), y, col = barcol, 
                lwd = lwd, lty = lty, angle = 90, length = smidge, 
                code = 1)
        if (!is.null(ui)) 
            arrows(ui, y, pmin(x + gap, ui), y, col = barcol, 
                lwd = lwd, lty = lty, angle = 90, length = smidge, 
                code = 1)
    }
    points(x, y, col = col, lwd = lwd, bg = pt.bg, type = type, 
        ...)
    invisible(list(x = x, y = y))
}

    devAskNewPage(ask = TRUE)
   c.level <- qnorm(1 - alpha/2)
    for (i in 1:length(unique(x$teach.effects$teacher_year))) {
      if(x$persistence=="GP"){
        for (j in i:length(unique(x$teach.effects$teacher_year))) {
            temp.df <- x$teach.effects[x$teach.effects$teacher_year == 
                i & x$teach.effects$effect_year == j, ]
            temp.df <- temp.df[order(temp.df$EBLUP), ]
            plotCI(temp.df$EBLUP, uiw = c.level * temp.df$std_error, labels=temp.df$teacher,
                barcol = 2, xlab = "Ranked Teachers", ylab = "Teacher Effect")
            title(paste("Year ", i, " Teachers' Year ", j," Effect\nwith ",(1-alpha)*100,"% Confidence Intervals", sep = ""))
            abline(h = 0)
        }
        }else   if(x$persistence=="rGP"){
        for (j in i:min((i+1),length(unique(x$teach.effects$teacher_year)))) {
            temp.df <- x$teach.effects[x$teach.effects$teacher_year == 
                i & x$teach.effects$effect_year == j, ]
            temp.df <- temp.df[order(temp.df$EBLUP), ]
            plotCI(temp.df$EBLUP, uiw = c.level * temp.df$std_error, labels=temp.df$teacher,
                barcol = 2, xlab = "Ranked Teachers", ylab = "Teacher Effect")
            title(paste("Year ", i, " Teachers' Year ", j," Effect\nwith ",(1-alpha)*100,"% Confidence Intervals", sep = ""))
            abline(h = 0)
        }
        }else{
        for (j in i) {
            temp.df <- x$teach.effects[x$teach.effects$teacher_year == 
                i & x$teach.effects$effect_year == j, ]
            temp.df <- temp.df[order(temp.df$EBLUP), ]
            plotCI(temp.df$EBLUP, uiw = c.level * temp.df$std_error,   labels=temp.df$teacher,
                barcol = 2, xlab = "Ranked Teachers", ylab = "Teacher Effect")
            title(paste("Year ", i, " Teachers' Year ", j," Effect\nwith ",(1-alpha)*100,"% Confidence Intervals", sep = ""))
            abline(h = 0)
        }
        }
        
    }
    qqnorm(x$cresid, main = "Normal Q-Q Plot\n for raw conditional residuals")
    qqline(x$cresid)
    qqnorm(x$sresid, main = "Normal Q-Q Plot\n for scaled conditional residuals")
    qqline(x$sresid)
    plot(x$yhat.s, x$sresid, main = "Scaled Conditional Residuals\n(by inverse Cholesky root\nof conditional error matrix)", 
        xlab = "Predicted", ylab = "Residuals")
    plot(x$yhat, x$cresid, main = "Raw conditional residuals", 
        xlab = "Predicted", ylab = "Residuals")
    plot(x$yhat.m, x$mresid, main = "Raw marginal residuals", 
        xlab = "Predicted", ylab = "Residuals")
    devAskNewPage(ask = FALSE)
    invisible(x)
}
