plot.gantt <- function (x,
                        xlim,
                        time.format = NULL,
                        time.labels.by, # = "2 months",
                        time.lines.by,  # = "1 month",
                        event.time = NULL,
                        event.label = NULL,
                        event.side = 3,
                        col.done = gray(0.3),
                        col.notdone = gray(0.9),
                        col.event = gray(0.1),
                        col.connector = "black",
                        grid.col = "lightgray", grid.lty = "dotted",
                        main = "",
                        cex=par("cex"),
                        debug=FALSE,
                        ...)
{
    if (!inherits(x, "gantt")) stop("method is only for gantt objects")
    opar <- par(no.readonly = TRUE)

    mgp <- c(3, 1, 0)
    half.height <- 0.33
    t0 <- as.POSIXct("1970-01-01 00:00:00")
    ndescriptions <- length(x$description)
    charheight <- strheight("M", units = "inches")
    maxwidth <- max(strwidth(x$description, units = "inches")) * 1.1

    ## Get around some problems with autoscaling of POSIXt values
    r <- if (missing(xlim)) range(x$start, x$end) else xlim
    if (debug) {cat("range: ", as.character(r[1]), "to", as.character(r[2]), "\n")}
    s <- as.numeric(difftime(r[2], r[1], units="days"))
    r <- as.POSIXlt(r)
    subTics <- NULL
    if (s > 100) {
        if (is.null(time.format)) time.format <-  "%b %Y" # month/year
        r[2] <- r[2] + 86400
        r[1:2]$hour <- r[1:2]$min <- r[1:2]$sec <- 0
        if (debug){cat("range: ", as.character(r[1]), "to", as.character(r[2]), "\n")}
        ## monthly ticks
        lhs <- as.POSIXlt(r[1])
        lhs$mon <- 0
        lhs$mday <- 1
        rhs <- as.POSIXlt(r[2])
        rhs$mon <- 11
        rhs$mday <- 31
        subTics <- seq(lhs, rhs, by="month")
    } else {
        if (s > 10) {
            if (is.null(time.format)) time.format <-  "%d/%b" # day/month
            r[2] <- r[2] + 86400
            r[1:2]$hour <- r[1:2]$min <- r[1:2]$sec <- 0
            if(debug){cat("range: ", as.character(r[1]), "to", as.character(r[2]), "\n")}
        } else {
            if (s > 1) {
                if (is.null(time.format)) time.format <-  "%d/%b" # day/month
                r[2] <- r[2] + 86400
                r[1:2]$hour <- r[1:2]$min <- r[1:2]$sec <- 0
                if(debug){cat("range: ", as.character(r[1]), "to", as.character(r[2]), "\n")}
            } else {
                if (is.null(time.format)) time.format <-  "%d/%b" # day/month
            }
        }
    }

    bottom.margin <- 0.5
    par(mai = c(bottom.margin, maxwidth, charheight * 2, 0.1))
    par(omi = c(0.1, 0.1, 0.1, 0.1), xaxs = "i", yaxs = "i")
    plot(c(r[1], r[2]), c(1,2*ndescriptions),
         ylim=c(0.5, ndescriptions + 0.5),
         main = "", xlab = "", ylab = "", xaxs="r", type="n", axes = FALSE)
    xlim <- as.POSIXct(par("usr")[1:2] + t0)
    box()
    if (nchar(main)) mtext(main, 3, 2, cex=cex)
    if (missing(time.labels.by)) {
        xaxp <- par("xaxp")
        lines.at.0 <- axis.POSIXct(1,
                                   at=pretty(r, 10), #seq(xaxp[1], xaxp[2], length.out=xaxp[3]) + t0,
                                   format=time.format, cex.axis=cex, ...)
    } else {
        lines.at.0 <- axis.POSIXct(1,
                                   at=as.POSIXct(seq.POSIXt(as.POSIXct(xlim[1]), as.POSIXct(xlim[2]), by=time.labels.by)),
                                   format=time.format, cex.axis=cex, ...)
    }
    if (!is.null(subTics))
        rug(subTics)
    if (missing(time.lines.by)) {
        abline(v=lines.at.0, col = grid.col, lty=grid.lty)
    } else {
        abline(v = seq.POSIXt(as.POSIXct(xlim[1]), as.POSIXct(xlim[2]), by=time.lines.by), col = grid.col, lty=grid.lty)
    }
    if (FALSE) {
        mtext(paste(paste(format(xlim[1:2], format="%Y-%m-%d"), collapse=" to "),
                    attr(x$data$ts$time[1], "tzone")),
              side=3, cex=2/3, adj=0)
    }
    topdown <- seq(ndescriptions, 1)
    axis(2, at = topdown, labels = x$description, las = 2, tick=FALSE, cex.axis=cex)

    ## Connectors
    for (t in 1:ndescriptions) {
        nb <- x$neededBy[t][[1]]
        if (!is.na(nb)) {
            source.y <- topdown[t]
            source.t <- as.POSIXct(x$end[t])
            for (nbi in 1:length(nb)) {
                r <- as.numeric(nb[nbi])
                receiver.t <- as.POSIXct(x$start[r])
                receiver.y <- topdown[r]
                lines(c(source.t,receiver.t), c(source.y,receiver.y),col=col.connector)
            }
        }
    }
    ## Events
    if (!is.null(event.time)) {
        ne <- length(event.time)
        for (e in 1:ne) {
            t <- as.POSIXct(event.time[e])
            abline(v=t, col=col.event)
            mtext(event.label[e], side=event.side, col=col.event, at=t)
        }
    }
    ## Description
    for (i in 1:ndescriptions) {
        mid <- as.POSIXct(x$start[i]) +
            x$done[i] * as.numeric(difftime(as.POSIXct(x$end[i]),
                                            as.POSIXct(x$start[i]),
                                            units="secs")) / 100
        if (debug) {cat(as.character(x$description[i]),"takes", as.numeric(difftime(as.POSIXct(x$end[i]), as.POSIXct(x$start[i]), units="secs")), "s\n")}

        bottom <- topdown[i] - half.height
        top <- topdown[i] + half.height
        left <- as.POSIXct(x$start[i])
        right <- as.POSIXct(x$end[i])

        if (debug){cat(as.character(x$description[i]));cat(" done=",x$done[i]," mid=");print(mid);cat(" left=");print(left);cat("right=");print(right);cat("\n")}

        rect(left, bottom, right, top, col = col.notdone,   border = FALSE)
        rect(left, bottom, mid,   top, col = col.done,      border = FALSE)
        rect(left, bottom, right, top, col = "transparent", border = TRUE)
    }
    abline(h = (topdown[1:(ndescriptions - 1)] + topdown[2:ndescriptions])/2,  col = grid.col, lty=grid.lty)
    invisible(x)
}

summary.gantt <- function(object, ...)
{
    if (!inherits(object, "gantt")) stop("method is only for ganttobjects")
    res <- object
    class(res) <- "summary.gantt"
    res
}

print.summary.gantt <- function(x, ...)
{
    if (!inherits(x, "summary.gantt")) stop("method is only for summary.gantt objects")
    max.description.width <- max(nchar(as.character(x$description)))
    num.descriptions <- length(x$description)
    cat("Key, Description,", paste(rep(" ", max.description.width-12), collapse=""), "Start,      End,        Done, NeededBy\n")
    for (t in 1:num.descriptions) {
        spacer <- paste(rep(" ", 1 + max.description.width - nchar(as.character(x$description[t]))),
                        collapse="")
        cat(paste(format(x$key[t], width=3, justify="right"), ",", sep=""),
            paste(as.character(x$description[t]), ",",
                  spacer,
                  format(x$start[t]), ", ",
                  x$end[t],  ", ",
                  format(x$done[t], width=4, justify="right"), sep = ""))
        nb <- x$neededBy[t][[1]]
        if (!is.null(nb) && !is.na(nb[1])) {
            cat(", ")
            for (nbi in 1:length(nb)) {
                cat(x$description[as.numeric(nb[nbi])], " ")
            }
        }
        cat("\n")
    }
}

as.gantt <- function(key, description, start, end, done, neededBy)
{
    if (missing(key))
        stop("must give 'key'")
    if (missing(description))
        stop("must give 'description'")
    if (missing(start))
        stop("must give 'start'")
    if (missing(end))
        stop("must give 'end'")
    n <- length(key)
    if (missing(done))
        done <- rep(0, n)
    if (missing(neededBy))
        neededBy <- rep(NA, n)
    rval <- list(key=key,
                 description=as.character(description),
                 start=as.POSIXct(start),
                 end=as.POSIXct(end),
                 done=done,
                 neededBy=neededBy)
    class(rval) <- "gantt"
    rval
}

read.gantt <- function(file, debug=FALSE)
{
    if (is.character(file)) {
        file <- file(file, "r")
        on.exit(close(file))
    }
    if (!inherits(file, "connection")) stop("argument `file' must be a character string or connection")
    if (!isOpen(file)) {
        open(file, "r")
        on.exit(close(file))
    }
    quiet <- !debug
    tokens <- trim.whitespace(scan(file,what='char',sep=",",nlines=1,quiet=quiet))
    check.tokens(tokens, c("Key", "Description", "Start", "End", "Done", "NeededBy"))
    key <- description <- start <- end <- done <- neededBy <- c()
    while (TRUE) {
        tokens <- trim.whitespace(scan(file, what=character(0), nlines=1,
                                       blank.lines.skip=FALSE, quiet=quiet, sep=","))
        ni <- length(tokens)
        if (ni > 1) {
            if (ni < 3) stop("need at least 3 items per line")
            key <- c(key, as.numeric(tokens[1]))
            description <- c(description, tokens[2])
            start <- c(start, tokens[3])
            end <- c(end, tokens[4])
            done <- c(done, if (ni >= 5) as.numeric(tokens[5]) else NA)
            neededBy <- c(neededBy, if (ni >= 6) as.numeric(tokens[6:ni]) else NA)
        } else {
            break
        }
    }
    as.gantt(key=key,
             description=as.character(description),
             start=as.POSIXct(start),
             end=as.POSIXct(end),
             done=done,
             neededBy=neededBy)
}
