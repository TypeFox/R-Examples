## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

### Automatic playwith support

autoplay <-
    function(on=NA, lattice.on=on, base.on=on, grid.on=on, ask=FALSE)
{
    if (all(is.na(c(lattice.on, base.on, grid.on))))
        message("No action taken.")
    if (!is.na(lattice.on)) {
        library("lattice")
        if (packageDescription("lattice")$Version < package_version("0.17-1"))
            stop("this requires lattice package version >= 0.17")
        newFoo <- if (lattice.on) playwith.trellis else NULL
        lattice.options(print.function = newFoo)
        message("Automatic `playwith` for Lattice graphics is now ",
                if (lattice.on) "ON." else "OFF.")
    }
    if (!is.na(base.on)) {
        newFoo <- if (base.on) list(playwith.plot.new) else NULL
        setHook("before.plot.new", newFoo, "replace")
        message("Automatic `playwith` for base graphics is now ",
                if (base.on) "ON." else "OFF.")
        if (base.on) {
            StateEnv$.autoplay.ask <- ask
            if (ask) message("The high-level plot call will be chosen from a list.")
        }
    }
    if (!is.na(grid.on)) {
        newFoo <- if (grid.on) list(playwith.grid.newpage) else NULL
        setHook("grid.newpage", newFoo, "replace")
        message("Automatic `playwith` for grid graphics is now ",
                if (grid.on) "ON." else "OFF.")
        if (grid.on) {
            StateEnv$.autoplay.ask <- ask
            if (ask) message("The high-level plot call will be chosen from a list.")
        }
    }
    invisible()
}

## not to be called directly
playwith.trellis <-
    function(x, position = NULL, split = NULL, more = FALSE, newpage = TRUE,
             packet.panel = packet.panel.default, draw.in = NULL, ...)
{
    dev.interactive2 <- function(orNone)
    {
        dev.interactive(orNone) ||
        (interactive() && .Device == "null device" &&
         getOption("device") == "Cairo")
    }
    new <- (newpage && is.null(draw.in) &&
            !lattice:::lattice.getStatus("print.more"))
    if (dev.interactive2(TRUE) && new) {
        ## starting a new plot on an interactive device
        eval(call("playwith", x$call, envir=parent.frame(2)))
        return(invisible())
    }
    ## call `plot.trellis` from lattice package, as usual
    ocall <- sys.call()
    ocall[[1]] <- quote(plot)
    eval.parent(ocall)
}

## not to be called directly
playwith.plot.new <- function(...)
{
    sysCallNames <- sapply(sys.calls(), function(x)
                           ifelse(is.symbol(x[[1]]), toString(x[[1]]), ""))
    playing <- any(c("playReplot", "playNewPlot")
                   %in% sysCallNames)
                                        #multifig <- !isTRUE(all.equal(par("mfrow"), c(1,1)))
    rcmdr <- (any(grepl("Rcmdr", search())) &&
              any(c("doItAndPrint", "justDoIt") %in% sysCallNames))
    first <- isTRUE(all.equal(par("mfg")[1:2], c(1,1)))
    opar <- par(no.readonly=TRUE)
    if (dev.interactive() && !playing && first && !par("new")) {
        ## starting a new plot on an interactive device
        frameNum <- which(sysCallNames != "")[1]
        ## find plot command sent from Rcmdr
        if (rcmdr && any(sysCallNames == "eval")) {
            frameNum <- tail(which(sysCallNames == "eval"), 1) + 1
        }
        if (any(StateEnv$.autoplay.ask)) {
            items <- make.unique(sapply(sys.calls(), function(x)
                                        toString(deparseOneLine(x), width=34)))
            frameNum <- menu(items, title="Choose plot call for playwith:",
                             graphics = TRUE)
            if (frameNum == 0) {
                #.C(do_interrupt)
                return()
            }
        }
        parentFrame <- sys.frame(sys.parents()[frameNum])
        newCall <- call("playwith", sys.call(frameNum),
                        envir = parentFrame)
                                        # eval.args=FALSE ?
        eval(newCall)
        par(opar) ## on the new device (redrawing now)
    }
    return()
}

## not to be called directly
playwith.grid.newpage <- function(...)
{
    sysCallNames <- sapply(sys.calls(), function(x)
                           ifelse(is.symbol(x[[1]]), toString(x[[1]]), ""))
    playing <- any(c("playReplot", "playNewPlot", "plot.trellis", "print.trellis")
                   %in% sysCallNames)
    if (dev.interactive() && !playing) {
        ## starting a new plot on an interactive device
        frameNum <- which(sysCallNames != "")[1]
        if (any(StateEnv$.autoplay.ask)) {
            items <- make.unique(sapply(sys.calls(), function(x)
                                        toString(deparseOneLine(x), width=34)))
            frameNum <- menu(items, title="Choose plot call for playwith:",
                             graphics = TRUE)
            if (frameNum == 0)
                return()
        }
        dev.off() ## close screen device (from `grid.newpage`)
        parentFrame <- sys.frame(sys.parents()[frameNum])
        newCall <- call("playwith", sys.call(frameNum),
			envir=parentFrame)
                                        # eval.args=FALSE ?
        eval(newCall)
        ## on the new device (redrawing now)
        ##grid.newpage() ## recursive nightmares
    }
    return()
}
