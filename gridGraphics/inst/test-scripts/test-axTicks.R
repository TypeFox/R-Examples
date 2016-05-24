
library(gridGraphics)

axTicks1 <- function() {
    plot(1:7, 10*21:27)
    axTicks(1)
    axTicks(2)
}

axTicks2 <- function() {
    ## Show how axTicks() and axis() correspond :
    par(mfrow = c(3, 1), mar = c(4, 4, 1, 1))
    for(x in 9999 * c(1, 2, 8)) {
        plot(x, 9, log = "x")
        cat(formatC(par("xaxp"), width = 5),";", T <- axTicks(1),"\n")
        rug(T, col =  adjustcolor("red", 0.5), lwd = 4)
    }
}

axTicks3 <- function() {
    x <- 9.9*10^(-3:10)
    plot(x, 1:14, log = "x")
    axTicks(1) # now length 5, in R <= 2.13.x gave the following
    axTicks(1, nintLog = Inf) # rather too many
}

axTicks4 <- function() {
    ## An example using axTicks() without reference to an existing plot
    ## (copying R's internal procedures for setting axis ranges etc.),
    ## You do need to supply _all_ of axp, usr, log, nintLog
    ## standard logarithmic y axis labels
    ylims <- c(0.2, 88)
    get_axp <- function(x) 10^c(ceiling(x[1]), floor(x[2]))
    ## mimic par("yaxs") == "i"
    usr.i <- log10(ylims)
    (aT.i <- axTicks(side = 2, usr = usr.i,
                     axp = c(get_axp(usr.i), n = 3), log = TRUE, nintLog = 5))
    ## mimic (default) par("yaxs") == "r"
    usr.r <- extendrange(r = log10(ylims), f = 0.04)
    (aT.r <- axTicks(side = 2, usr = usr.r,
                     axp = c(get_axp(usr.r), 3), log = TRUE, nintLog = 5))

    ## Prove that we got it right :
    plot(0:1, ylims, log = "y", yaxs = "i")
}

axTicks5 <- function() {
    ylims <- c(0.2, 88)
    get_axp <- function(x) 10^c(ceiling(x[1]), floor(x[2]))
    ## mimic par("yaxs") == "i"
    usr.i <- log10(ylims)
    (aT.i <- axTicks(side = 2, usr = usr.i,
                     axp = c(get_axp(usr.i), n = 3), log = TRUE, nintLog = 5))
    ## mimic (default) par("yaxs") == "r"
    usr.r <- extendrange(r = log10(ylims), f = 0.04)
    (aT.r <- axTicks(side = 2, usr = usr.r,
                     axp = c(get_axp(usr.r), 3), log = TRUE, nintLog = 5))
    plot(0:1, ylims, log = "y", yaxs = "r")
}

plotdiff(expression(axTicks1()), "axTicks-1")
plotdiff(expression(axTicks2()), "axTicks-2", width=8, height=8)
plotdiff(expression(axTicks3()), "axTicks-3")
plotdiff(expression(axTicks4()), "axTicks-4")
plotdiff(expression(axTicks5()), "axTicks-5")

plotdiffResult()
