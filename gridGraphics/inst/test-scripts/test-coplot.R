
library(gridGraphics)

coplot1 <- function() {
    ## Tonga Trench Earthquakes
    coplot(lat ~ long | depth, data = quakes)
}

coplot2 <- function() {
    given.depth <- co.intervals(quakes$depth, number = 4, overlap = .1)
    coplot(lat ~ long | depth, data = quakes, given.v = given.depth, rows = 1)
}

ll.dm <- lat ~ long | depth * mag

coplot3 <- function() {
    ## Conditioning on 2 variables:
    coplot(ll.dm, data = quakes)
}

coplot4 <- function() {
    coplot(ll.dm, data = quakes, number = c(4, 7), show.given = c(TRUE, FALSE))
}

coplot5 <- function() {
    coplot(ll.dm, data = quakes, number = c(3, 7),
           overlap = c(-.5, .1)) # negative overlap DROPS values
}

Index <- seq(length = nrow(warpbreaks)) # to get nicer default labels

coplot6 <- function() {
    ## given two factors
    coplot(breaks ~ Index | wool * tension, data = warpbreaks,
           show.given = 0:1)
}

coplot7 <- function() {
    coplot(breaks ~ Index | wool * tension, data = warpbreaks,
           col = "red", bg = "pink", pch = 21,
           bar.bg = c(fac = "light blue"))
}

coplot8 <- function() {
    ## Example with empty panels:
    with(data.frame(state.x77), {
        coplot(Life.Exp ~ Income | Illiteracy * state.region, number = 3,
               panel = function(x, y, ...) panel.smooth(x, y, span = .8, ...))
    })
}

coplot9 <- function() {
    ## y ~ factor -- not really sensible, but 'show off':
    with(data.frame(state.x77), {
        coplot(Life.Exp ~ state.region | Income * state.division,
               panel = panel.smooth)
    })
}

plotdiff(expression(coplot1()), "coplot-1")
plotdiff(expression(coplot2()), "coplot-2")
plotdiff(expression(coplot3()), "coplot-3", height=14)
plotdiff(expression(coplot4()), "coplot-4", height=10)
plotdiff(expression(coplot5()), "coplot-5", height=15)
plotdiff(expression(coplot6()), "coplot-6")
plotdiff(expression(coplot7()), "coplot-7")
plotdiff(expression(coplot8()), "coplot-8", height=10)
plotdiff(expression(coplot9()), "coplot-9", width=10, height=15)

plotdiffResult()

