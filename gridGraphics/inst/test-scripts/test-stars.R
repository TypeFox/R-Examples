
library(gridGraphics)

require(grDevices)

stars1 <- function() {
    stars(mtcars[, 1:7], key.loc = c(14, 2),
          main = "Motor Trend Cars : stars(*, full = F)", full = FALSE)
}

stars2 <- function() {
    stars(mtcars[, 1:7], key.loc = c(14, 1.5),
          main = "Motor Trend Cars : full stars()", flip.labels = FALSE)
}

stars3 <- function() {
    ## 'Spider' or 'Radar' plot:
    stars(mtcars[, 1:7], locations = c(0, 0), radius = FALSE,
          key.loc = c(0, 0), main = "Motor Trend Cars", lty = 2)
}

stars4 <- function() {
    ## Segment Diagrams:
    palette(rainbow(12, s = 0.6, v = 0.75))
    stars(mtcars[, 1:7], len = 0.8, key.loc = c(12, 1.5),
          main = "Motor Trend Cars", draw.segments = TRUE)
}

stars5 <- function() {
    stars(mtcars[, 1:7], len = 0.6, key.loc = c(1.5, 0),
          main = "Motor Trend Cars", draw.segments = TRUE,
          frame.plot = TRUE, nrow = 4, cex = .7)
}

## scale linearly (not affinely) to [0, 1]
USJudge <- apply(USJudgeRatings, 2, function(x) x/max(x))
Jnam <- row.names(USJudgeRatings)
Snam <- abbreviate(substring(Jnam, 1, regexpr("[,.]",Jnam) - 1), 7)

stars6 <- function() {
    stars(USJudge, labels = Jnam, scale = FALSE,
          key.loc = c(13, 1.5), main = "Judge not ...", len = 0.8)
}

stars7 <- function() {
    stars(USJudge, labels = Snam, scale = FALSE,
          key.loc = c(13, 1.5), radius = FALSE)
}

stars8 <- function() {
    loc <- stars(USJudge, labels = NULL, scale = FALSE,
                 radius = FALSE, frame.plot = TRUE,
                 key.loc = c(13, 1.5), main = "Judge not ...", len = 1.2)
    text(loc, Snam, col = "blue", cex = 0.8, xpd = TRUE)
}

stars9 <- function() {
    ## 'Segments':
    stars(USJudge, draw.segments = TRUE, scale = FALSE, key.loc = c(13,1.5))
}

stars10 <- function() {
    ## 'Spider':
    stars(USJudgeRatings, locations = c(0, 0), scale = FALSE, radius  =  FALSE,
          col.stars = 1:10, key.loc = c(0, 0), main = "US Judges rated")
}

stars11 <- function() {
    ## Same as above, but with colored lines instead of filled polygons.
    stars(USJudgeRatings, locations = c(0, 0), scale = FALSE, radius  =  FALSE,
          col.lines = 1:10, key.loc = c(0, 0), main = "US Judges rated")
}

stars12 <- function() {
    ## 'Radar-Segments'
    stars(USJudgeRatings[1:10,], locations = 0:1, scale = FALSE,
          draw.segments = TRUE, col.segments = 0, col.stars = 1:10,
          key.loc =  0:1,
          main = "US Judges 1-10 ")
}

stars13 <- function() {
    stars(cbind(1:16, 10*(16:1)), draw.segments = TRUE,
          main = "A Joke -- do *not* use symbols on 2D data!")
}

plotdiff(expression(stars1()), "stars-1")
plotdiff(expression(stars2()), "stars-2")
plotdiff(expression(stars3()), "stars-3")
plotdiff(expression(stars4()), "stars-4")
plotdiff(expression(stars5()), "stars-5")
plotdiff(expression(stars6()), "stars-6")
plotdiff(expression(stars7()), "stars-7")
plotdiff(expression(stars8()), "stars-8")
plotdiff(expression(stars9()), "stars-9")
plotdiff(expression(stars10()), "stars-10")
plotdiff(expression(stars11()), "stars-11")
plotdiff(expression(stars12()), "stars-12")
plotdiff(expression(stars13()), "stars-13")

plotdiffResult()
