RPalpha <-
function(RP.out # the result of a call to RPParellelEnsembleOut
        ,Y # Traning set classes if samplesplit = FALSE, Validation set classes if samplesplit = TRUE
        ,p1 # prior prob estimate
         )
   {
    n <- length(Y)
    Train.Class <- RP.out[1:n, ]
    vote1 <- rowMeans(Train.Class[Y == 1, ], na.rm = TRUE)
    vote2 <- rowMeans(Train.Class[Y == 2, ], na.rm = TRUE)
    errecdfm <- function(x) {
        p1 * ecdf(vote1)(x) + (1 - p1) * (1 - ecdf(vote2)(x))
    }
    errecdfM <- function(x) {
        p1 * ecdf(vote1)(-x) + (1 - p1) * (1 - ecdf(vote2)(-x))
    }
    alpham <- optimise(errecdfm, c(1, 2), maximum = TRUE)$maximum
    alphaM <- optimise(errecdfM, c(-2, -1), maximum = TRUE)$maximum
    alpha <- (alpham - alphaM)/2
    return(alpha)
}
