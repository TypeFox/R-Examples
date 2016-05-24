"SeasonalMannKendall" <-
function(x)
{
    if(!is.ts(x))
        stop("error: input must be ts object")
    fr <- frequency(x)
    Score <- 0
    Denom <- 0
    varScore <- 0
    for(i in 1:fr) {
        ans <- MannKendall(x[cycle(x) == i])
        Score <- Score + ans$S
        Denom <- Denom + ans$D
        varScore <- varScore + ans$varS
    }
    sl <- 2 * (1 - pnorm(abs(Score/sqrt(varScore))))
    ans <- list(tau = Score/Denom, sl = sl, S=Score, D=Denom, varS=varScore)
    oldClass(ans) <- "Kendall"
    ans
}

