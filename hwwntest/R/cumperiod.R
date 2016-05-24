cumperiod <-
function (x)
{
    data = x
    N = length(data)
    FT = abs(fft(data)/sqrt(N))^2
    q <- N/2 + 1
    periodogram = FT[2:q]

    answer <- cumsum(periodogram)/sum(periodogram)

    wp <- (2:q)/(2*q)

    return(list(wp=wp, cumperiod=answer))
}
