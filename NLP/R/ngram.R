ngrams <-
function(x, n)
{
    N <- length(x)
    n <- n[(n >= 1L) & (n <= N)]
    lapply(unlist(lapply(n,
                         function(k) {
                             pos <- seq_len(k)
                             lapply(seq.int(0, N - k),
                                    `+`, pos)
                         }),
                  recursive = FALSE),
           function(e) x[e])
}
