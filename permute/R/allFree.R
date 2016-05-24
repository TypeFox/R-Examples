## Modified version of allFree() provided by Doug Bates
## via personal email on 19 Jan 2012
`allFree` <- function(n, v = seq_len(n)) {
    if(n == 1L) return(array(v, c(1L, 1L)))
    do.call(rbind,
            lapply(seq_len(n),
                   function(i) cbind(v[i], allFree(n - 1L, v[-i]))))
}
