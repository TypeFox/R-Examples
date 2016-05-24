hypoexpTvaru <- function(process) {
    function(prob) {
        quant <- hypoexpVaru(process)(prob)
        return(rev(cumsum(rev(int.multi(f     = hypoexpRuinprob(process)[['psi']],
                                        nodes = c(quant, Inf))))) / prob + quant)
    }
}
