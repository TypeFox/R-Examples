H600 <-
function () 
{
    Rcmdr <- options()$Rcmdr
    output.height <- if (!is.null(Rcmdr$suppress.menus) &&
                         Rcmdr$suppress.menus) 18 else 15
    H600.options <- list(log.font.size = 9, log.width = 80,
                         log.height = 5, output.height = output.height,
                         scale.factor = 1.4,
                         placement="-0+0")
    Rcmdr[names(H600.options)] <- H600.options
    options(Rcmdr = Rcmdr)
    trellis.par.set(list(superpose.symbol = list(pch = rep(16, 
        length(trellis.par.get("superpose.symbol")$pch))), plot.symbol = list(pch = 16)))
    par(pch = 16)
    putRcmdr("autoRestart", TRUE)
    closeCommander(ask = FALSE)
    Commander()
}
## H600()
## source("~/HH-R.package/RcmdrPlugin.HH/R/H600.R")
