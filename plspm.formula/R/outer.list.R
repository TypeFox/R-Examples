outer.list <-
function (res.def, res.inter, Data, ovect) 
{
    varsm <- names(Data)
    VLendo <- res.inter[[1]]
    mlist <- res.def[[2]]
    names(mlist) <- res.def[[1]]
    olist <- mlist[ovect]
    outer.l <- list()
    calc.kf <- function(kf) {
        fvars <- olist[[kf]]
        findex <- which(varsm %in% fvars)
        outer.l <<- c(outer.l, list(findex))
    }
    calc.kf <- Vectorize(calc.kf)
    calc.kf(names(olist))
    return(outer.l)
}
