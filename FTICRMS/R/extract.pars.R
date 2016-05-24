`extract.pars` <-
function(par.file = "parameters.RData", root.dir="."){
    load(paste(root.dir, "/", par.file, sep=""))
    ret <- list()
    for(i in setdiff(ls(),c("i","ret"))){
        ret[[i]] <- get(i)
    }
    ret
}

