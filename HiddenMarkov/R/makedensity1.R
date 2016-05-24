makedensity1 <- function(distn){
    dname <- paste("d", distn, sep="")
    header <- "function(x, pm, pn=NULL, log=FALSE, prt=FALSE){"
    arglist <- "c(list(x=x), pm, pn, list(log=log))"
    prt <- paste("    if (prt) print(", arglist, ")", sep="")
    fcall <- paste("    do.call(\"", dname, "\", ", arglist, ")",
                   sep="")
    finish <- "}"
    func <- paste(header, prt, fcall, finish, sep="\n")
    eval(parse(text=func))
}

