`correl` <-
function (x,keep.rares.correl=FALSE){

    boa.init()
    chain.import(x,keep.rares.chain=keep.rares.correl)
    work <- boa.chain("work")
    boa.par()

    cat("\n", "CROSS-CORRELATION MATRIX:\n", "=========================\n", 
        sep = "")
    for (i in names(work)) {
       
        corr <- round(cor(work[[i]]), digits = options()$digits)
        lt <- lower.tri(corr, diag = TRUE)
        corr[!lt] <- ""
        print(corr, quote = FALSE)
        cat("\nPress <ENTER> to continue")
        readline()
    }
   
}

