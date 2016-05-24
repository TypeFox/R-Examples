# last modified 2012-09-17 by J. Fox

ram <- function(object, digits=getOption("digits"), startvalues=FALSE){
    var.names <- rownames(object$A)
    ram <- object$ram
    if (!startvalues) colnames(ram) <- c(colnames(ram)[1:4], "estimate")
    par <- object$coeff
    par.names <- rep(" ", nrow(ram))
    t <- object$t
    for (i in 1:t) {
        which.par <- ram[,4] == i
        if (!startvalues)  ram[which.par, 5] <- par[i]
        par.names[which.par] <- names(par)[i]
        }
    par.code <- paste(var.names[ram[,2]], c('<---', '<-->')[ram[,1]],
                    var.names[ram[,3]])
    ram <- data.frame(ram, arrow = par.code)
    print(ram, rowlab=par.names, digits=digits)
    }
