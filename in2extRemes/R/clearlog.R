clearlog <- function(base.txt) {

    out <- paste("Cleared log file: ", date(), sep="")
    print(out)
    write(out, file="in2extRemes.log", append=FALSE)
    invisible()

}
