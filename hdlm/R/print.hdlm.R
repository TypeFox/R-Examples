print.hdlm <-
function(x, ...){
    digits = max(3, getOption("digits") - 3)
    z <- x
    cat("\nCall:\n", paste(deparse(z$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("Parameters:\n", " Observations = ", z$rank[1], ", Variables = ", z$rank[2],"\n", sep="") 
    cat("\n\n")
    invisible(z)
}

