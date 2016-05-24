diffproc <-
function (text) 
{
    if (length(text) != 4) 
        stop("object of length 4 is required")
    if (!inherits(text, "list")) 
        text <- as.list(text)
    errors <- !(unlist(lapply(text, is.character)))
    n.errors <- sum(errors)
    if (n.errors > 0) 
        stop(paste("element", rep("s", n.errors > 1), " [", paste((1:4)[errors], 
            collapse = ","), "] in object ", rep(c("is", "are"), 
            c(n.errors == 1, n.errors > 1)), " not character", 
            sep = ""))
    errors <- unlist(lapply(text, function(l) inherits(try(parse(text = l), 
        silent = TRUE), "try-error")))
    n.errors <- sum(errors)
    if (n.errors > 0) 
        stop(paste("the mathematical expression of element", 
            rep("s", n.errors > 1), " [", paste((1:4)[errors], 
                collapse = ","), "] in object show", rep("s", 
                n.errors < 2), " syntax errors", sep = ""))
    errors <- unlist(lapply(text[2:3], function(l) inherits(try(D(parse(text = l), 
        "x"), silent = TRUE), "try-error")))
    n.errors <- sum(errors)
    if (n.errors > 0) 
        stop(paste("R can not compute the symbolic derivative with respect to 'x' of the mathematical expression of element", 
            rep("s", n.errors > 1), " [", paste((2:3)[errors], 
                collapse = ","), "] in object. \n  (use another equivalent expression", 
            rep("s", n.errors > 1), " from simple functions in table of derivatives)", 
            sep = ""))
    out <- lapply(text, function(l) as.character(parse(text = l)))
    attr(out, "names") <- c("mean", "var", "tpdf", "tpdF")
    class(out) <- c("diffproc", "list")
    return(out)
}
