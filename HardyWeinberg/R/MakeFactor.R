MakeFactor <- function (x, coding = c(0, 1, 2))
{
    if (is.factor(x)) {
        x <- as.numeric(x)
        coding <- c(1, 2, 3)
    }
    if (!all(x[!is.na(x)] %in% coding))
        stop("Error: incorrect coding of the genotype data.\n There exist genotypes not specified in the coding.")
    labs <- names(table(x))
    if (length(labs) == 3)
        xf <- factor(x, levels = coding, labels = c("AA", "AB",
            "BB"))
    if (length(labs) == 2) {
        if (all(labs == c(as.character(coding[1]), as.character(coding[2]))))
            xf <- factor(x, levels = c(coding[1], coding[2]),
                labels = c("AA", "AB"))
        if (all(labs == c(as.character(coding[1]), as.character(coding[3]))))
            xf <- factor(x, levels = c(coding[1], coding[3]),
                labels = c("AA", "BB"))
        if (all(labs == c(as.character(coding[2]), as.character(coding[3]))))
            xf <- factor(x, levels = c(coding[2], coding[3]),
                labels = c("AB", "BB"))
    }
    if (length(labs) == 1) {
        if (labs == as.character(coding[1])) xf <- factor(x, levels = coding[1], labels = "AA")
        if (labs == as.character(coding[2])) xf <- factor(x, levels = coding[2], labels = "AB")
        if (labs == as.character(coding[3])) xf <- factor(x, levels = coding[3], labels = "BB")
    }
    return(xf)
}

