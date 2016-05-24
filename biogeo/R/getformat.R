getformat <-
function (x) 
{
    n <- length(x)
    fm <- rep("_", n)
    for (j in 1:n) {
        a <- as.character(x[j])
        aa <- {
        }
        for (i in 1:nchar(a)) {
            a1 <- substr(a, i, i)
            suppressWarnings(b <- as.numeric(substr(a, i, i)))
            b1 <- ifelse(is.na(b), "*", "0")
            if (b1 == "*") {
                b1 <- ifelse(a1 == "N" | a1 == "S" | a1 == "W" | 
                  a1 == "E", "L", "*")
                b1 <- ifelse(a1 == ".", ".", b1)
            }
            aa <- paste(aa, b1, sep = "")
        }
        fm[j] <- aa
    }
    return(fm)
}
