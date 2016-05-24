`print.intervals` <-
function (x, len = 6, d = 2, exclude.intercept = TRUE, pval = TRUE, 
    ...) 
{
    n <- x
    dd <- dim(n)
    mx <- 10^(len - (d + 1))
    n[n > mx] <- Inf
    a <- formatC(n, d, len, format = "f")
    dim(a) <- dd
    if (length(dd) == 1) {
        dd <- c(1, dd)
        dim(a) <- dd
        lab <- " "
    }
    else lab <- dimnames(n)[[1]]
    if (!pval) {
        mx <- max(nchar(lab)) + 1
        cat(paste(rep(" ", mx), collapse = ""), paste(" ", dimnames(n)[[2]]), 
            "\n")
        for (i in (1 + exclude.intercept):dd[1]) {
            lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]), 
                collapse = "")
            if (i == (1 + exclude.intercept)) 
                cat(lab[i], formatC(n[i, 1], 4, 6, format = "f"), 
                  a[i, 2], "Reference haplotype", "\n")
            else cat(lab[i], ifelse(is.na(n[i, 1]),"      ",formatC(n[i, 1], 4, 6, format = "f")), 
                a[i, 2], "(", a[i, 3], "-", a[i, 4], ") \n")
        }
    }
    else {
        mx <- max(nchar(lab)) + 1
        cat(paste(rep(" ", mx), collapse = ""), paste(" ", dimnames(n)[[2]]), 
            "\n")
        for (i in (1 + exclude.intercept):dd[1]) {
            lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]), 
                collapse = "")
            if (i == (1 + exclude.intercept)) 
                cat(lab[i], formatC(n[i, 1], 4, 6, format = "f"), 
                  a[i, 2], "Reference haplotype", "\n")
            else cat(lab[i], ifelse(is.na(n[i, 1]),"      ",formatC(n[i, 1], 4, 6, format = "f")), 
                a[i, 2], "(", a[i, 3], "-", a[i, 4], ") ", formatC(n[i, 
                  5], 4, 6, format = "f"), "\n")
        }
    }
}

