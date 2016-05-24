`print.snpOut` <-
function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    temp <- attr(x, "label.snp")
    temp2 <- gsub("(snp\\()", "", temp)
    attr(x, "label.snp") <- gsub("\\)", "", temp2)
    cat("\n")
    if (!attr(x, "Interaction")) {
        if (is.null(attr(x, "strata"))) {
            cat("SNP:", attr(x, "label.snp"), " adjusted by:", 
                attr(x, "varAdj"), "\n")
            class(x) <- NULL
            attr(x, "varAdj") <- attr(x, "label.snp") <- attr(x, 
                "BigTable") <- attr(x, "Interaction") <- NULL
            if (!is.Monomorphic(x)) {
                print(x, na.print = "", digits = digits, quote = FALSE)
            }
            else cat("Monomorphic\n")
        }
        else {
            nstrat <- length(attr(x, "strata"))
            for (i in 1:nstrat) {
                cat("       strata:", attr(x, "strata")[i], "\n")
                cat("SNP:", attr(x, "label.snp"), " adjusted by:", 
                  attr(x, "varAdj"), "\n")
                class(x) <- NULL
                attr(x, "varAdj") <- attr(x, "label.snp") <- attr(x, 
                  "BigTable") <- attr(x, "Interaction") <- NULL
                if (!is.Monomorphic(x)) {
                  print(x[[i]], na.print = "", digits = digits, 
                    quote = FALSE)
                }
                else cat("Monomorphic\n")
                cat("\n")
            }
        }
    }
    else {
        cat("      SNP:", attr(x, "label.snp"), " adjusted by:", 
            attr(x, "varAdj"), "\n")
        cat(" Interaction \n")
        cat("---------------------\n")
        print(x[[1]], digits = digits)
        cat("\n")
        cat("p interaction:", x[[4]], "\n")
        cat("\n", paste(attr(x, "label.int"), "within",attr(x, "label.snp")), "\n")
        cat("---------------------\n")
        for (i in 1:length(x[[2]])) {
            cat(names(x[[2]])[i], "\n")
            print(x[[2]][[i]], digits = digits)
            cat("\n")
        }
        cat("p trend:", x[[5]], "\n")
        cat("\n", paste(attr(x, "label.snp"),"within",attr(x, "label.int")), "\n")
        cat("---------------------\n")
        for (i in 1:length(x[[3]])) {
            cat(names(x[[3]])[i], "\n")
            print(x[[3]][[i]], digits = digits)
            cat("\n")
        }
        cat("p trend:", x[[6]], "\n")
    }
}

